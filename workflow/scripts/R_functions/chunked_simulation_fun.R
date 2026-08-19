## Memory-bounded (gene-chunked) versions of the power-simulation count generator.
##
## WHY: simulate_tapseq_counts() materialises THREE dense n_genes x n_cells
## doubles per rep (`mu`, the effect-size matrix, and the rnbinom draw). For the
## NPC scaffold (5,009 genes x 46,000 cells) that is ~1.84 GB each, ~5.5 GB peak
## per split, and it scales linearly with the gene count -- a full transcriptome
## would need ~40 GB per split.
##
## The generation is independent across genes: gene i's counts depend only on
## gene i's mean/dispersion, the per-cell size factors, and gene i's own
## effect-size draws. So we generate in gene blocks and accumulate into ONE
## sparse matrix. Peak dense allocation becomes O(chunk_size x n_cells) and is
## independent of the total gene count.
##
## RESULT EQUIVALENCE: statistically identical, not bit-identical -- chunking
## reorders the RNG stream. Note the pipeline has no set.seed() anywhere in the
## simulation path, so consecutive runs of the ORIGINAL code already differ by
## Monte Carlo noise; this does not add a new class of variability.
##
## The full-height sparse result keeps every row the caller expects, so the
## sceptre_object stays valid and run_discovery_analysis() is still called once
## on the complete matrix -- no sceptre internals are touched.

# Effect-size matrix for a block of genes. Mirrors create_effect_size_matrix()
# exactly, but only for `gene_effect_sizes_chunk`.
create_effect_size_matrix_chunk <- function(grna_pert_status, pert_guides,
                                            gene_effect_sizes_chunk, guide_sd) {

  n_pert_guides <- length(pert_guides)
  n_ctrl_guides <- max(grna_pert_status) - n_pert_guides

  guide_effect_sizes_pert <- vapply(gene_effect_sizes_chunk, FUN = rnorm,
                                    n = n_pert_guides, sd = guide_sd,
                                    FUN.VALUE = numeric(n_pert_guides))
  guide_effect_sizes_ctrl <- vapply(rlang::rep_along(gene_effect_sizes_chunk, 1),
                                    FUN = rnorm, n = n_ctrl_guides, sd = guide_sd,
                                    FUN.VALUE = numeric(n_ctrl_guides))
  guide_effect_sizes <- rbind(guide_effect_sizes_pert, guide_effect_sizes_ctrl)

  # set negative guide effect sizes to 0
  guide_effect_sizes[guide_effect_sizes < 0] <- 0

  # add row with no effect for non-perturbed cells
  guide_effect_sizes <- rbind(1, guide_effect_sizes)

  # pick correct effect sizes for every cell based on its gRNA perturbation status
  es_mat <- t(guide_effect_sizes[grna_pert_status + 1, , drop = FALSE])
  colnames(es_mat) <- names(grna_pert_status)

  es_mat
}

# Centering, applied within a gene block. center_effect_size_matrix() works
# row-wise (per gene), so restricting to a block is exact.
center_effect_size_matrix_chunk <- function(effect_size_mat, pert_status,
                                            gene_effect_sizes_chunk) {
  mean_es_pert <- rowMeans(effect_size_mat[, pert_status == 1, drop = FALSE])
  mean_es_ctrl <- rowMeans(effect_size_mat[, pert_status == 0, drop = FALSE])

  pert_shift <- gene_effect_sizes_chunk - mean_es_pert
  ctrl_shift <- 1 - mean_es_ctrl

  effect_size_mat[, pert_status == 1] <- effect_size_mat[, pert_status == 1] + pert_shift
  effect_size_mat[, pert_status == 0] <- effect_size_mat[, pert_status == 0] + ctrl_shift

  effect_size_mat[effect_size_mat < 0] <- 0
  effect_size_mat
}

#' Simulate the response matrix in gene blocks.
#'
#' Drop-in replacement for the
#'   create_effect_size_matrix -> center_effect_size_matrix -> sim_tapseq_sce
#'   -> as(..., "RsparseMatrix")
#' sequence in sceptre_power_analysis.R.
#'
#' @return dgRMatrix, genes x cells, same dimnames as the baseline path.
simulate_response_matrix_chunked <- function(pert_object, grna_pert_status,
                                             pert_guides, effect_sizes,
                                             guide_sd, chunk_size = 500L) {

  cell_ids  <- colnames(assay(pert_object, "counts"))
  pert_status <- colData(pert_object)$pert

  gene_means <- rowData(pert_object)[, "mean"]
  gene_disps <- rowData(pert_object)[, "dispersion"]
  size_factors <- colData(pert_object)[, "size_factors"]
  gene_ids <- rownames(pert_object)

  n_genes <- length(gene_ids)
  n_cells <- length(size_factors)
  chunk_size <- max(1L, as.integer(chunk_size))
  starts <- seq(1L, n_genes, by = chunk_size)

  message("Simulating counts in ", length(starts), " gene chunk(s) of <= ",
          chunk_size, " genes (", n_genes, " genes x ", n_cells, " cells)")

  blocks <- vector("list", length(starts))

  for (bi in seq_along(starts)) {
    idx <- starts[bi]:min(starts[bi] + chunk_size - 1L, n_genes)

    # --- effect sizes for this gene block ---
    es_chunk <- create_effect_size_matrix_chunk(
      grna_pert_status, pert_guides = pert_guides,
      gene_effect_sizes_chunk = effect_sizes[idx], guide_sd = guide_sd
    )
    es_chunk <- center_effect_size_matrix_chunk(
      es_chunk, pert_status = pert_status,
      gene_effect_sizes_chunk = effect_sizes[idx]
    )
    # align cells to the counts assay, exactly as the baseline does
    es_chunk <- es_chunk[, cell_ids, drop = FALSE]

    # --- mu = gene_mean x cell_size_factor x effect_size ---
    mu <- matrix(rep(gene_means[idx], n_cells), ncol = n_cells)
    mu <- sweep(mu, 2, size_factors, "*")
    mu <- mu * es_chunk
    rm(es_chunk)

    # --- negative-binomial draw for this block ---
    counts_chunk <- matrix(
      rnbinom(length(mu), mu = mu, size = 1 / gene_disps[idx]),
      ncol = n_cells
    )
    rm(mu)

    dimnames(counts_chunk) <- list(gene_ids[idx], cell_ids)
    # collapse to sparse immediately -- this is where the dense memory is released
    blocks[[bi]] <- as(as(counts_chunk, "Matrix"), "CsparseMatrix")
    rm(counts_chunk)
    gc(verbose = FALSE)
  }

  out <- do.call(rbind, blocks)
  rm(blocks); gc(verbose = FALSE)

  as(out, "RsparseMatrix")
}
