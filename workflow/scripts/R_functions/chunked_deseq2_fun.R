## Memory-bounded (gene-chunked, sparse-aware) replacement for fit_negbinom_deseq2().
##
## WHY: fit_negbinom_deseq2() calls DESeqDataSetFromMatrix(), which DENSIFIES the
## raw counts. For a full transcriptome that is 36,000 x 160,799 x 8 B = ~46 TB.
## Even the current 5,009-gene NPC matrix costs 6.4 GB dense. The gene axis cannot
## be trimmed instead: the dispersion trend is fit ACROSS genes and the per-cell
## covariates are library-size sums, so both need the full gene set.
##
## WHAT IS EXACT HERE. DESeq2's estimateDispersions() is a three-stage pipeline
##   estimateDispersionsGeneEst -> estimateDispersionsFit -> estimateDispersionsMAP
## in which only the middle stage is cross-gene. Reading the DESeq2 source
## (v1.42-era) the complete list of quantities that couple genes together is:
##
##   1. the dispersion trend, from (baseMean, dispGeneEst) SCALARS only
##      -- estimateDispersionsFit never touches counts();
##   2. varLogDispEsts, an attribute of the dispersion function, used by MAP
##      for the dispOutlier flag;
##   3. dispPriorVar, from estimateDispersionsPriorVar(), which needs
##      varLogDispEsts and m = nrow(modelMatrix) = n_cells;
##   4. maxDisp = max(10, ncol(object)) -- n_cells, identical in every chunk;
##   5. the poscounts size factors, which are a per-cell median over genes.
##
## Everything else -- baseMean/baseVar, mu, dispGeneEst, and the MAP fit itself --
## is per-gene. So we compute (1)-(4) once, globally, from per-gene scalars and
## inject them into each chunk. The result is mathematically identical to the
## one-shot call, not merely statistically equivalent. `validate_chunked_deseq2.R`
## checks that empirically against the one-shot path on a gene subset.
##
## Peak memory becomes O(chunk_size x n_cells) plus the sparse input, independent
## of the total gene count.

## Only Matrix + DESeq2 are needed here. The rowData/colData/metadata accessors
## used below are SummarizedExperiment generics (loaded by DESeq2), so this file
## works on any SummarizedExperiment -- including the pipeline's
## SingleCellExperiment -- without depending on SingleCellExperiment itself.
## That matters locally: the `deseq2` conda env has DESeq2 but not
## SingleCellExperiment, and `potc` has the reverse.
suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
})

#' poscounts size factors, computed on the sparse matrix without densifying.
#'
#' Reproduces exactly, for `type = "poscounts"`:
#'   geoMeanNZ(x) = if (all(x == 0)) 0 else exp(sum(log(x[x > 0])) / length(x))
#'   sf[j]        = exp(median over {i : x[i,j] > 0, loggeomeans[i] finite}
#'                               of (log(x[i,j]) - loggeomeans[i]))
#'   sf           = sf / exp(mean(log(sf)))          # incomingGeoMeans branch
#'
#' Both statistics only ever look at NONZERO entries, so the sparse
#' representation carries all the information needed; the zeros contribute
#' nothing but the divisor `length(x) = ncol`.
#'
#' @param counts sparse genes x cells matrix (any Matrix class; coerced to CSC)
#' @param cell_chunk number of cells to process per block
sparse_poscounts_size_factors <- function(counts, cell_chunk = 20000L) {

  # column access dominates, so work in CSC even if handed a dgRMatrix
  if (!inherits(counts, "CsparseMatrix")) counts <- as(counts, "CsparseMatrix")

  n_genes <- nrow(counts)
  n_cells <- ncol(counts)

  # --- log geometric means over positive counts only -----------------------
  # sum over nonzero entries of log(x), per gene, divided by n_cells.
  lognz <- counts
  lognz@x <- log(lognz@x)
  loggeomeans <- Matrix::rowSums(lognz) / n_cells
  rm(lognz)

  all_zero <- Matrix::rowSums(counts) == 0
  loggeomeans[all_zero] <- -Inf          # geoMeanNZ returns 0 -> log(0) = -Inf
  finite_gm <- is.finite(loggeomeans)

  if (all(!finite_gm)) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }

  # --- per-cell median of (log count - loggeomean) over that cell's nonzeros -
  sf <- numeric(n_cells)
  starts <- seq(1L, n_cells, by = as.integer(cell_chunk))
  for (s in starts) {
    e <- min(s + as.integer(cell_chunk) - 1L, n_cells)
    blk <- counts[, s:e, drop = FALSE]
    p <- blk@p
    i <- blk@i + 1L                      # 0-based -> 1-based gene index
    x <- blk@x
    for (k in seq_len(e - s + 1L)) {
      rng <- seq_len(p[k + 1L] - p[k]) + p[k]
      if (length(rng) == 0L) { sf[s + k - 1L] <- NA_real_; next }
      gi <- i[rng]
      keep <- finite_gm[gi] & x[rng] > 0
      sf[s + k - 1L] <- exp(stats::median((log(x[rng]) - loggeomeans[gi])[keep]))
    }
    rm(blk); gc(verbose = FALSE)
  }

  sf / exp(mean(log(sf)))                # incomingGeoMeans = TRUE
}

#' Chunked, exact replacement for fit_negbinom_deseq2(..., size_factors =
#' "poscounts", fit_type = "parametric", disp_type = "dispersion").
#'
#' Signature and return value match the original: the simulated_sce with
#' rowData mean/dispersion/disp_outlier_deseq2, colData size_factors, and
#' metadata dispersionFunction.
#'
#' @param gene_chunk genes per block. Dense working set is gene_chunk x n_cells
#'   doubles (2,000 x 160,799 x 8 B ~= 2.6 GB; 500 -> ~640 MB).
#' @param trend_genes integer/character gene selector for stage 1, whose ONLY
#'   product is the cross-gene dispersion trend plus varLogDispEsts and
#'   dispPriorVar. Those are a 2-parameter curve fit and two scalars, so they do
#'   not need all ~36,000 genes -- they need an UNBIASED sample spanning the
#'   expression range. (The objection to an HVG-selected matrix is bias, not
#'   sample size: HVGs are high-variance by construction, which drags the trend
#'   up and biases simulated power DOWN.) A random sample of ~10k genes is
#'   statistically equivalent; verify with validate_levers.R before relying on
#'   it. NULL = all genes.
#' @param map_genes integer/character gene selector for stage 3. Only genes that
#'   will actually be SIMULATED AND TESTED need a final MAP dispersion, and the
#'   MAP fit is per-gene given the global trend -- so restricting it is exact for
#'   the genes retained. Genes outside `map_genes` come back with dispersion NA,
#'   which is already how upstream marks a gene as untestable:
#'   split_target_response_pairs.R drops `is.na(dispersion)` responses from the
#'   pair list. NULL = all genes.
fit_negbinom_deseq2_chunked <- function(response_matrix,
                                        simulated_sce,
                                        coldata,
                                        size_factors = "poscounts",
                                        fit_type = "parametric",
                                        disp_type = "dispersion",
                                        gene_chunk = 500L,
                                        cell_chunk = 20000L,
                                        trend_genes = NULL,
                                        map_genes = NULL,
                                        n_cores = 1L,
                                        quiet = FALSE) {

  stopifnot(identical(size_factors, "poscounts"))
  stopifnot(identical(fit_type, "parametric"))

  n_genes <- nrow(response_matrix)
  n_cells <- ncol(response_matrix)
  gene_chunk <- max(1L, as.integer(gene_chunk))

  # resolve the two gene selectors to sorted integer indices
  as_idx <- function(sel) {
    if (is.null(sel)) return(seq_len(n_genes))
    if (is.character(sel)) sel <- match(sel, rownames(response_matrix))
    if (is.logical(sel)) sel <- which(sel)
    sel <- sort(unique(as.integer(sel)))
    stopifnot(!anyNA(sel), all(sel >= 1L), all(sel <= n_genes))
    sel
  }
  trend_idx <- as_idx(trend_genes)
  map_idx   <- as_idx(map_genes)

  blocks_of <- function(idx) split(idx, ceiling(seq_along(idx) / gene_chunk))
  trend_chunks <- blocks_of(trend_idx)
  map_chunks   <- blocks_of(map_idx)

  if (!quiet) message(sprintf(
    "chunked DESeq2: %d genes x %d cells | trend on %d gene(s), MAP on %d gene(s) | chunk <= %d",
    n_genes, n_cells, length(trend_idx), length(map_idx), gene_chunk))

  ## ---- stage 0: size factors (sparse, exact) ------------------------------
  if (!quiet) message("size factors (poscounts, sparse)")
  sf <- sparse_poscounts_size_factors(response_matrix, cell_chunk = cell_chunk)

  # design is ~ 1, so the model matrix is a single intercept column over cells.
  # Passed explicitly everywhere, since our per-chunk objects must agree with
  # the real cell count for dispPriorVar (m - p) and maxDisp.
  model_matrix <- stats::model.matrix(~ 1, data = coldata)

  # Chunks are pure functions of (their genes, sf, model_matrix) plus -- for
  # stage 3 -- the global trend scalars. No RNG is involved in either stage, so
  # running them in parallel is EXACTLY result-preserving, not merely
  # statistically equivalent. Without this the rule is one single-threaded job:
  # ~12,000 gene-passes x ~1.7 s at 160,799 cells is ~5.7 h.
  # Each worker densifies gene_chunk x n_cells (500 x 160,799 x 8 B = 643 MB)
  # and DESeq2 holds a few copies -- budget ~3 GB per core.
  run_chunks <- function(chunks, label, fn) {
    if (n_cores <= 1L || length(chunks) == 1L) {
      out <- vector("list", length(chunks))
      for (ci in seq_along(chunks)) {
        out[[ci]] <- fn(chunks[[ci]])
        if (!quiet) message(sprintf("  %s chunk %d/%d", label, ci, length(chunks)))
      }
      return(out)
    }
    gc(verbose = FALSE)   # shrink the parent before forking
    if (!quiet) message(sprintf("  %s: %d chunks across %d cores",
                                label, length(chunks), n_cores))
    out <- parallel::mclapply(chunks, fn, mc.cores = n_cores, mc.preschedule = FALSE)
    bad <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(bad)) stop(label, ": ", sum(bad), " chunk(s) failed; first: ",
                       conditionMessage(attr(out[[which(bad)[1]]], "condition")))
    out
  }

  build_chunk <- function(idx) {
    sub <- as.matrix(response_matrix[idx, , drop = FALSE])
    mode(sub) <- "integer"
    dds <- DESeqDataSetFromMatrix(countData = sub, colData = coldata, design = ~ 1)
    sizeFactors(dds) <- sf
    dds
  }

  ## ---- stage 1: per-gene dispersion estimates (chunked) -------------------
  if (!quiet) message("gene-wise dispersion estimates")
  base_mean <- numeric(n_genes)
  base_var  <- numeric(n_genes)
  all_zero  <- logical(n_genes)
  disp_gene_est <- numeric(n_genes)

  trend_res <- run_chunks(trend_chunks, "trend", function(idx) {
    dds <- build_chunk(idx)
    dds <- estimateDispersionsGeneEst(dds, modelMatrix = model_matrix, quiet = TRUE)
    m <- mcols(dds)
    out <- list(baseMean = m$baseMean, baseVar = m$baseVar,
                allZero = m$allZero, dispGeneEst = m$dispGeneEst)
    rm(dds, m); gc(verbose = FALSE)
    out
  })
  for (ci in seq_along(trend_chunks)) {
    idx <- trend_chunks[[ci]]; r <- trend_res[[ci]]
    base_mean[idx] <- r$baseMean; base_var[idx] <- r$baseVar
    all_zero[idx]  <- r$allZero;  disp_gene_est[idx] <- r$dispGeneEst
  }
  rm(trend_res); gc(verbose = FALSE)

  ## ---- stage 2: the cross-gene trend, on SCALARS only ---------------------
  # estimateDispersionsFit reads only mcols(baseMean, dispGeneEst, allZero) --
  # it never calls counts(). So a 2-column shell carrying the real per-gene
  # scalars gives a bit-identical dispersion function at negligible cost.
  if (!quiet) message("mean-dispersion relationship")
  # The shell's counts are never read -- every mcols column the fit consumes is
  # overwritten below. Two DIFFERENT dummy values only to dodge DESeqDataSet's
  # "all genes have equal values for all samples" warning, which would be
  # misleading noise in the log.
  shell <- DESeqDataSetFromMatrix(
    countData = matrix(c(1L, 2L), nrow = length(trend_idx), ncol = 2, byrow = TRUE,
                       dimnames = list(rownames(response_matrix)[trend_idx], c("a", "b"))),
    colData = data.frame(row.names = c("a", "b")), design = ~ 1)
  mcols(shell)$baseMean    <- base_mean[trend_idx]
  mcols(shell)$baseVar     <- base_var[trend_idx]
  mcols(shell)$allZero     <- all_zero[trend_idx]
  mcols(shell)$dispGeneEst <- disp_gene_est[trend_idx]

  shell <- estimateDispersionsFit(shell, fitType = fit_type, quiet = quiet)
  disp_function <- dispersionFunction(shell)      # carries varLogDispEsts + fitType

  # dispPriorVar is global: it needs varLogDispEsts and m = n_cells (NOT the
  # shell's 2 columns), hence the explicit modelMatrix.
  disp_prior_var <- estimateDispersionsPriorVar(shell, modelMatrix = model_matrix)
  var_log_disp_ests <- attr(disp_function, "varLogDispEsts")
  if (!quiet) message(sprintf("  dispPriorVar = %.6f  varLogDispEsts = %.6f",
                              disp_prior_var, var_log_disp_ests))
  rm(shell); gc(verbose = FALSE)

  ## ---- stage 3: MAP shrinkage per chunk, with the GLOBAL prior ------------
  if (!quiet) message("final dispersion estimates")
  dispersion  <- rep(NA_real_, n_genes)
  disp_outlier <- rep(NA, n_genes)

  disp_fit <- rep(NA_real_, n_genes)

  map_res <- run_chunks(map_chunks, "MAP", function(idx) {
    dds <- build_chunk(idx)
    # geneEst must be re-run to regenerate the per-gene `mu` assay that MAP
    # consumes; it is per-gene and deterministic, so this reproduces stage 1.
    # Caching mu from stage 1 instead would mean holding an n_genes x n_cells
    # dense assay -- precisely the allocation this whole file exists to avoid.
    # Hence the ~2x compute cost: bounded memory bought with a second pass.
    # (It also means map_genes need NOT be a subset of trend_genes -- their
    # baseMean/dispGeneEst are produced right here.)
    dds <- estimateDispersionsGeneEst(dds, modelMatrix = model_matrix, quiet = TRUE)

    # Inject the global trend. The setter derives dispFit per gene as
    # disp_function(baseMean) -- already exactly what we want -- but ALSO
    # recomputes varLogDispEsts as mad(residuals)^2 over whatever genes happen
    # to be in this chunk. That one is cross-gene, and MAP uses it for the
    # dispOutlier threshold, so it is overwritten with the global value.
    suppressMessages(dispersionFunction(dds) <- disp_function)
    dfn <- dds@dispersionFunction
    attr(dfn, "varLogDispEsts") <- var_log_disp_ests
    dds@dispersionFunction <- dfn

    dds <- estimateDispersionsMAP(dds, dispPriorVar = disp_prior_var,
                                  modelMatrix = model_matrix, quiet = TRUE)
    m <- mcols(dds)
    out <- list(baseMean = m$baseMean, dispGeneEst = m$dispGeneEst,
                dispFit = m$dispFit, dispersion = m$dispersion,
                dispOutlier = m$dispOutlier)
    rm(dds, m, dfn); gc(verbose = FALSE)
    out
  })
  for (ci in seq_along(map_chunks)) {
    idx <- map_chunks[[ci]]; r <- map_res[[ci]]
    base_mean[idx]     <- r$baseMean
    disp_gene_est[idx] <- r$dispGeneEst
    disp_fit[idx]      <- r$dispFit
    dispersion[idx]    <- r$dispersion
    disp_outlier[idx]  <- r$dispOutlier
  }
  rm(map_res); gc(verbose = FALSE)

  ## ---- assemble the SCE exactly as fit_negbinom_deseq2() does -------------
  # Genes outside map_genes keep dispersion NA, which is precisely how upstream
  # marks a response as untestable -- split_target_response_pairs.R drops
  # is.na(dispersion) rows from the pair list. So restricting map_genes both
  # skips the work AND selects the tested gene set, with no extra patch.
  mean_out <- rep(NA_real_, n_genes)
  mean_out[map_idx] <- base_mean[map_idx]
  rowData(simulated_sce)[, "mean"] <- mean_out
  rowData(simulated_sce)[, "dispersion"] <- switch(
    disp_type,
    dispersion  = dispersion,
    dispFit     = disp_fit,
    dispGeneEst = disp_gene_est,
    stop("unsupported disp_type for the chunked path: ", disp_type))
  rowData(simulated_sce)[, "disp_outlier_deseq2"] <- disp_outlier

  # size factors are resampled with replacement into the scaffold's cells, so
  # only their empirical distribution matters (upstream behaviour, unchanged).
  colData(simulated_sce)[, "size_factors"] <- sample(sf, size = ncol(simulated_sce),
                                                     replace = TRUE)
  # match the one-shot path, which stores the function AFTER MAP has stamped
  # dispPriorVar onto it (stored for provenance only; nothing reads it back).
  attr(disp_function, "dispPriorVar") <- disp_prior_var
  metadata(simulated_sce)[["dispersionFunction"]] <- disp_function

  simulated_sce
}

#' Choose the gene set to actually simulate and test (the `map_genes` argument).
#'
#' The deliverable is a power-vs-expression curve with one line per
#' cells-per-perturbation tier, so the gene axis only has to DEFINE the curve,
#' not enumerate the transcriptome. What matters is even coverage of the
#' expression range -- especially the low end, where the 80%-power crossing
#' lives -- rather than gene count.
#'
#' Sampling equal numbers per log-expression bin (rather than at random) is the
#' point: gene abundance is heavily skewed, so a random sample concentrates in
#' the low-expression pile-up and leaves the transition region thin.
#'
#' Genes with rowSum < 2 are dropped -- DESeq2 returns NA dispersion for them,
#' so upstream would exclude them from the pair list anyway.
#'
#' @return integer gene indices into response_matrix.
stratified_expression_genes <- function(response_matrix, n_genes = 2000L,
                                        n_bins = 40L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  row_sums  <- Matrix::rowSums(response_matrix)
  gene_mean <- row_sums / ncol(response_matrix)
  eligible  <- which(row_sums >= 2 & is.finite(gene_mean) & gene_mean > 0)
  if (length(eligible) <= n_genes) return(sort(eligible))

  bins <- cut(log10(gene_mean[eligible]), breaks = as.integer(n_bins),
              labels = FALSE, include.lowest = TRUE)
  per_bin <- ceiling(n_genes / length(unique(bins)))

  picked <- unlist(lapply(split(eligible, bins), function(g) {
    if (length(g) <= per_bin) g else sample(g, per_bin)
  }), use.names = FALSE)

  # sparse bins can leave us short of the target; top up at random from the rest
  if (length(picked) < n_genes) {
    rest <- setdiff(eligible, picked)
    picked <- c(picked, sample(rest, min(n_genes - length(picked), length(rest))))
  }
  sort(unique(picked))
}
