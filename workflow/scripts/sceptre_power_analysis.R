# Script to create a sceptre object based off some simulated counts for downstream power simulations

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
# save.image("RDA_objects/sceptre_power_analysis.rda")
# message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES ============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(sceptre)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/chunked_simulation_fun.R"))
})

# Download Sceptre (skip if already installed -- avoids N concurrent jobs each
# racing a from-source recompile of the same package into the same conda env)
library(devtools)
if (!requireNamespace("sceptre", quietly = TRUE)) {
  message("Installing Sceptre")
  devtools::install_github("katsevich-lab/sceptre")
  message("Sceptre Installation Complete")
} else {
  message("Sceptre already installed (", as.character(utils::packageVersion("sceptre")), "); skipping install_github")
}

# Load inputs
simulated_sceptre_object <- readRDS(snakemake@input$simulated_sceptre_object)
simulated_sce_disp <- readRDS(snakemake@input$simulated_sce_disp)
discovery_pairs_split <- read_tsv(snakemake@input$discovery_pairs_split, col_names = c("grna_group", "response_id"))
grna_target_data_frame <- read_tsv(snakemake@input$grna_target_data_frame)
response_matrix <- readRDS(snakemake@input$raw_counts)

# Reduce raw_counts to the one thing this script actually needs from it -- the
# per-gene mean, for the average_expression_all_cells column -- and release it
# immediately. It is otherwise held for the whole rep loop in EVERY split: at a
# full transcriptome that is 6.7 GB x 10 concurrent splits = ~67 GB of duplicated
# memory, for a single rowMeans() call.
#
# Indexing BY NAME also keeps this correct when simulated_sce_disp has been
# subset to the tested genes while raw_counts still spans the transcriptome --
# the two no longer have the same number of rows, and positional use would
# silently mismatch (or, as it did, error out after the reps have all run).
gene_means_all <- rowMeans(response_matrix)
gene_means <- gene_means_all[rownames(simulated_sce_disp)]
stopifnot(!anyNA(gene_means))
rm(response_matrix, gene_means_all); gc(verbose = FALSE)


# Load params
effect_size <- 1 - as.numeric(snakemake@params$effect_size)
reps <- snakemake@params$reps
guide_sd <- 0.13
# genes per block for the chunked count generator (Fix B). Dense working set is
# gene_chunk_size x n_cells doubles, x3 concurrent -- 1,000 x 46,000 ~= 368 MB
# each. Independent of the total gene count.
gene_chunk_size <- as.integer(snakemake@params$gene_chunk_size)


### PRECOMPUTATIONS TO RUN SIMULATIONS ======================================

# Get all the perts to be run
perts <- unique(discovery_pairs_split$grna_group)

# Create an empty results data frame
discovery_results <- data.frame()


### RUN POWER SIMULATIONS ===================================================

for (pert in perts) {
 
  # Initialize the pert object with the given pert and sce
  pert_object <- pert_input(pert, sce = simulated_sce_disp, pert_level = "cre_perts")
  
  # Get all the guides that target the current `pert`
  pert_guides <- grna_target_data_frame %>%
    filter(grna_target == pert) %>%
    pull(grna_id)
  
  # Get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  # Convert to a sparse matrix, so the sampling function works in `create_guide_pert_status`
  grna_perts <- as(grna_perts, "CsparseMatrix")
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts, pert_guides = pert_guides)
  
  # Create effect size matrix (sampled from negative binomial distribution around effect_size or 1)
  effect_sizes <- structure(rep(effect_size, nrow(pert_object)), names = rownames(pert_object))
  
  # Loop through each rep
  for (rep in seq(reps)) {
    
    # FIX B: generate the counts in gene blocks rather than materialising three
    # dense n_genes x n_cells doubles (mu, the effect-size matrix, and the
    # rnbinom draw). Peak becomes O(chunk_size x n_cells), independent of the
    # gene count. Statistically identical, not bit-identical -- chunking
    # reorders the RNG stream, but the simulation path has no set.seed() at all,
    # so two runs of the ORIGINAL code already differ by the same amount.
    # Validated 2026-08-19: a 3-arm test (baseline x2 for the noise floor, vs
    # chunked) put the chunked arm exactly on the noise floor.
    message("Simulating Counts (gene-chunked)")
    simulated_response_matrix <- simulate_response_matrix_chunked(
      pert_object, grna_pert_status = grna_pert_status, pert_guides = pert_guides,
      effect_sizes = effect_sizes, guide_sd = guide_sd, chunk_size = gene_chunk_size
    )


    # Save a temp sceptre object
    sceptre_object_use <- simulated_sceptre_object
    
    # Assign the simulated response matrix to the sceptre object
    sceptre_object_use@response_matrix <- list(simulated_response_matrix)
    # Subset and assign the discovery pairs relevant to the current perturbation to the sceptre object
    sceptre_object_use@discovery_pairs <- discovery_pairs_split[discovery_pairs_split$grna_group == pert,]
    
    # Run the discovery analysis
    message("Running discovery analysis")
    sceptre_object_use <- run_discovery_analysis(
      sceptre_object = sceptre_object_use,
      parallel = FALSE
    )

    # Get the discovery analysis results
    message("Returning discovery results")
    discovery_result <- get_result(
      sceptre_object = sceptre_object_use,
      analysis = "run_discovery_analysis"
    )

    # Add the number of perturbed cells and the rep to each pair
    n_pert_cells <- length(sceptre_object_use@grna_assignments$grna_group_idxs[[pert]])
    discovery_result$num_pert_cells <- n_pert_cells
    discovery_result$rep <- rep
    
    # Save the results
    discovery_results <- data.frame(rbind(discovery_results, discovery_result))

  }
}


### RUN POST COMPUTATIONS ===================================================

message("Processing output.")

# Construct a combined data frame
combined_data <- data.frame(
  response_id = rownames(simulated_sce_disp),
  disp_outlier_deseq2 = rowData(simulated_sce_disp)[, "disp_outlier_deseq2"],
  dispersion = rowData(simulated_sce_disp)[, "dispersion"],
  average_expression_all_cells = gene_means,
  stringsAsFactors = FALSE
)

# Merge the combined data with discovery_results
discovery_results <- left_join(discovery_results, combined_data, by = "response_id")


### SAVE OUTPUT =============================================================

# save simulation output
message("Saving output to file.")
write_tsv(discovery_results, file = snakemake@output[[1]])


# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)