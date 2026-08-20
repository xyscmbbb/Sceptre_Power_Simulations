# Script to create a sceptre object based off some simulated counts for downstream power simulations

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/create_simulated_sceptre_object.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES ============================================================

# Download Sceptre (skip if already installed -- avoids N concurrent jobs each
# racing a from-source recompile of the same package into the same conda env)
if (!requireNamespace("sceptre", quietly = TRUE)) {
  message("Installing Sceptre")
  devtools::install_github("katsevich-lab/sceptre")
  message("Sceptre Installation Complete")
} else {
  message("Sceptre already installed (", as.character(utils::packageVersion("sceptre")), "); skipping install_github")
}

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(sceptre)
  library(tidyverse)
  library(Matrix)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/simulate_perturbations.R"))
})


# Load inputs
simulated_sce_disp <- readRDS(snakemake@input$simulated_sce_disp)
grna_target_data_frame <- read_tsv(snakemake@input$grna_target_data_frame)
discovery_pairs <- read_tsv(snakemake@input$discovery_pairs)
raw_counts <- readRDS(snakemake@input$raw_counts)


### SIMULATE COUNTS ========================================================

# Use an example pert and create the pert_object 
result_matrix <- simulate_sample_perturbation_response(
  pert = unique(discovery_pairs$grna_target)[[1]],
  sce = simulated_sce_disp,
  effect_size = 0.95,
  guide_sd = 0.13,
  grna_target_data_frame = grna_target_data_frame
)


### CREATE SCEPTRE OBJECT ==================================================

# Initialize object
sceptre_object <- import_data(
  response_matrix = result_matrix,
  grna_matrix = assay(altExp(simulated_sce_disp, "grna_perts"), "perts"), 
  grna_target_data_frame = grna_target_data_frame,
  moi = "high"
)


### DROP THE SPURIOUS p_mito COVARIATE =====================================

# `import_data()` derives `response_p_mito` whenever it sees MT- genes, and
# `run_discovery_analysis()` later auto-runs `run_qc()`, whose default
# `p_mito_threshold = 0.2` then deletes every cell above it.
#
# On this pipeline that filter is meaningless AND actively harmful, for two
# reasons:
#
#   1. It is computed on SIMULATED counts, so it measures nothing about the
#      real cells' mitochondrial content.
#   2. The denominator is only the genes being simulated, not the
#      transcriptome. Those genes are a stratified sample across log-expression
#      -- deliberately enriched for low-expression genes -- while MT genes are
#      among the highest expressed in any scRNA-seq dataset. So MT survives the
#      subsetting at full strength while the denominator collapses, and the
#      ratio inflates by roughly the reciprocal of the retained UMI fraction.
#
# Measured on NPC (2026-08-20): real MT fraction over all 34,597 genes = 7.2%,
# but 19.8% within the 2,021 simulated genes (which hold 36% of all UMIs).
# Median simulated cell = 0.196 against a 0.2 threshold, so 45% of cells were
# silently discarded -- every cells-per-perturbation tier came out at ~53% of
# its configured size and the top tier (8,000) collapsed into the one below it.
# Job B, whose 25 target genes contain no MT gene, was untouched: the covariate
# was never created, so the same run_qc() call removed only 4% of cells.
#
# Removing the covariate is what makes the two comparable. Cellwise QC on the
# REAL data belongs upstream, where p_mito is measured against the whole
# transcriptome; it must not be re-derived from a simulated gene subset.
if ("response_p_mito" %in% names(sceptre_object@covariate_data_frame)) {
  pm <- sceptre_object@covariate_data_frame$response_p_mito
  n_mt <- length(grep("^MT-", rownames(result_matrix), ignore.case = TRUE))
  message(sprintf(
    paste0("Dropping simulated response_p_mito (%d MT- genes among %d simulated ",
           "genes; median %.3f, %.1f%% of cells above the 0.2 default threshold ",
           "that run_qc() would otherwise have removed)"),
    n_mt, nrow(result_matrix), median(pm), 100 * mean(pm > 0.2)))
  sceptre_object@covariate_data_frame$response_p_mito <- NULL
}


### RESAMPLE THE COVARIATES ================================================

# Calculate the number of response and the number of umis for each cell
response_n_nonzero <- colSums(raw_counts > 0)
response_n_umis <- colSums(raw_counts)

# Sample the covariates with replacement from the true covariates
# Covariate data frame
sceptre_object@covariate_data_frame$response_n_nonzero <- sample(response_n_nonzero, size = ncol(result_matrix), replace = TRUE)
sceptre_object@covariate_data_frame$response_n_umis <- sample(response_n_umis, size = ncol(result_matrix), replace = TRUE)

### CORRECT ERROR IN COVARIATES ============================================

# Because we have ~1 gRNA per column, we have to add 1 to grna_n_nonzero and grna_n_umis (because the log is taken)
sceptre_object@covariate_data_frame$grna_n_nonzero <- sceptre_object@covariate_data_frame$grna_n_nonzero + 1
sceptre_object@covariate_data_frame$grna_n_umis <- sceptre_object@covariate_data_frame$grna_n_umis + 1

# We have to set the formulat to not include grna_n_nonzero
# grna_n_nonzero and grna_n_umis are not linearly independent, so they can't both be used as covariates
formula <- formula(~ log(response_n_nonzero) + log(response_n_umis))


### OTHER SCEPTRE PARAMETERS ===============================================

# Set the analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  side = "left",
  grna_integration_strategy = "union",
  control_group = "complement",
  formula_object = formula,
  resampling_mechanism = "permutations", 
  multiple_testing_method = "BH",
  multiple_testing_alpha = 0.05,
  resampling_approximation = "skew_normal"
)

# Assign gRNAs
sceptre_object <- assign_grnas(
  sceptre_object = sceptre_object,
  method = "thresholding",
  threshold = 1
)


### SAVE OUTPUT =============================================================

# save simulation output
message("Saving output to file.")
saveRDS(sceptre_object, snakemake@output$simulated_sceptre_object)

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
