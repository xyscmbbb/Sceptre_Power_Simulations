# Script to simulate guide assignments given parameters in config file

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/simulate_guide_assignments.rda")
message("Saved Image")
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
  source(file.path(snakemake@scriptdir, "R_functions/simulate_perturbations.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
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
library(sceptre)

# Load in input and the response matrix
message("Loading input and response matrix")
response_matrix <- readRDS(snakemake@input$raw_counts)
colnames(response_matrix) <- NULL

# Load in snakemake parameters
message("Loading parameters")
num_cells_per_pert <- snakemake@params$num_cells_per_pert
num_guides_per_pert <- 15


### SIMULATE GUIDE COUNTS ====================================================

message("Creating the sce object")
# Estimate how many cells we need in our sce object
num_cells <- num_cells_per_pert[length(num_cells_per_pert)] * 5 + 6000 # The 6000 is to make sure there are enough n_ctrl cells but this calculation is arbitrary otherwise

# Create the raw matrix with the new number of cells.
# FIX A: sparse, not a dense all-zero base R matrix. This assay is NEVER read
# for its values -- the only consumer is
# colnames(assay(pert_object, "counts")) in sceptre_power_analysis.R -- so it is
# a pure dimnames placeholder. Dense it costs nrow x num_cells x 8 B, which is
# serialised into simulated_sce_disp.rds and re-loaded by EVERY split:
# 13.2 GB per split at a full transcriptome (34,597 x 46,000). Sparse is
# bit-identical for every downstream consumer.
counts <- Matrix::Matrix(0, nrow = nrow(response_matrix), ncol = num_cells, sparse = TRUE)
colnames(counts) <- paste("cell", 1:num_cells, sep="")
rownames(counts) <- rownames(response_matrix)

# Create the sce object
sce <- SingleCellExperiment(assays=list(counts=counts))

# Simulate perturbations
message("Simulating perturbations")
sce <- simulate_perturbations(sce, 
                              cells_per_pert=num_cells_per_pert, 
                              guides_per_pert=num_guides_per_pert)


# Fill in row data names for grna_perts
rowData(altExp(sce, "grna_perts"))$name <- rownames(altExp(sce, "grna_perts"))
rowData(altExp(sce, "grna_perts"))$target_name <- sub("_.*", "", rownames(altExp(sce, "grna_perts")))

# save simulation output
message("Saving output to file.")
saveRDS(sce, file = snakemake@output$simulated_sce)

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)