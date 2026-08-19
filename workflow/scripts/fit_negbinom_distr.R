## Fit negative binomial distributions for each gene using DESeq2

# Saving image for debugging
save.image("RDA_objects/fit_dispersions.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DESeq2)
  # chunked_deseq2_fun.R is self-contained (Matrix + DESeq2 only) and provides
  # everything this rule needs. power_simulations_fun.R is sourced only for the
  # ORIGINAL fit_negbinom_deseq2(), which the chunked path replaces -- and it
  # opens with library(tidyverse), which the DESeq2 conda env does not carry.
  # So source it opportunistically rather than making it a hard dependency.
  source(file.path(snakemake@scriptdir, "R_functions/chunked_deseq2_fun.R"))
  if (requireNamespace("tidyverse", quietly = TRUE)) {
    source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
  } else {
    message("tidyverse absent; skipping power_simulations_fun.R (chunked path does not need it)")
  }
})

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
response_matrix <- readRDS(snakemake@input$raw_counts)
colnames(response_matrix) <- NULL
simulated_sce <- readRDS(snakemake@input$simulated_sce)

# Calculate total_umis and detected_genes for Deseq2 object creation
coldata <- data.frame(
  total_umis = colSums(response_matrix),
  detected_genes = colSums(response_matrix > 0)
)

### CHOOSE THE TWO GENE SETS ================================================
# FIX C decomposes estimateDispersions() into its three stages, which lets the
# two expensive passes be restricted independently:
#
#   trend_genes -- stage 1, whose only product is the cross-gene trend plus
#                  varLogDispEsts/dispPriorVar. A 2-parameter curve and two
#                  scalars need an UNBIASED sample, not every gene.
#   map_genes   -- stage 3, needed only for genes that will actually be
#                  simulated and tested. Exact for the genes retained; genes
#                  left out get dispersion NA, which is already how
#                  split_target_response_pairs.R marks a response untestable.
#
# Measured 2026-08-19 (validate_levers.R): map_genes is BIT-IDENTICAL,
# trend_genes costs 0.0003 median relative dispersion error. A third lever --
# subsampling CELLS for the fit -- was REJECTED: 0.098 median dispersion error
# and baseMean shifts up to 1.94, which feed straight into rnbinom().
n_trend <- as.integer(snakemake@params$n_trend_genes)
n_map   <- as.integer(snakemake@params$n_map_genes)

set.seed(as.integer(snakemake@params$gene_selection_seed))
trend_genes <- if (n_trend > 0 && n_trend < nrow(response_matrix)) {
  sort(sample(nrow(response_matrix), n_trend))
} else NULL

map_genes <- if (n_map > 0 && n_map < nrow(response_matrix)) {
  sel <- stratified_expression_genes(response_matrix, n_genes = n_map,
                                     seed = as.integer(snakemake@params$gene_selection_seed))
  # force the screen's own target genes in -- they are the whole point of the
  # run, and stratified sampling would otherwise include them only by luck
  force_file <- snakemake@params$force_genes_file
  if (!is.null(force_file) && nzchar(force_file) && file.exists(force_file)) {
    forced <- match(readLines(force_file), rownames(simulated_sce))
    forced <- forced[!is.na(forced)]
    message("forcing ", length(forced), " target genes into map_genes")
    sel <- sort(unique(c(sel, forced)))
  }
  sel
} else NULL

message(sprintf("trend_genes = %s, map_genes = %s (of %d)",
                if (is.null(trend_genes)) "ALL" else length(trend_genes),
                if (is.null(map_genes)) "ALL" else length(map_genes),
                nrow(response_matrix)))

# fit negative binomial distributions to estimate gene-level dispersion
message("Estimate dispersion using DESeq2 (chunked):")
sce <- fit_negbinom_deseq2_chunked(response_matrix,
                                   simulated_sce,
                                   coldata,
                                   size_factors = "poscounts",
                                   fit_type = "parametric",
                                   disp_type = "dispersion",
                                   gene_chunk = as.integer(snakemake@params$gene_chunk),
                                   trend_genes = trend_genes,
                                   map_genes = map_genes,
                                   n_cores = as.integer(snakemake@threads))


### SAVE PLAIN PARAMETERS, NOT AN SCE =======================================
# This rule runs in the DESeq2 conda env, which is a newer R/Bioconductor than
# the sceptre env that consumes the result. Writing a SingleCellExperiment here
# bakes in that env's S4 class definitions (e.g. the Seqinfo package, split out
# of GenomeInfoDb in recent Bioc), and the older env then cannot deserialize it:
#   "unable to find required package 'Seqinfo'"
# Reading old -> new is fine; writing new -> old is not. So emit plain atomic
# vectors, which carry no class baggage, and let attach_dispersions.R (running
# in the sceptre env) fold them back onto the SCE that env created itself.
#
# The dispersionFunction is deliberately NOT carried across: it is a closure
# over the DESeq2 namespace, nothing downstream reads it (it is only ever
# written), and serialising it would reintroduce the same coupling.
saveRDS(list(
  gene_ids            = rownames(response_matrix),
  mean                = rowData(sce)[, "mean"],
  dispersion          = rowData(sce)[, "dispersion"],
  disp_outlier_deseq2 = rowData(sce)[, "disp_outlier_deseq2"],
  size_factors        = colData(sce)[, "size_factors"]
), file = snakemake@output$dispersion_params)

# close log file connection
sink()
sink(type = "message")
close(log)





