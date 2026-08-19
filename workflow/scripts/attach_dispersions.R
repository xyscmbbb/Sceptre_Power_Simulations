## Fold the DESeq2 parameters back onto the SingleCellExperiment.
##
## Exists purely to keep S4 objects inside one conda env. fit_negbinom_distr.R
## runs under the DESeq2 env (newer R/Bioconductor) and emits plain vectors;
## this script runs under the sceptre env, which created simulated_sce.rds in
## the first place and can therefore both read and write it safely.

message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(SingleCellExperiment))

simulated_sce <- readRDS(snakemake@input$simulated_sce)
params <- readRDS(snakemake@input$dispersion_params)

stopifnot(identical(params$gene_ids, rownames(simulated_sce)))
stopifnot(length(params$size_factors) == ncol(simulated_sce))

rowData(simulated_sce)[, "mean"] <- params$mean
rowData(simulated_sce)[, "dispersion"] <- params$dispersion
rowData(simulated_sce)[, "disp_outlier_deseq2"] <- params$disp_outlier_deseq2
colData(simulated_sce)[, "size_factors"] <- params$size_factors

### SUBSET TO THE TESTED GENES ==============================================
# NA dispersion already excludes a gene from the PAIR LIST
# (split_target_response_pairs.R drops is.na(dispersion) responses), but it does
# NOT exclude it from the SIMULATION: create_simulated_sceptre_object.R and
# pert_input() both walk every row of this object, so an untested gene is still
# drawn for, at full cost, via rnbinom(mu = NA, size = 1/NA) -> a column of NAs.
#
# With map_genes restricting the fit, that is the difference between simulating
# 2,021 genes and 34,597 -- ~17x the work per rep, plus a response matrix that is
# mostly NA. So drop those rows here, which makes the simulated gene set and the
# tested gene set the same thing by construction.
#
# Safe because nothing downstream needs the dropped rows: the sceptre covariates
# (response_n_umis / response_n_nonzero) are computed from the FULL raw_counts
# in create_simulated_sceptre_object.R, not from this object, and
# average_expression_all_cells is left_join()ed by response_id, so extra rows in
# raw_counts are simply unused.
keep <- !is.na(params$dispersion)
message(sprintf("attached: %d genes, %d with a fitted dispersion", nrow(simulated_sce), sum(keep)))
if (isTRUE(as.logical(snakemake@params$subset_to_tested)) && any(!keep)) {
  simulated_sce <- simulated_sce[keep, ]
  message(sprintf("subset to %d tested genes (dropped %d without a fitted dispersion)",
                  nrow(simulated_sce), sum(!keep)))
}

saveRDS(simulated_sce, file = snakemake@output$simulated_sce_disp)

sink(); sink(type = "message"); close(log)
