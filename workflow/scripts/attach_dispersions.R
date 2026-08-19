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

message(sprintf("attached: %d genes, %d with a fitted dispersion (these become the pair list)",
                nrow(simulated_sce), sum(!is.na(params$dispersion))))

saveRDS(simulated_sce, file = snakemake@output$simulated_sce_disp)

sink(); sink(type = "message"); close(log)
