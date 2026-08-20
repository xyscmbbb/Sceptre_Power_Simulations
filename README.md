# Sceptre-based Power Simulations

> ### About this fork
>
> This is a **lighter** fork of
> [jamesgalante/Sceptre_Power_Simulations](https://github.com/jamesgalante/Sceptre_Power_Simulations).
>
> Lighter, not smaller: it does the same work and answers the same question, but
> its peak memory no longer scales with the number of genes, so it will run a
> **full transcriptome** on an ordinary VM instead of needing a
> subsetted counts matrix. Nothing was removed or approximated to get there —
> every change is either exactly result-preserving or validated to sit on the
> Monte-Carlo noise floor (see *What changed*, below). Upstream's
> documentation from *About The Project* onward still applies unchanged.
>
> Reference run: **160,799 cells x 34,597 genes** (NPC), all genes retained.
>
> ---
>
> ### What the workflow actually does
>
> The question is: *how many cells per perturbation do I need before I can
> detect a given knockdown?* The pipeline answers it by simulation, using two
> different tools for the two halves of the problem:
>
> - a **DESeq2-based count simulator** (from the TAP-seq power analysis) as the
>   *generator* — it learns each gene's mean and dispersion from your real
>   counts, then draws new negative-binomial counts with a knockdown injected
>   into the perturbed cells;
> - **sceptre** as the *tester* — the same method you would use on the real
>   screen, so the estimated power reflects the analysis you will actually run.
>
> Power for a gene is then just the fraction of simulation replicates in which
> sceptre calls that gene significant.
>
> The steps, in order:
>
> | Step | What it does |
> |---|---|
> | `simulate_guide_assignments` | Builds a synthetic cell population with one perturbation per `num_cells_per_pert` tier, so a single run sweeps the whole cells-per-perturbation axis at once |
> | `fit_dispersions` / `attach_dispersions` | Fits the DESeq2 mean-dispersion trend on the real counts — the expensive step, and the only place the real data is used quantitatively |
> | `split_target_response_pairs` | Shards the target x gene pair list so the simulation can run in parallel |
> | `create_simulated_sceptre_object` | Builds the sceptre object once: guide assignments, covariates, analysis parameters |
> | `sceptre_power_analysis` | The rep loop — simulate counts with the effect injected, run sceptre's discovery analysis, record per-pair significance |
> | `compute_power_from_simulations` | Collapses reps into a power per (perturbation, gene) pair |
> | `visualize_power_results` | The report, including the canonical power-vs-gene-expression plot |
>
> **Two things worth knowing before you read any output:**
>
> 1. `effect_size` is the **knockdown fraction, not the residual**. `0.5` means
>    expression is halved (median log2FC ~ -1), not reduced to 50% of a 50%
>    knockdown. Larger `effect_size` = stronger perturbation = more power.
> 2. Power is driven far more by **gene expression** than by cell count. Cells
>    mostly buy you access to *lower-expressed* genes; a well-expressed gene is
>    detectable with very few. Read the power-vs-expression plot, not just the
>    aggregate curve — an aggregate over a full transcriptome may never reach
>    80% at any cell count, simply because the low-expression tail cannot be
>    reached at all.
>
> ---
>
> ### What changed
>
> **Memory, so a full transcriptome fits.** The generator used to materialise
> several dense `n_genes x n_cells` doubles at once, which is what made a full
> transcriptome impossible. Now:
>
> - the all-zero scaffold matrix is **sparse** rather than dense;
> - counts are generated in **gene blocks**, so the dense working set is
>   `chunk_size x n_cells` and is independent of the total gene count;
> - the DESeq2 dispersion fit is likewise **chunked**, with sparse
>   `poscounts` size factors;
> - `raw_counts` is reduced to the per-gene means it is actually needed for and
>   released, instead of being held for the whole rep loop in every parallel
>   split;
> - the SCE is subset to the genes being tested, so simulation cost matches the
>   pair list rather than the input matrix.
>
> **Correctness.** Two bugs that silently changed results:
>
> - `response_p_mito` is dropped before QC. sceptre derives it whenever MT-
>   genes are present, and the implicit `run_qc()` then deletes every cell above
>   `p_mito_threshold = 0.2`. On simulated counts that ratio is meaningless, and
>   its denominator is only the *simulated* genes rather than the transcriptome —
>   so a stratified gene set inflates it several-fold. Measured on NPC: a real
>   7.2% mitochondrial fraction read as 19.8%, straddling the threshold, and
>   **45% of cells were silently discarded** — every cell tier realised at ~53%
>   of its configured size. A targets-only run, whose genes contain no MT gene,
>   was untouched, which is how the divergence surfaced.
> - Post-computation rows are matched **by gene name**, not position, which
>   matters once the SCE and `raw_counts` no longer have the same row count.
>
> **Validation.** Chunking reorders the RNG stream, so results are statistically
> identical rather than bit-identical — and the original code has no `set.seed()`
> in the simulation path, so two runs of *upstream* already differ by the same
> amount. A three-arm test (baseline twice, to measure the noise floor, versus
> chunked) put the chunked arm exactly on that floor:
>
> | | baseline vs baseline | baseline vs chunked |
> |---|---|---|
> | mean power delta | +0.0038 | +0.0042 |
> | per-gene mean abs. difference | 0.0395 | 0.0414 |
> | power correlation | 0.9539 | 0.9518 |
> | KS D (pooled p-values) | 0.0290 (p=0.31) | 0.0290 (p=0.31) |
>
> The chunked DESeq2 fit is exact, not merely equivalent: size factors agree to
> 2.8e-14, baseMean to 2.1e-14, dispersions to 2.2e-09 (inside the solver's own
> `dispTol`), with zero dispersion-outlier disagreements.
>
> ---

## About The Project

The **Sceptre-based Power Simulations** pipeline is designed to help determine the number of cells needed in CRISPRi or CRISPRko experiments to achieve a desired statistical power (e.g., 80%) for detecting a specific percentage decrease in gene expression. By simulating different effect sizes and cell numbers, this tool provides valuable insights for experimental design and helps determine the feasibility of detecting small changes in gene expression.

This pipeline leverages the **sceptre** R package, which offers a statistically rigorous and scalable approach for single-cell CRISPR screen data analysis. It calculates mean expression and dispersion from your input data, simulates a specified percentage decrease in gene expression, and computes the statistical power to detect that change across varying numbers of cells.

## Installation

To get started with the Sceptre-based Power Simulations pipeline, follow these steps:

1. **Clone the Repository**

   ```bash
   git clone https://github.com/jamesgalante/Sceptre_Power_Simulations.git
   cd Sceptre_Power_Simulations
   ```

2. **Create a Conda Environment with Snakemake**

   It's recommended to create a new conda environment specifically for this project. We'll install Snakemake version `7.32.4`.

   ```bash
   conda create -n sceptre_power_sim python=3.9
   conda activate sceptre_power_sim
   conda install -c bioconda snakemake=7.32.4
   ```

   *Note:* All other required packages will be dynamically downloaded when the Snakemake pipeline is run.

## Usage

### Input Data Preparation

The pipeline requires a raw counts matrix in the form of a `.rds` file (created with saveRDS). This file should be a `dgRMatrix` format. An example of the expected format can be found in `resources/test_data/raw_counts.rds`.

If you wish to visualize the power simulations results by TPM instead of UMIs, provide a TPM file formatted as in `resources/tpm_per_gene.tsv`. The TPM file is only used to match your input genes with the power simulation results during plotting. If you are designing a screen based on TPM values, this might be a good option. Since single-cell UMIs may differ from TPM calculated values in bulk screens, the resulting plots may look more jagged.

### Configuring the Pipeline

Before running the pipeline, you need to set up the configuration file `config.yaml`. This file controls various parameters of the simulation and analysis.

Here is an example of the `config.yaml` file:

```yaml
# Main Config File Format for Sceptre-based Power Simulations Pipeline

samples:
  your_sample_name

simulate_guide_assignments:
  num_cells_per_pert: [50, 75, 100, 150, 250, 500, 750, 1000, 1400, 1800, 2500, 4000, 7500, 10000]

sceptre_power_analysis:
  effect_size: 0.15
  reps: 20

compute_power_from_simulations:
  pval_adj_thresh: 0.1
  positive_proportion: 0.05

visualize_power_results:
  tpm_per_gene: "resources/tpm_per_gene.tsv"  # Set to NULL if not using TPM
  gene_format: "ensembl"  # Must be "ensembl" or "symbol"
```

#### Configuration Parameters Explained

- **samples**: A list containing the names of the samples you wish to analyze. Each sample name should correspond to a directory in `resources/` that contains the `raw_counts.rds` file. Replace `your_sample_name` with your actual sample name.

- **simulate_guide_assignments**:
  - **num_cells_per_pert**: A list of numbers representing the different counts of cells per perturbation you want to simulate. The pipeline will assess the statistical power for each specified cell count.

- **sceptre_power_analysis**:
  - **effect_size**: The percent decrease in gene expression you wish to simulate. For example, `0.15` represents a 15% decrease.
  - **reps**: The number of simulation replicates to run for each condition. More replicates reduce noise but increase computational time. A value of `20` is typically sufficient.

- **compute_power_from_simulations**:
  - **pval_adj_thresh**: The adjusted p-value threshold for determining statistical significance after applying Benjamini-Hochberg (BH) correction. The default Sceptre recommendation is `0.1`.
  - **positive_proportion**: The expected proportion of true positives (genes with actual expression changes) in your data. This parameter is important because the BH correction assumes a certain proportion of null hypotheses. This number is hard to estimate before running an experiment, but for CRISPR screens, the proportion of positive tests is typically low. A value of `0.05` corresponds to 5% of the resulting tests being positives.

- **visualize_power_results**:
  - **tpm_per_gene**: Path to the TPM per gene file. If you wish to use UMI counts instead, set this to `NULL`. See the example `resources/tpm_per_gene.tsv` for formatting requirements.
  - **gene_format**: Specifies the format of gene names used in your data. Must be either `"ensembl"` or `"symbol"` or refer to a specific column in the TPM file.

### Running the Pipeline

To execute the pipeline, run the following commands from the root directory of the project:

1. **Dry Run**

   It's good practice to perform a dry run to ensure that the pipeline is configured correctly and all files are in place.

   ```bash
   snakemake --use-conda all -np
   ```

   This command will show you what steps the pipeline will perform without actually executing them.

2. **Full Run**

   If the dry run looks correct, execute the pipeline:

   ```bash
   snakemake --use-conda all
   ```

   *Note:* The `--use-conda` flag tells Snakemake to create and use the necessary conda environments for each rule.
   *Note:* See the relevant documentation [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/index.html) to understand what flags might be necessary. If using slurm, it might be easiest to create a snakemake profile in `~/.config/snakemake/profile_name/config.yml`. Here is an example profile:
   
   ```yaml
   jobs: 100 
   slurm: True 
   retries: 1
   use-conda: True 
   notemp: True 
   default-resources: 
       - slurm_account=slurm_account_name
       - slurm_partition=slurm_account_name,owners,normal
       - runtime="6h"
       - slurm_extra="--nice"
   ```
   The pipeline can then be run with `snakemake --profile profile_name all`

#### Troubleshooting

- If you encounter issues, you can test the pipeline with the provided test data. Rename or delete the `results/test_data/` directory and run the pipeline again to ensure everything works as expected.

- Ensure that your `config.yaml` file is properly configured and that the `raw_counts.rds` object matches the format of the `test_data`.


## How It Works

The pipeline simulates gene expression data to assess the statistical power of detecting specified effect sizes across varying numbers of cells. It begins by calculating mean expression and dispersion for each gene from your input data. Cells are then assigned to control or perturbation groups, simulating a decrease in expression based on the specified effect size. Using the sceptre package, differential expression analysis is performed, and p-values are adjusted using the Benjamini-Hochberg procedure, considering the expected proportion of true positives. Finally, the pipeline computes the power for each condition and generates visualizations to help interpret the results.

## Notes

- **Data Format**: Ensure that your input data is correctly formatted. The gene names should match between your raw counts matrix and the TPM file (if used).

- **Sceptre Package**: The sceptre package is maintained by another team. This pipeline utilizes it for its efficiency and statistical robustness in single-cell CRISPR screen analysis. For more details, refer to their documentation.

- **Understanding FDR Correction**: The Benjamini-Hochberg procedure assumes a certain proportion of true positives. The `positive_proportion` parameter allows you to adjust this assumption based on your experimental context.

- **TO DO**: Use dispersion estimates calculated by Sceptre

## References

- **Sceptre Package Documentation**: [Hands-On Single-Cell CRISPR Screen Analysis](https://timothy-barry.github.io/sceptre-book/)

---
