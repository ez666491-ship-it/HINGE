# HINGE

Code accompanying **“Direct identification of genetically central traits in shared genetic architecture from GWAS summary statistics.”**

## Files

- `simulation.R`: data generation, simulation comparisons, and F1/TPR/precision plots.
- `real_data.R`: GWAS summary-data analysis, LD clumping, and combined trait-score/spectral plots.
- `data/`: example trait files named `example1.tsv`, `example2.tsv`, etc. Example files are not yet included.

## Running the code

Use the repository root as the working directory. All paths are relative.

Install the required R packages once:

```r
install.packages(c(
  "MASS", "glasso", "data.table", "R.utils", "ieugwasr",
  "tidyverse", "RColorBrewer", "cowplot", "ggh4x",
  "ggrepel", "patchwork"
))
```

Run either script:

```r
source("simulation.R")
source("real_data.R")
```

Each script performs analysis followed by plotting. Set `run_analysis <- FALSE` at the top to plot existing results only; set `run_plots <- FALSE` to run analysis only. Results and figures are saved under `results/`. Existing output files are overwritten.

## Simulation parameters

- `p`: number of traits (100, 200).
- `n`: GWAS sample size (30,000, 50,000, 100,000).
- `m`: number of SNPs (500, 1,000, 2,000).
- `T0_prop`: proportion of traits in the hub-enriched block (0.5, 0.75, 1).
- `ph`: hub-edge probability (0.4, 0.8).
- `UHP`: absence/presence of pleiotropic contamination (0, 1).
- `r`: number of true hubs (5).
- `pnh` / `pneff`: non-hub/background edge probabilities (0.05 / 0.01).
- `diagonal_shift`: minimum population precision eigenvalue (2).
- Edge weights: random signs and absolute magnitudes drawn from Uniform(4, 5).
- `UHPratio` / `UHP_n_traits` / `UHPshift`: contaminated SNP fraction (0.1), trait count (20), and shift coefficient (1).
- `w_mult`: additional error-vector count as a multiple of m (10).
- `ipchd_perc`: off-diagonal correlation threshold quantile (0.7).
- `egg33_subtime`: EGG subsampling repetitions in the full run (15).
- `R` / `seed0`: simulation repetitions (100) and seed base (123).

The complete grid contains 216 scenarios. It is computationally intensive. The execution section defaults to one CPU core; `HINGE_CORES` can increase this on Linux/macOS. Windows uses one core. The supplied plots require the complete grid and do not include the separate overall sorted F1-distribution figure.

## Real-data inputs

- `data_dir`: input directory, default `data/`.
- One trait per `exampleN.tsv` file; trait count is detected automatically, with no manual list. The current spectral rule requires at least six traits.
- Required columns: `variant`, `beta`, `se`, `pval`, `minor_AF`, `low_confidence_variant`.
- `variant`: chromosome:position:allele1:allele2; SNP identifiers, genome build and effect-allele orientation must agree across files.
- `se` must be positive; `pval` is the ordinary association P value.
- Random 500-row subsets may not contain enough jointly significant SNPs to complete analysis.

### PLINK and reference panel

Users download these separately; neither is uploaded with this repository:

- [PLINK 1.9](https://www.cog-genomics.org/plink/).
- [1000 Genomes LD reference and local LD instructions](https://mrcieu.r-universe.dev/ieugwasr/doc/local_ld.html), including the [reference archive](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz).

Place the executable at `tools/plink.exe` (Windows) or `tools/plink` (Linux/macOS), and the selected population's files at `reference/ld_reference.bed`, `.bim`, and `.fam`. Alternatively edit the relative paths near the start of `real_data.R`. Select a reference population and genome build appropriate to the data.

SNP IDs are matched through the reference BIM file; no separate mapping file is required. Strand conversion and liftover are not performed.

## Implementation

- Real-data input uses Z statistics (beta/SE).
- Null SNPs satisfy P > 0.05 in every trait; error covariance averages 300 estimates from 10% subsamples.
- Signal SNPs satisfy joint-test P < 5e-8 and undergo PLINK LD clumping (r2 = 0.001; 100 kb). Null SNPs are not LD-clumped.
- Both workflows use corrected correlation matrices. Simulation additionally thresholds off-diagonal entries; real data do not.
- Real-data eigenvalues are floored at 1e-6; simulation eigenvalues below 0.001 are set to zero.
- Hub classification uses mean score + 2 SD; scores describe panel-specific centrality, not causality.
- Internal output labels `raw`, `raw_deg`, `ip_thr`, `egg_theta_alpha`, `egg_theta_deg` correspond to Raw_alpha, Raw_deg, HINGE, EGG_alpha and EGG_deg.

The manuscript reports R 4.3.2. Scripts save session information. This consolidated package has not yet been execution-tested or plot-rendered in the packaging environment.
