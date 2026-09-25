# HINGE real-data analysis and plotting. Run from the repository root.
# Set run_analysis <- FALSE to redraw figures from existing results.
run_analysis <- TRUE
run_plots <- TRUE

# ---- Analysis: automatic trait discovery and normal LD clumping ----
if (run_analysis) {
suppressPackageStartupMessages({
  library(data.table)
})

# Real-data workflow: traits are discovered automatically from input files.

# Only the input directory is user-supplied; all default paths are relative.
data_dir <- Sys.getenv("HINGE_DATA_DIR", "data")
out_dir <- "results/real_data"
# External LD resources have conventional relative locations; no config file.
plink_bin <- if (.Platform$OS.type == "windows") "tools/plink.exe" else "tools/plink"
ld_bfile <- "reference/ld_reference"
clump_kb <- 100
clump_r2 <- 0.001
clump_p <- 1
if (!dir.exists(data_dir)) stop("Input directory not found: ", data_dir)
if (!file.exists(plink_bin)) stop("Place PLINK at: ", plink_bin)
if (!all(file.exists(paste0(ld_bfile, c(".bed", ".bim", ".fam"))))) stop("Missing reference/ld_reference.bed/.bim/.fam")
if (!requireNamespace("ieugwasr", quietly = TRUE)) stop("Install ieugwasr before LD clumping.")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

seed_global <- 123

# ---- stage 1 candidate SNP filter ----
maf_thresh <- 0.05
remove_mhc <- TRUE

# ---- stage 2 / EGG-style screening ----
null_min_pval     <- 0.05   # null SNP must have p > 0.05 for all traits
joint_p_threshold <- 5e-8

# ---- Sigma_omega estimation ----
subsample_prop <- 0.10
n_rep_omega    <- 300

# ---- hub threshold ----
hub_sd_multiplier <- 2

set.seed(seed_global)

# =========================================================
# 1. Trait files are discovered below; no manual trait list is required.
# =========================================================

# =========================================================
# 2. Helpers
# =========================================================
log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf(...))
  cat("\n")
  flush.console()
}

mad_sd <- function(x) {
  1.483 * median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
}

spearmancov <- function(A) {
  p <- ncol(A)
  s <- apply(A, 2, mad_sd)
  R_rho <- suppressWarnings(cor(A, method = "spearman", use = "pairwise.complete.obs"))
  R_gauss <- 2 * sin(pi * R_rho / 6)
  D <- diag(s, p, p)
  S <- D %*% R_gauss %*% D
  (S + t(S)) / 2
}

make_spd <- function(M, eps = 1e-6) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  vals <- eig$values
  vals[vals < eps] <- eps
  M_spd <- eig$vectors %*% diag(vals, length(vals)) %*% t(eig$vectors)
  (M_spd + t(M_spd)) / 2
}

cov2cor_safe <- function(S) {
  trait_names <- colnames(S)
  d <- sqrt(diag(S))
  d[d <= 0] <- 1
  Dinv <- diag(1 / d, length(d))
  R <- Dinv %*% S %*% Dinv
  R <- (R + t(R)) / 2
  diag(R) <- 1
  rownames(R) <- trait_names
  colnames(R) <- trait_names
  R
}

ridge_tuning <- function(eigvals) {
  eig_min <- min(eigvals)
  tuning_par <- max(0, -2 * eig_min)
  tuning_par + eigvals
}

eigen_ratios <- function(eigvals, cutoff = NULL) {
  d <- length(eigvals)
  if (is.null(cutoff)) cutoff <- floor(d / 2) - 1
  vals <- 1 / eigvals[d:1]
  vals[1:cutoff] / vals[2:(cutoff + 1)]
}

optimal_gap <- function(eigen_rat, overest_type, d) {
  sortv <- sort(eigen_rat, decreasing = TRUE)
  max1 <- sortv[1]
  max2 <- sortv[2]
  
  if (is.numeric(overest_type)) return(list(shat = overest_type, overest = TRUE))
  if (overest_type == "true")   return(list(shat = floor(sqrt(d)), overest = TRUE))
  if (max1 > 1.5 * max2)        return(list(shat = which.max(eigen_rat), overest = FALSE))
  if (overest_type == "frac")   return(list(shat = floor(d / 5), overest = TRUE))
  if (overest_type == "sqrt")   return(list(shat = floor(sqrt(d)), overest = TRUE))
  
  list(shat = floor(d / 5), overest = TRUE)
}

influence_measures <- function(eig_mat, shat) {
  d <- length(eig_mat$values)
  if (shat == 1) {
    eig_mat$vectors[, d]^2
  } else {
    vectors <- eig_mat$vectors[, d:(d - shat + 1), drop = FALSE]
    apply(vectors^2, 1, sum)
  }
}

ipchd_score_from_mat <- function(mat, overest_type = "frac") {
  eig_mat <- eigen(mat, symmetric = TRUE)
  ridge_eig <- ridge_tuning(eig_mat$values)
  eigen_rat <- eigen_ratios(ridge_eig)
  shat <- optimal_gap(eigen_rat, overest_type, ncol(mat))$shat
  score <- influence_measures(eig_mat, shat)
  
  list(
    score = score,
    shat = shat,
    eigvals = eig_mat$values,
    eigen_ratios = eigen_rat
  )
}

read_one_trait_subset <- function(file, snp_index) {
  dat <- fread(
    file,
    select = c("variant", "beta", "se", "pval"),
    showProgress = FALSE
  )
  setkey(dat, variant)
  dat2 <- dat[snp_index, nomatch = 0L]
  dat2[, .(idx, beta, se, pval)]
}

build_effect_matrices <- function(snp_vec, file_map, prefix, out_dir) {
  pheno_codes_local <- file_map$phenotype
  n_traits <- length(pheno_codes_local)
  m <- length(snp_vec)
  
  beta_file     <- file.path(out_dir, paste0("beta_mat_", prefix, ".rds"))
  se_file       <- file.path(out_dir, paste0("se_mat_", prefix, ".rds"))
  pval_file     <- file.path(out_dir, paste0("pval_mat_", prefix, ".rds"))
  progress_file <- file.path(out_dir, paste0(prefix, "_matrix_progress.csv"))
  info_file     <- file.path(out_dir, paste0(prefix, "_matrix_info.csv"))
  
  log_msg("Building %s matrices: m = %d SNPs, p = %d traits", prefix, m, n_traits)
  
  snp_index <- data.table(variant = snp_vec, idx = seq_len(m))
  setkey(snp_index, variant)
  
  beta_mat <- matrix(NA_real_, nrow = m, ncol = n_traits,
                     dimnames = list(snp_vec, pheno_codes_local))
  se_mat   <- matrix(NA_real_, nrow = m, ncol = n_traits,
                     dimnames = list(snp_vec, pheno_codes_local))
  pval_mat <- matrix(NA_real_, nrow = m, ncol = n_traits,
                     dimnames = list(snp_vec, pheno_codes_local))
  
  progress_dt <- data.table(
    phenotype = pheno_codes_local,
    file_path = file_map$file_path,
    n_rows_matched = NA_integer_
  )
  
  for (j in seq_len(n_traits)) {
    code <- pheno_codes_local[j]
    file <- file_map$file_path[j]
    
    log_msg("[%s matrix] trait %d / %d: %s", prefix, j, n_traits, code)
    
    dtj <- read_one_trait_subset(file, snp_index)
    setorder(dtj, idx)
    
    if (nrow(dtj) != m) {
      stop(sprintf("[%s matrix] trait %s matched %d SNPs, expected %d",
                   prefix, code, nrow(dtj), m))
    }
    
    beta_mat[, j] <- dtj$beta
    se_mat[, j]   <- dtj$se
    pval_mat[, j] <- dtj$pval
    
    progress_dt[j, n_rows_matched := nrow(dtj)]
    fwrite(progress_dt, progress_file)
    
    rm(dtj)
    gc(verbose = FALSE)
  }
  
  saveRDS(beta_mat, beta_file)
  saveRDS(se_mat, se_file)
  saveRDS(pval_mat, pval_file)
  
  info_dt <- data.table(
    prefix = prefix,
    n_snps = m,
    n_traits = n_traits,
    beta_dim1 = nrow(beta_mat),
    beta_dim2 = ncol(beta_mat),
    se_dim1 = nrow(se_mat),
    se_dim2 = ncol(se_mat),
    pval_dim1 = nrow(pval_mat),
    pval_dim2 = ncol(pval_mat)
  )
  fwrite(info_dt, info_file)
  
  log_msg("Saved %s matrices", prefix)
  log_msg("  %s", beta_file)
  log_msg("  %s", se_file)
  log_msg("  %s", pval_file)
  log_msg("  %s", progress_file)
  log_msg("  %s", info_file)
  
  invisible(list(
    beta_file = beta_file,
    se_file = se_file,
    pval_file = pval_file,
    progress_file = progress_file,
    info_file = info_file
  ))
}

build_z_matrix <- function(snp_vec, file_map, prefix, out_dir) {
  pheno_codes_local <- file_map$phenotype
  n_traits <- length(pheno_codes_local)
  m <- length(snp_vec)
  
  z_file        <- file.path(out_dir, paste0("z_mat_", prefix, ".rds"))
  progress_file <- file.path(out_dir, paste0(prefix, "_z_matrix_progress.csv"))
  info_file     <- file.path(out_dir, paste0(prefix, "_z_matrix_info.csv"))
  
  log_msg("Building %s z matrix: m = %d SNPs, p = %d traits", prefix, m, n_traits)
  
  snp_index <- data.table(variant = snp_vec, idx = seq_len(m))
  setkey(snp_index, variant)
  
  z_mat <- matrix(NA_real_, nrow = m, ncol = n_traits,
                  dimnames = list(snp_vec, pheno_codes_local))
  
  progress_dt <- data.table(
    phenotype = pheno_codes_local,
    file_path = file_map$file_path,
    n_rows_matched = NA_integer_
  )
  
  for (j in seq_len(n_traits)) {
    code <- pheno_codes_local[j]
    file <- file_map$file_path[j]
    
    log_msg("[%s z matrix] trait %d / %d: %s", prefix, j, n_traits, code)
    
    dtj <- read_one_trait_subset(file, snp_index)
    setorder(dtj, idx)
    
    if (nrow(dtj) != m) {
      stop(sprintf("[%s z matrix] trait %s matched %d SNPs, expected %d",
                   prefix, code, nrow(dtj), m))
    }
    
    z_mat[, j] <- dtj$beta / dtj$se
    
    progress_dt[j, n_rows_matched := nrow(dtj)]
    fwrite(progress_dt, progress_file)
    
    rm(dtj)
    gc(verbose = FALSE)
  }
  
  saveRDS(z_mat, z_file)
  
  info_dt <- data.table(
    prefix = prefix,
    n_snps = m,
    n_traits = n_traits,
    z_dim1 = nrow(z_mat),
    z_dim2 = ncol(z_mat),
    all_finite = all(is.finite(z_mat))
  )
  fwrite(info_dt, info_file)
  
  log_msg("Saved %s z matrix", prefix)
  log_msg("  %s", z_file)
  log_msg("  %s", progress_file)
  log_msg("  %s", info_file)
  
  invisible(list(
    z_file = z_file,
    progress_file = progress_file,
    info_file = info_file
  ))
}

# =========================================================
# 3. File map
# =========================================================
log_msg("Step 1/12: discovering trait files")
# Discover example1.tsv, example2.tsv, ...; no manually entered trait count.
input_files <- list.files(data_dir, pattern = "^example[1-9][0-9]*\\.tsv$", full.names = TRUE)
input_files <- input_files[!file.info(input_files)$isdir]
example_index <- as.integer(sub("^example([0-9]+)\\.tsv$", "\\1", basename(input_files)))
input_files <- input_files[order(example_index)]
if (length(input_files) < 6L) stop("At least six exampleN.tsv trait files are required; found ", length(input_files))
pheno_codes <- sub("\\.tsv$", "", basename(input_files))
if (anyDuplicated(pheno_codes)) stop("Duplicate trait names after removing extensions. Keep exactly one file per trait.")
files <- setNames(input_files, pheno_codes)
file_map <- data.table(phenotype = pheno_codes, file_path = unname(files))
file_map_file <- file.path(out_dir, "file_map.csv")
fwrite(file_map, file_map_file)
required_columns <- c("variant", "beta", "se", "pval", "minor_AF", "low_confidence_variant")
for (f in input_files) {
  header <- names(fread(f, nrows = 0L, showProgress = FALSE))
  missing_columns <- setdiff(required_columns, header)
  if (length(missing_columns)) stop(basename(f), ": missing columns: ", paste(missing_columns, collapse = ", "))
}
log_msg("Automatically detected %d traits", length(pheno_codes))
log_msg("Step 2/12: finding common valid SNPs across detected traits")

get_valid_variants <- function(file) {
  dat <- fread(
    file,
    select = c("variant", "beta", "se", "pval"),
    showProgress = FALSE
  )
  
  keep <- !is.na(dat$beta) &
    !is.na(dat$se) &
    !is.na(dat$pval) &
    is.finite(dat$beta) &
    is.finite(dat$se) &
    is.finite(dat$pval) & dat$se > 0 & dat$pval >= 0 & dat$pval <= 1
  
  if (anyDuplicated(dat$variant)) stop("Duplicate variant IDs in ", file)
  unique(dat$variant[keep])
}

variant_list <- vector("list", length(pheno_codes))
names(variant_list) <- pheno_codes

qc_summary <- data.table(
  phenotype = pheno_codes,
  file_path = unname(files),
  n_valid_snps = NA_integer_
)

qc_counts_file   <- file.path(out_dir, "qc_valid_snp_counts.csv")
common_snps_csv  <- file.path(out_dir, "common_snps.csv")
common_snps_rds  <- file.path(out_dir, "common_snps.rds")
qc_overview_file <- file.path(out_dir, "qc_overview.csv")

for (i in seq_along(pheno_codes)) {
  code <- pheno_codes[i]
  log_msg("[common SNP] trait %d / %d: %s", i, length(pheno_codes), code)
  
  vars <- get_valid_variants(files[[code]])
  variant_list[[code]] <- vars
  qc_summary[i, n_valid_snps := length(vars)]
  fwrite(qc_summary, qc_counts_file)
  
  log_msg("  valid SNPs = %d", length(vars))
  gc(verbose = FALSE)
}

common_snps <- Reduce(intersect, variant_list)
if (!length(common_snps)) stop("No common valid SNPs across trait files.")

fwrite(qc_summary, qc_counts_file)
fwrite(data.table(variant = common_snps), common_snps_csv)
saveRDS(common_snps, common_snps_rds)

qc_overview <- data.table(
  n_traits = length(pheno_codes),
  min_valid_snps = min(qc_summary$n_valid_snps, na.rm = TRUE),
  median_valid_snps = median(qc_summary$n_valid_snps, na.rm = TRUE),
  max_valid_snps = max(qc_summary$n_valid_snps, na.rm = TRUE),
  n_common_snps = length(common_snps)
)
fwrite(qc_overview, qc_overview_file)

log_msg("Final common SNP count = %d", length(common_snps))
log_msg("Saved common SNP outputs")

# =========================================================
# 5. Stage 1 candidate SNPs
# =========================================================
log_msg("Step 3/12: stage 1 candidate SNP filtering")

ref_file <- file_map$file_path[1]

candidate_csv  <- file.path(out_dir, "candidate_snps_stage1.csv")
candidate_rds  <- file.path(out_dir, "candidate_snps_stage1.rds")
candidate_info <- file.path(out_dir, "candidate_snps_stage1_info.csv")

common_dt <- data.table(variant = common_snps)
setkey(common_dt, variant)

ref_dt <- fread(
  ref_file,
  select = c("variant", "minor_AF", "low_confidence_variant"),
  showProgress = FALSE
)
setkey(ref_dt, variant)

dt <- ref_dt[common_dt, nomatch = 0L]
tmp <- tstrsplit(dt$variant, ":", fixed = TRUE)
dt[, chr := tmp[[1]]]
if (length(tmp) < 2L) stop("variant IDs must begin chromosome:position, e.g. 1:100000:A:G")
dt[, pos := as.integer(tmp[[2]])]
if (anyNA(dt$pos) || anyNA(dt$chr)) stop("Invalid chromosome:position variant IDs.")

dt_filt <- dt[
  !is.na(minor_AF) &
    is.finite(minor_AF) &
    minor_AF >= maf_thresh &
    low_confidence_variant == FALSE
]

log_msg("After MAF + low_confidence filter: %d", nrow(dt_filt))

if (remove_mhc) {
  before_mhc <- nrow(dt_filt)
  dt_filt <- dt_filt[!(chr == "6" & pos >= 25000000 & pos <= 34000000)]
  log_msg("Removed MHC SNPs: %d", before_mhc - nrow(dt_filt))
}

setorder(dt_filt, chr, pos)

fwrite(dt_filt, candidate_csv)
saveRDS(dt_filt$variant, candidate_rds)

info_dt <- data.table(
  common_snps = nrow(common_dt),
  matched_in_reference = nrow(dt),
  maf_threshold = maf_thresh,
  removed_low_confidence = TRUE,
  removed_mhc = remove_mhc,
  final_candidate_snps = nrow(dt_filt)
)
fwrite(info_dt, candidate_info)

log_msg("Final candidate SNP count after stage 1 = %d", nrow(dt_filt))

# =========================================================
# 6. Stage 2 SNP summary across detected traits
# =========================================================
log_msg("Step 4/12: computing stage 2 SNP-level summary across all traits")

candidate_snps <- readRDS(candidate_rds)
n_snps <- length(candidate_snps)
if (n_snps < 3L) stop("Too few common SNPs remain after QC.")
n_traits <- nrow(file_map)

summary_csv   <- file.path(out_dir, "candidate_snps_stage2_summary.csv")
summary_rds   <- file.path(out_dir, "candidate_snps_stage2_summary.rds")
progress_file <- file.path(out_dir, "candidate_snps_stage2_progress_by_trait.csv")
info_file     <- file.path(out_dir, "candidate_snps_stage2_summary_info.csv")

candidate_index <- data.table(
  variant = candidate_snps,
  idx = seq_len(n_snps)
)
setkey(candidate_index, variant)

all_finite <- rep(TRUE, n_snps)
min_pval   <- rep(1, n_snps)

progress_dt <- data.table(
  phenotype = file_map$phenotype,
  file_path = file_map$file_path,
  n_rows_matched = NA_integer_
)

for (j in seq_len(n_traits)) {
  code <- file_map$phenotype[j]
  file <- file_map$file_path[j]
  
  log_msg("[stage2 summary] trait %d / %d: %s", j, n_traits, code)
  
  dtj <- read_one_trait_subset(file, candidate_index)
  progress_dt[j, n_rows_matched := nrow(dtj)]
  fwrite(progress_dt, progress_file)
  
  finite_now <- is.finite(dtj$beta) & is.finite(dtj$se) & is.finite(dtj$pval)
  finite_trait <- rep(FALSE, n_snps)
  finite_trait[dtj$idx] <- finite_now
  all_finite <- all_finite & finite_trait
  
  valid_p <- is.finite(dtj$pval)
  idx_p <- dtj$idx[valid_p]
  if (length(idx_p) > 0) {
    min_pval[idx_p] <- pmin(min_pval[idx_p], dtj$pval[valid_p])
  }
  
  log_msg("  matched = %d | finite this trait = %d", nrow(dtj), sum(finite_now))
  
  rm(dtj, finite_now, finite_trait, valid_p, idx_p)
  gc(verbose = FALSE)
}

summary_dt <- data.table(
  variant = candidate_snps,
  all_finite = all_finite,
  min_pval = min_pval
)

fwrite(summary_dt, summary_csv)
saveRDS(summary_dt, summary_rds)

info_dt <- data.table(
  n_candidate_snps = n_snps,
  n_traits = n_traits,
  n_all_finite = sum(summary_dt$all_finite, na.rm = TRUE),
  min_of_min_pval = min(summary_dt$min_pval, na.rm = TRUE),
  median_of_min_pval = median(summary_dt$min_pval, na.rm = TRUE),
  max_of_min_pval = max(summary_dt$min_pval, na.rm = TRUE)
)
fwrite(info_dt, info_file)

log_msg("Saved stage 2 SNP summary")

# =========================================================
# 7. Final null SNPs
# =========================================================
log_msg("Step 5/12: selecting final null SNPs")

null_csv     <- file.path(out_dir, "null_snps_final.csv")
null_rds     <- file.path(out_dir, "null_snps_final.rds")
info_csv     <- file.path(out_dir, "final_null_snp_set_info.csv")

sum_dt <- fread(summary_csv)

null_final <- sum_dt[
  all_finite == TRUE &
    min_pval > null_min_pval,
  .(variant)
]
setorder(null_final, variant)

fwrite(null_final, null_csv)
saveRDS(null_final$variant, null_rds)

info_dt <- data.table(
  n_total_stage2 = nrow(sum_dt),
  null_rule_min_pval_gt = null_min_pval,
  null_final_n = nrow(null_final)
)
fwrite(info_dt, info_csv)

log_msg("null_final_n = %d", nrow(null_final))
if (nrow(null_final) == 0) {
  stop("No null SNPs left after applying all-trait p > 0.05 criterion.")
}

# =========================================================
# 8. Build null effect matrices
# =========================================================
log_msg("Step 6/12: building null effect matrices")
null_snps <- readRDS(null_rds)
build_effect_matrices(
  snp_vec = null_snps,
  file_map = file_map,
  prefix = "null",
  out_dir = out_dir
)

# =========================================================
# 9. Estimate Sigma_omega from null z matrix
# =========================================================
log_msg("Step 7/12: estimating Sigma_omega from null SNP z-score matrix")

beta_file <- file.path(out_dir, "beta_mat_null.rds")
se_file   <- file.path(out_dir, "se_mat_null.rds")

omega_rds_file  <- file.path(out_dir, "Sigma_omega_z_hat.rds")
omega_csv_file  <- file.path(out_dir, "Sigma_omega_z_hat.csv")
omega_info_file <- file.path(out_dir, "Sigma_omega_z_hat_info.csv")
omega_diag_file <- file.path(out_dir, "Sigma_omega_z_hat_diag.csv")

beta_mat <- readRDS(beta_file)
se_mat   <- readRDS(se_file)

stopifnot(all(dim(beta_mat) == dim(se_mat)))

m_null <- nrow(beta_mat)
p <- ncol(beta_mat)

log_msg("Number of null SNPs = %d", m_null)
log_msg("Number of traits = %d", p)

z_mat <- beta_mat / se_mat

finite_all <- all(is.finite(z_mat))
log_msg("All null z finite = %s", finite_all)
if (!finite_all) stop("z_mat contains non-finite values.")

if (m_null < 10L) stop("At least 10 common null SNPs are required for the 10% subsampling demonstration.")
n_sub <- ceiling(m_null * subsample_prop)
log_msg("Subsample size per repetition = %d", n_sub)

Sigma_sum <- matrix(0, nrow = p, ncol = p)

set.seed(seed_global)
for (b in seq_len(n_rep_omega)) {
  if (b %% 10 == 0 || b == 1 || b == n_rep_omega) {
    log_msg("[Sigma_omega] repetition %d / %d", b, n_rep_omega)
  }
  
  idx <- sample.int(m_null, size = n_sub, replace = FALSE)
  Z_sub <- z_mat[idx, , drop = FALSE]
  Sigma_b <- crossprod(Z_sub) / n_sub
  Sigma_sum <- Sigma_sum + Sigma_b
}

Sigma_omega_hat <- Sigma_sum / n_rep_omega
Sigma_omega_hat <- (Sigma_omega_hat + t(Sigma_omega_hat)) / 2

trait_names <- colnames(z_mat)
rownames(Sigma_omega_hat) <- trait_names
colnames(Sigma_omega_hat) <- trait_names

diag_vals <- diag(Sigma_omega_hat)
eig_vals  <- eigen(Sigma_omega_hat, symmetric = TRUE, only.values = TRUE)$values

omega_info <- data.table(
  n_null_snps = m_null,
  n_traits = p,
  subsample_prop = subsample_prop,
  n_sub = n_sub,
  n_rep = n_rep_omega,
  diag_min = min(diag_vals),
  diag_median = median(diag_vals),
  diag_max = max(diag_vals),
  eig_min = min(eig_vals),
  eig_median = median(eig_vals),
  eig_max = max(eig_vals)
)
omega_diag <- data.table(
  trait = trait_names,
  diag_value = diag_vals
)

saveRDS(Sigma_omega_hat, omega_rds_file)
fwrite(as.data.table(Sigma_omega_hat, keep.rownames = "trait"), omega_csv_file)
fwrite(omega_info, omega_info_file)
fwrite(omega_diag, omega_diag_file)

log_msg("Saved Sigma_omega outputs")

# =========================================================
# 10. Joint chi-square screening of signal SNPs with LD clumping
# =========================================================
log_msg("Step 8/12: joint chi-square screening")

candidate_z_prefix <- "candidate"
candidate_z_rds    <- file.path(out_dir, paste0("z_mat_", candidate_z_prefix, ".rds"))
joint_csv          <- file.path(out_dir, "joint_chisq_screening_results.csv")
joint_rds          <- file.path(out_dir, "joint_chisq_screening_results.rds")
joint_sig_csv      <- file.path(out_dir, "joint_significant_snps.csv")
clump_out_csv      <- file.path(out_dir, "ld_clump_output.csv")
clump_input_csv    <- file.path(out_dir, "ld_clump_input_rsid_pval.csv")
analysis_csv       <- file.path(out_dir, "analysis_snps_final.csv")
analysis_rds       <- file.path(out_dir, "analysis_snps_final.rds")
analysis_info_csv  <- file.path(out_dir, "analysis_snps_final_info.csv")

candidate_snps <- readRDS(candidate_rds)
build_z_matrix(
  snp_vec = candidate_snps,
  file_map = file_map,
  prefix = candidate_z_prefix,
  out_dir = out_dir
)

z_mat_candidate <- readRDS(candidate_z_rds)
if (!all(is.finite(z_mat_candidate))) {
  stop("candidate z matrix contains non-finite values.")
}

Sigma_omega_for_test <- make_spd(Sigma_omega_hat)
Omega_inv <- solve(Sigma_omega_for_test)

joint_stat <- rowSums((z_mat_candidate %*% Omega_inv) * z_mat_candidate)
joint_pval <- pchisq(joint_stat, df = ncol(z_mat_candidate), lower.tail = FALSE)

joint_dt <- data.table(
  variant = rownames(z_mat_candidate),
  joint_chisq = joint_stat,
  joint_pval = joint_pval
)
setorder(joint_dt, joint_pval)

fwrite(joint_dt, joint_csv)
saveRDS(joint_dt, joint_rds)

joint_sig <- joint_dt[joint_pval < joint_p_threshold]
fwrite(joint_sig, joint_sig_csv)

log_msg("joint significant SNPs before final selection = %d", nrow(joint_sig))
if (nrow(joint_sig) == 0) {
  stop("No SNPs passed the joint test. Check input signal strength and sample size; thresholds are not relaxed automatically.")
}

# Map coordinate/allele variant IDs to the reference BIM IDs.
# Exact or swapped allele pairs only; no strand or genome-build conversion.
bim <- fread(paste0(ld_bfile, ".bim"), header = FALSE, showProgress = FALSE)
if (ncol(bim) != 6L) stop("LD reference BIM must have six columns.")
setnames(bim, c("chr", "rsid", "cm", "pos", "a1", "a2"))
norm_chr <- function(x) {
  x <- toupper(sub("^chr", "", as.character(x), ignore.case = TRUE))
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x %in% c("M", "MT")] <- "26"
  x
}
make_key <- function(chr, pos, a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  paste(norm_chr(chr), as.integer(pos), pmin(a1, a2), pmax(a1, a2), sep = ":")
}
bim[, key := make_key(chr, pos, a1, a2)]
bim <- bim[!is.na(rsid) & rsid != "" & rsid != "."]
# Exclude ambiguous loci rather than selecting an arbitrary reference ID.
bim <- unique(bim[, .(key, rsid)])
ambiguous <- bim[, .N, by = key][N > 1L, key]
bim <- bim[!key %in% ambiguous]
parts <- tstrsplit(joint_sig$variant, ":", fixed = TRUE)
if (length(parts) != 4L) stop("variant must be chromosome:position:allele1:allele2 for BIM matching.")
query <- data.table(variant = joint_sig$variant,
                    key = make_key(parts[[1]], parts[[2]], parts[[3]], parts[[4]]))
variant_dict <- merge(query, bim, by = "key", all.x = TRUE)[, .(variant, rsid)]
setkey(variant_dict, variant)
setDT(joint_sig)
setkey(joint_sig, variant)

joint_sig_map <- variant_dict[joint_sig]
setnames(joint_sig_map, c("variant", "rsid", "joint_chisq", "joint_pval"))
joint_sig_map <- joint_sig_map[!is.na(rsid) & rsid != ""]

if (nrow(joint_sig_map) == 0) {
  stop("No jointly significant SNPs could be mapped to rsid.")
}

joint_sig_map <- joint_sig_map[order(joint_pval)]
clump_input <- joint_sig_map[, .SD[1], by = rsid]
clump_input <- clump_input[, .(rsid, pval = joint_pval, variant)]
fwrite(clump_input, clump_input_csv)

if (!file.exists(plink_bin)) {
  stop("plink executable not found: ", plink_bin)
}
if (!file.exists(paste0(ld_bfile, ".bed")) ||
    !file.exists(paste0(ld_bfile, ".bim")) ||
    !file.exists(paste0(ld_bfile, ".fam"))) {
  stop("LD reference panel files not found for prefix: ", ld_bfile)
}
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  stop("Package 'ieugwasr' is required for ld_clump_local().")
}

clump_res <- try(
  ieugwasr::ld_clump_local(
    dat = clump_input[, .(rsid, pval)],
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = clump_p,
    bfile = ld_bfile,
    plink_bin = plink_bin
  ),
  silent = TRUE
)

if (inherits(clump_res, "try-error")) {
  stop("LD clumping failed: ", as.character(clump_res))
}

clump_res <- as.data.table(clump_res)
if (!"rsid" %in% names(clump_res)) {
  stop("ld_clump_local output does not contain column 'rsid'.")
}

analysis_final <- merge(
  clump_res[, .(rsid)],
  clump_input[, .(rsid, variant, pval)],
  by = "rsid",
  all.x = TRUE,
  sort = FALSE
)
analysis_final <- unique(analysis_final[, .(variant, rsid, joint_pval = pval)])
setorder(analysis_final, joint_pval, variant)

fwrite(clump_res, clump_out_csv)
fwrite(analysis_final, analysis_csv)
saveRDS(analysis_final$variant, analysis_rds)

analysis_info <- data.table(
  ld_mode = "plink",
  ld_clumping_performed = TRUE,
  n_candidate_stage1 = length(candidate_snps),
  n_joint_significant = nrow(joint_sig),
  n_analysis_final = nrow(analysis_final),
  joint_p_threshold = joint_p_threshold,
  clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = clump_p
)
fwrite(analysis_info, analysis_info_csv)

log_msg("analysis_final_n after LD clumping = %d", nrow(analysis_final))
if (nrow(analysis_final) < 3L) {
  stop("At least three analysis SNPs are required for robust covariance estimation.")
}

# =========================================================
# 11. Build analysis effect matrices and estimate Sigma_beta
# =========================================================
log_msg("Step 9/12: building analysis effect matrices")
analysis_snps <- readRDS(analysis_rds)
build_effect_matrices(
  snp_vec = analysis_snps,
  file_map = file_map,
  prefix = "analysis",
  out_dir = out_dir
)

log_msg("Step 10/12: estimating Sigma_beta")

beta_file  <- file.path(out_dir, "beta_mat_analysis.rds")
se_file    <- file.path(out_dir, "se_mat_analysis.rds")
omega_file <- file.path(out_dir, "Sigma_omega_z_hat.rds")

sigma_hatbeta_file <- file.path(out_dir, "Sigma_hatbeta_robust_z.rds")
sigma_beta_file    <- file.path(out_dir, "Sigma_beta_hat_final.rds")
info_file_sigma    <- file.path(out_dir, "Sigma_beta_hat_final_info.csv")

beta_mat <- readRDS(beta_file)
se_mat   <- readRDS(se_file)
Sigma_omega_hat <- readRDS(omega_file)

stopifnot(all(dim(beta_mat) == dim(se_mat)))
stopifnot(ncol(beta_mat) == nrow(Sigma_omega_hat))
stopifnot(ncol(beta_mat) == ncol(Sigma_omega_hat))

z_mat_analysis <- beta_mat / se_mat
if (!all(is.finite(z_mat_analysis))) {
  stop("z_mat_analysis contains non-finite values.")
}

if (any(apply(z_mat_analysis, 2, mad_sd) <= 0)) stop("At least one trait has zero MAD among selected signal SNPs; provide variable signal data.")
log_msg("Computing robust observed covariance from analysis z matrix")
Sigma_hatbeta_rob <- spearmancov(z_mat_analysis)
Sigma_hatbeta_rob <- make_spd(Sigma_hatbeta_rob)

log_msg("Subtracting Sigma_omega and projecting to SPD")
Sigma_beta_hat <- Sigma_hatbeta_rob - Sigma_omega_hat
Sigma_beta_hat <- make_spd(Sigma_beta_hat)

trait_names <- colnames(z_mat_analysis)
rownames(Sigma_hatbeta_rob) <- trait_names
colnames(Sigma_hatbeta_rob) <- trait_names
rownames(Sigma_beta_hat) <- trait_names
colnames(Sigma_beta_hat) <- trait_names

saveRDS(Sigma_hatbeta_rob, sigma_hatbeta_file)
saveRDS(Sigma_beta_hat, sigma_beta_file)

info_dt <- data.table(
  n_analysis_snps = nrow(z_mat_analysis),
  n_traits = ncol(z_mat_analysis),
  hatbeta_eig_min = min(eigen(Sigma_hatbeta_rob, symmetric = TRUE, only.values = TRUE)$values),
  hatbeta_eig_max = max(eigen(Sigma_hatbeta_rob, symmetric = TRUE, only.values = TRUE)$values),
  beta_eig_min = min(eigen(Sigma_beta_hat, symmetric = TRUE, only.values = TRUE)$values),
  beta_eig_max = max(eigen(Sigma_beta_hat, symmetric = TRUE, only.values = TRUE)$values)
)
fwrite(info_dt, info_file_sigma)

log_msg("Saved Sigma_beta outputs")

# =========================================================
# 12. IPC-HD scores and hub calling
# =========================================================
log_msg("Step 11/12: IPC-HD hub scoring and hub selection")

hub_score_file  <- file.path(out_dir, "ipchd_hub_scores_realdata_with_traits.csv")
hub_result_file <- file.path(out_dir, "ipchd_hub_scores_realdata_with_traits_and_hubflag.csv")
hub_info_file   <- file.path(out_dir, "ipchd_hub_summary_info.csv")

Sigma_beta_hat <- readRDS(sigma_beta_file)

R_beta <- cov2cor_safe(Sigma_beta_hat)
res <- ipchd_score_from_mat(R_beta, overest_type = "frac")
saveRDS(R_beta, file.path(out_dir, "genetic_correlation.rds"))
eigenvalues <- eigen(R_beta, symmetric = TRUE, only.values = TRUE)$values
shifted <- ridge_tuning(eigenvalues)
ratios <- eigen_ratios(shifted)
fwrite(data.table(index = seq_along(eigenvalues), inverse_corr_eigenvalue = 1 / rev(eigenvalues)), file.path(out_dir, "spectrum.csv"))
fwrite(data.table(index = seq_along(ratios), eigen_ratio = ratios), file.path(out_dir, "eigenvalue_ratios.csv"))
fwrite(data.table(ld_mode = "plink", n_traits = ncol(R_beta), ld_clumping_performed = TRUE), file.path(out_dir, "run_metadata.csv"))

hub_score <- data.table(
  trait = colnames(R_beta),
  ipchd_score = res$score
)
setorder(hub_score, -ipchd_score)

thr <- mean(hub_score$ipchd_score) + hub_sd_multiplier * sd(hub_score$ipchd_score)
hub_score[, is_hub := ipchd_score > thr]

fwrite(hub_score[, .(trait, ipchd_score)], hub_score_file)
fwrite(hub_score, hub_result_file)

hub_info <- data.table(
  estimated_shat = res$shat,
  score_mean = mean(hub_score$ipchd_score),
  score_sd = sd(hub_score$ipchd_score),
  threshold = thr,
  n_hubs = sum(hub_score$is_hub)
)
fwrite(hub_info, hub_info_file)

log_msg("Estimated shat = %d", res$shat)
log_msg("Hub threshold = %.8f", thr)
log_msg("Number of hubs = %d", sum(hub_score$is_hub))
log_msg("Top 20 traits by IPC-HD score:")
print(head(hub_score, 20L))

log_msg("Step 12/12: all results saved under: %s", out_dir)
log_msg("Pipeline completed successfully")

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
}

# ---- Combined trait-score, spectrum and ratio figure ----
if (run_plots) {
# Combined adaptation of code_traits.R and code_ratio.R.
# One automatically sized panel; no embedded scores, eigenvalues or trait counts.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
result_dir <- "results/real_data"
read_result <- function(name) {
  path <- file.path(result_dir, name)
  if (!file.exists(path)) stop("Missing result: ", path, ". Run real_data.R first.")
  fread(path)
}
d <- read_result("ipchd_hub_scores_realdata_with_traits_and_hubflag.csv")
info <- read_result("ipchd_hub_summary_info.csv")
spectrum <- read_result("spectrum.csv")
ratios <- read_result("eigenvalue_ratios.csv")
metadata <- read_result("run_metadata.csv")
if (!all(c("trait", "ipchd_score", "is_hub") %in% names(d))) stop("Invalid score table.")
if (nrow(info) != 1L || nrow(metadata) != 1L) stop("Invalid run metadata.")
if (anyDuplicated(d$trait) || any(!is.finite(d$ipchd_score)) || any(d$ipchd_score < 0)) stop("Invalid scores.")
p <- nrow(d)
if (p != metadata$n_traits || nrow(spectrum) != p) stop("Result dimensions disagree; do not mix outputs from different runs.")
if (any(!is.finite(spectrum$inverse_corr_eigenvalue)) || any(spectrum$inverse_corr_eigenvalue <= 0) || any(!is.finite(ratios$eigen_ratio))) stop("Invalid spectrum or ratios.")
cutoff <- info$threshold
shat <- info$estimated_shat
if (!is.finite(cutoff) || cutoff <= 0 || shat < 1 || shat > p) stop("Invalid cutoff or selected dimension.")
if (anyNA(d$is_hub) || any(d$is_hub != (d$ipchd_score > cutoff))) stop("Hub flags and stored cutoff disagree.")

# Preserve the original log-score display, but label its +12 offset explicitly.
# Exact zeros cannot be logged; floor for display only and state this in caption.
score_floor <- 1e-12
d[, display_score := log(pmax(ipchd_score, score_floor)) + 12]
setorder(d, -ipchd_score, trait)
d[, Trait := factor(trait, levels = trait)]
d[, group := factor(ifelse(is_hub, "High", "Normal"), levels = c("High", "Normal"))]
base_y <- min(0, min(d$display_score) - 0.5)
if (!isTRUE(metadata$ld_clumping_performed)) stop("Expected results from the normal LD-clumping workflow; rerun analysis.")
caption <- "Signal SNPs were LD-clumped. Panel-specific centrality is not a causal effect estimate."
if (any(d$ipchd_score == 0)) caption <- paste(caption, "Zero scores displayed at 1e-12 before logarithmic transformation.")
theme_common <- theme_minimal(base_size = 11) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(color = "#DDE6EB"),
  axis.line = element_line(color = "#4A4A4A"),
  plot.title = element_text(hjust = 0.5), legend.title = element_blank())

score_plot <- ggplot(d, aes(Trait, display_score, color = group)) +
  geom_segment(aes(xend = Trait, y = base_y, yend = display_score), linewidth = 1.6, alpha = 0.5) +
  geom_hline(yintercept = log(cutoff) + 12, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  geom_text_repel(data = d[group == "High"], aes(label = trait), direction = "y", seed = 123, max.overlaps = Inf, show.legend = FALSE, size = 3) +
  scale_color_manual(values = c(High = "#F26B5E", Normal = "#3D5A80"), drop = FALSE) +
  labs(x = NULL, y = "ln(HINGE score) + 12", title = paste("Trait scores (", p, " traits)", sep = "")) +
  theme_common + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

spectrum[, selected := index <= shat]
spectrum_plot <- ggplot(spectrum, aes(index, inverse_corr_eigenvalue)) +
  geom_line(color = "#5B7C99") +
  geom_point(aes(color = selected), size = 1.8) +
  geom_vline(xintercept = shat + 0.5, linetype = "dashed", color = "#ABABAB") +
  scale_color_manual(values = c("TRUE" = "#EF9F27", "FALSE" = "#5B7C99"), guide = "none") +
  scale_y_log10() +
  labs(x = "Eigenvalue index", y = "Inverse-correlation eigenvalue (log10 axis)", title = "Inverse correlation spectrum", subtitle = paste("Selected dimension:", shat)) + theme_common

gap_index <- ratios$index[which.max(ratios$eigen_ratio)]
ratio_plot <- ggplot(ratios, aes(index, eigen_ratio)) +
  geom_line(color = "#5B7C99") + geom_point(color = "#5B7C99") +
  geom_point(data = ratios[index == gap_index], color = "#EF9F27", size = 3) +
  geom_vline(xintercept = shat, color = "#ABABAB", linetype = "dashed") +
  labs(x = "Tail-subspace index", y = "Ridge-stabilized eigenvalue ratio", title = "Eigenvalue ratios",
       subtitle = paste("Largest ratio at", gap_index, "| selected dimension", shat)) + theme_common
# Largest-ratio index and selected dimension can differ when fallback is used.
combined <- (score_plot / (spectrum_plot | ratio_plot)) +
  plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A", caption = caption)
plot_dir <- "results/figures/real_data"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
width <- min(24, max(10, p * 0.20))
ggsave(file.path(plot_dir, "real_data_combined.pdf"), combined, width = width, height = 10)
ggsave(file.path(plot_dir, "real_data_combined.png"), combined, width = width, height = 10, dpi = 300)
capture.output(sessionInfo(), file = file.path(plot_dir, "sessionInfo.txt"))
message("Combined real-data plots saved to ", plot_dir)

}
