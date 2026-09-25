# HINGE simulation study and plotting. Run from the repository root.
# Set run_analysis <- FALSE to redraw figures from existing complete results.
run_analysis <- TRUE
run_plots <- TRUE

############################################################
## FULL PIPELINE SCRIPT (GRID SIMULATION, PARALLEL on Linux)
## Simulation functions, execution and plotting in one file.
##
## Changes made:
## 1) UHP rows: randomly sample 10% of m rows (not first 10%)
## 2) UHP traits: randomly sample 20 traits each replicate
## 3) Numerical matrix values: Uniform([-5,-4] U [4,5]) for all 3 edge types
## 4) Remove T0_prop = 0.25 from grid
## 5) UHPshift = 1
##
## GRID: 216 combos
##   T0_prop : 0.5, 0.75, 1
##   n       : 30000, 50000, 100000
##   m       : 500, 1000, 2000
##   UHP     : 0/1
##   ph      : 0.4/0.8
##   p       : 100/200
##
## Full run: 100 replicates per scenario
## Methods: raw, raw_deg, ip_thr, egg_theta_alpha, egg_theta_deg
############################################################

suppressPackageStartupMessages({
  library(MASS)
  library(glasso)
  library(parallel)
})

############################################################
# Small utilities
############################################################
vec   <- function(a) as.vector(a)
trace <- function(A) sum(diag(A))
rowMedian <- function(A) apply(A, 1, median)

############################################################
# EGG-aligned CScov + Nover
############################################################
CScov <- function(p, rho){
  matrix(rho, p, p) + (1 - rho) * diag(p)
}

make_Nover <- function(p){
  if (p %% 5 != 0) stop("EGG Nover requires p divisible by 5.")
  kronecker(CScov(5, 0.5), matrix(1, p / 5, p / 5))
}

############################################################
# 0) Missing author base functions: matrixEigen + matrixMultiply
############################################################
matrixEigen <- function(A) {
  A <- (A + t(A)) / 2
  eg <- eigen(A, symmetric = TRUE)
  list(
    values  = eg$values,
    vectors = eg$vectors,
    value   = eg$values,
    vector  = eg$vectors
  )
}

matrixMultiply <- function(A, B) A %*% B

############################################################
# 1) Author PSD/PD helpers: positveadj, cov2cor1, matrixsqrt
############################################################
positveadj <- function(A, min.eps = 0.001){
  a <- matrixEigen(A)
  d <- c(a$values)
  d[d < min.eps] <- 0
  matrixMultiply(a$vectors, t(a$vectors) * d)
}

cov2cor1 <- function(A, kappa = 100){
  fit <- matrixEigen(A)
  d <- c(fit$values)
  eps <- max(d) / kappa
  d[d < eps] <- eps
  B <- matrixMultiply(fit$vectors, t(fit$vectors) * d)
  colnames(B) <- rownames(B) <- colnames(A)
  cov2cor(B)
}

matrixsqrt <- function(A){
  fit <- matrixEigen(A)
  d <- c(fit$value)
  d[d < 0] <- 0
  d1 <- d * 0
  d1[d > 0] <- 1 / d[d > 0]
  d  <- sqrt(d)
  d1 <- sqrt(d1)
  W  <- matrixMultiply(fit$vector, t(fit$vector) * d)
  Wi <- matrixMultiply(fit$vector, t(fit$vector) * d1)
  list(w = W, wi = Wi)
}

############################################################
# 2) HubDetection 041: alphavals
############################################################
.alphavals <- function(Theta) apply(Theta^2, 2, sum)

############################################################
# 3) HubDetection 011: thresholding (needed by IPC-HD thr)
############################################################
thresholding <- function(mat, lambda) {
  nodiag <- mat - diag(diag(mat))
  matHat_thresh <- nodiag * (abs(nodiag) > lambda) + diag(diag(mat))
  (matHat_thresh + t(matHat_thresh)) / 2
}

matrix_threshold_perc <- function(mat, perc = 0.05) {
  mat_entries <- abs(mat[upper.tri(mat, diag = FALSE)])
  as.numeric(quantile(mat_entries, perc))
}

sta_thresholding_perc <- function(X, mat_type, mat, var_inds, perc = 0.05) {
  p <- ncol(X)
  mat_red <- mat[var_inds, var_inds]
  lambda_opt <- matrix_threshold_perc(mat_red, perc)
  mat_red_thr <- thresholding(mat_red, lambda_opt)
  mat_thr_full <- matrix(0, p, p)
  mat_thr_full[var_inds, var_inds] <- mat_red_thr
  list(X = X, mat_type = mat_type, mat = mat_thr_full, var_inds = var_inds)
}

############################################################
# 4) IPC-HD core (thr only)
############################################################
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

sta_ipchd <- function(X, mat_type, mat, var_inds, overest_type = "frac") {
  p <- ncol(X)
  mat_red <- mat[var_inds, var_inds]
  eig_mat_red <- eigen(mat_red)
  ridge_eig <- ridge_tuning(eig_mat_red$values)
  eigen_rat <- eigen_ratios(ridge_eig)
  shat <- optimal_gap(eigen_rat, overest_type, length(var_inds))$shat
  output_red <- influence_measures(eig_mat_red, shat)
  output <- rep(0, p)
  output[var_inds] <- output_red
  output
}

ipchd_thr_score <- function(X, Rhat, ipchd_perc = 0.7){
  input_cor <- list(X = X, mat_type = "cor", mat = Rhat, var_inds = 1:ncol(Rhat))
  tmp_thr <- do.call(sta_thresholding_perc, c(input_cor, list(perc = ipchd_perc)))
  do.call(sta_ipchd, c(tmp_thr, list(overest_type = "frac")))
}

############################################################
# 5) HubDetection pm generator (002): hubs -> PD precision
############################################################
.adjmat <- function(p, T0, r, ph, pnh, pneff) {
  A <- matrix(0, p, p)
  A[] <- rbinom(p * p, 1, pneff)
  A[1:T0, (r + 1):T0] <- rbinom(T0 * (T0 - r), 1, pnh)
  A[1:T0, 1:r] <- rbinom(T0 * r, 1, ph)
  A[upper.tri(A, TRUE)] <- 0
  A + t(A)
}

.rsign <- function(n = 1) 2 * rbinom(n, 1, 0.5) - 1

# Generates values in Uniform([-max,-min] U [min,max])
.rsign_unif <- function(n = 1, min = 4, max = 5) {
  runif(n, min, max) * .rsign(n)
}

.rsymmmatrix <- function(p, T0, r,
                         type = c("unif", "gaussian"),
                         hmin = 0.5, hmax = 0.8,
                         nhmin = 0.5, nhmax = 0.8,
                         neffmin = 0.5, neffmax = 0.8,
                         hsd = 1, nhsd = 1, neffsd = 1){
  type <- match.arg(type)
  theta <- matrix(0, p, p)
  
  if (type == "unif") {
    theta[] <- .rsign_unif(p * p, min = neffmin, max = neffmax)
    theta[1:T0, (r + 1):T0] <- .rsign_unif(T0 * (T0 - r), min = nhmin, max = nhmax)
    theta[1:T0, 1:r] <- .rsign_unif(T0 * r, min = hmin, max = hmax)
    theta[upper.tri(theta, TRUE)] <- 0
  } else {
    theta[] <- rnorm(p * p, sd = neffsd)
    theta[1:T0, (r + 1):T0] <- rnorm(T0 * (T0 - r), sd = nhsd)
    theta[1:T0, 1:r] <- rnorm(T0 * r, sd = hsd)
    theta[upper.tri(theta, TRUE)] <- 0
  }
  
  theta + t(theta)
}

.rhubmat <- function(p, T0, r, ph, pnh, pneff,
                     type = c("unif", "gaussian"),
                     hmin = 0.5, hmax = 0.8,
                     nhmin = 0.5, nhmax = 0.8,
                     neffmin = 0.5, neffmax = 0.8,
                     hsd = 1, nhsd = 1, neffsd = 1){
  A <- .adjmat(p, T0, r, ph, pnh, pneff)
  A * .rsymmmatrix(
    p, T0, r,
    type = type,
    hmin = hmin, hmax = hmax,
    nhmin = nhmin, nhmax = nhmax,
    neffmin = neffmin, neffmax = neffmax,
    hsd = hsd, nhsd = nhsd, neffsd = neffsd
  )
}

r.sparse.pdhubmat <- function(p, T0, r, ph, pnh, pneff,
                              diagonal_shift = 1,
                              type = c("unif", "gaussian"),
                              hmin = 0.5, hmax = 0.8,
                              nhmin = 0.5, nhmax = 0.8,
                              neffmin = 0.5, neffmax = 0.8,
                              hsd = 1, nhsd = 1, neffsd = 1){
  type <- match.arg(type)
  pm <- .rhubmat(
    p, T0, r, ph, pnh, pneff,
    type = type,
    hmin = hmin, hmax = hmax,
    nhmin = nhmin, nhmax = nhmax,
    neffmin = neffmin, neffmax = neffmax,
    hsd = hsd, nhsd = nhsd, neffsd = neffsd
  )
  lambda_min <- eigen(pm)$value[p]
  pm + (-lambda_min + diagonal_shift) * diag(p)
}

############################################################
# 6) EGG-style summary X generation
# Modified:
# - randomly sample 10% rows for UHP injection
# - randomly sample 20 traits for UHP injection
# - UHPshift default = 1
############################################################
simulate_X_EGG <- function(
    m, p, n,
    Sbb, Siguu,
    UHP = 0, UHPratio = 0.10, UHPshift = 1,
    UHP_n_traits = 20
){
  b <- MASS::mvrnorm(m, rep(0, p), Sbb) * sqrt(n) / sqrt(m)
  u <- MASS::mvrnorm(m, rep(0, p), Siguu)
  
  if (UHP == 1) {
    n_inject_rows <- max(1, round(m * UHPratio))
    idx <- sample(seq_len(m), size = n_inject_rows, replace = FALSE)
    
    n_inject_traits <- min(UHP_n_traits, p)
    jj <- sample(seq_len(p), size = n_inject_traits, replace = FALSE)
    
    b[idx, jj] <- b[idx, jj] + UHPshift / sqrt(m) * sqrt(n)
  }
  
  list(X = b + u)
}

############################################################
# 7) Robust covariance for EGG 3.1/3.2 (Spearman+MAD), Omega by w
############################################################
mad_sd <- function(x) 1.483 * median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)

spearmancov <- function(A) {
  p <- ncol(A)
  s <- apply(A, 2, mad_sd)
  R_rho <- suppressWarnings(cor(A, method = "spearman", use = "pairwise.complete.obs"))
  R_gauss <- 2 * sin(pi * R_rho / 6)
  D <- diag(s, p, p)
  S <- D %*% R_gauss %*% D
  (S + t(S)) / 2
}

estimate_omega_by_w <- function(Siguu, m, mult = 10){
  p <- ncol(Siguu)
  w <- MASS::mvrnorm(mult * m, rep(0, p), Siguu)
  cov(w)
}

estimate_genetic_cov_EGG31_32_w_author <- function(
    X, Siguu,
    m_for_w = NULL, w_mult = 10,
    min_eig = 0.001
){
  m <- nrow(X)
  if (is.null(m_for_w)) m_for_w <- m
  
  Omega_hat <- estimate_omega_by_w(Siguu = Siguu, m = m_for_w, mult = w_mult)
  Omega_hat <- (Omega_hat + t(Omega_hat)) / 2
  Omega_hat <- positveadj(Omega_hat, min.eps = min_eig)
  
  Sigma_hatbeta_rob <- spearmancov(X)
  Sigma_hatbeta_rob <- (Sigma_hatbeta_rob + t(Sigma_hatbeta_rob)) / 2
  Sigma_hatbeta_rob <- positveadj(Sigma_hatbeta_rob, min.eps = min_eig)
  
  Sigma_tilde_raw <- Sigma_hatbeta_rob - Omega_hat
  Sigma_tilde_raw <- (Sigma_tilde_raw + t(Sigma_tilde_raw)) / 2
  Sigma_tilde <- positveadj(Sigma_tilde_raw, min.eps = min_eig)
  
  list(
    Omega_hat = Omega_hat,
    Sigma_hatbeta_rob = Sigma_hatbeta_rob,
    Sigma_tilde = Sigma_tilde
  )
}

############################################################
# 8) RawInv scoring
############################################################
rawinv_score <- function(Rhat){
  Theta_hat <- MASS::ginv(Rhat)
  .alphavals(Theta_hat)
}

rawinv_degree_score <- function(Rhat, theta_perc = 0.7){
  Theta_hat <- MASS::ginv(Rhat)
  Theta_hat <- (Theta_hat + t(Theta_hat)) / 2
  
  off <- abs(Theta_hat[upper.tri(Theta_hat, diag = FALSE)])
  lam <- as.numeric(quantile(off, probs = theta_perc, na.rm = TRUE))
  
  Theta_thr <- Theta_hat
  diag_keep <- diag(Theta_thr)
  diag(Theta_thr) <- 0
  Theta_thr[abs(Theta_thr) <= lam] <- 0
  diag(Theta_thr) <- diag_keep
  Theta_thr <- (Theta_thr + t(Theta_thr)) / 2
  
  Theta_off <- Theta_thr
  diag(Theta_off) <- 0
  deg <- apply(Theta_off, 2, function(x) sum(abs(x) > 0))
  
  list(score = deg, lambda = lam)
}

############################################################
# 9) EGG 3.3 (MCP + stability selection)
############################################################
entropyloss <- function(A, B, eps = 1e-8){
  C <- A %*% B
  a <- trace(C)
  b <- Re(eigen(C, symmetric = FALSE)$values)
  ind <- which(b > eps)
  b <- sum(log(b[ind]))
  a - b - ncol(A)
}

soft <- function(a, b){
  c <- abs(a) - b
  c[c < 0] <- 0
  c * sign(a)
}

mcp <- function(a, lam, ga = 3.7){
  b <- abs(a)
  z <- soft(a, lam) / (1 - 1 / ga)
  z[b > (ga * lam)] <- a[b > (ga * lam)]
  z
}

MCPthreshold <- function(S, lam, ga = 3){
  d <- sqrt(diag(S))
  R <- cov2cor(S)
  s <- as.vector(R)
  s2 <- mcp(abs(s), lam, ga) * sign(s)
  S1 <- matrix(s2, ncol(S), ncol(S))
  diag(S1) <- 1
  diag(d) %*% S1 %*% diag(d)
}

entropy.mcp.spearman.sampling <- function(
    BETA, Rnoise,
    lamvec = seq(0.03, 0.080, length.out = 12),
    max.eps = 0.005, max.iter = 20,
    rho = 0.05, mineig = 0.01,
    subtime = 40, subfrac = 0.2, subthres = 0.98,
    alpha = 0
){
  m <- nrow(BETA)
  p <- ncol(BETA)
  alpha <- alpha * rho
  Const <- matrix(1, p, p)
  
  subvecerror <- matrix(0, length(lamvec), subtime)
  Thetalist <- array(NA_real_, c(length(lamvec), subtime, p, p))
  
  for (j in 1:subtime) {
    indsub <- sample(m, round(subfrac * m), replace = FALSE)
    BETAS  <- BETA[indsub, , drop = FALSE]
    
    S <- spearmancov(BETAS) * nrow(BETAS)
    M <- Rnoise * nrow(BETAS)
    S1 <- cov2cor1(S - M)
    
    Theta0 <- glasso::glasso(S1, 0.1)$wi
    
    BETAS2 <- BETA[-indsub, , drop = FALSE]
    S2 <- spearmancov(BETAS2) * nrow(BETAS2)
    M2 <- Rnoise * nrow(BETAS2)
    S2 <- cov2cor1(S2 - M2)
    S2 <- MCPthreshold(S2, 2 * sqrt(log(p) / m), ga = 3)
    
    for (i in seq_along(lamvec)) {
      Theta  <- Theta0
      Theta1 <- Theta0 * 0
      Delta1 <- Theta
      Gamma1 <- Delta1 * 0
      Delta2 <- Theta
      Gamma2 <- Delta2 * 0
      
      error <- norm(Theta - Theta1, "f") / sqrt(p)
      iter  <- 0
      
      while (error > max.eps && iter < max.iter) {
        Theta1 <- Theta
        
        Q <- S1 + Gamma1 - rho * Delta1
        Theta <- (-Q + matrixsqrt(Q %*% Q + 4 * rho * diag(p))$w) / (2 * rho + alpha)
        
        Delta1 <- mcp(vec(Theta + Gamma1 / rho), lamvec[i] / rho, ga = 3)
        Delta1 <- matrix(Delta1, p, p) * Const
        Gamma1 <- Gamma1 + rho * (Theta - Delta1)
        
        if (iter %% 5 == 0 && min(eigen((Theta + t(Theta)) / 2, symmetric = TRUE)$values) < mineig) {
          Q <- S1 + Gamma1 + Gamma2 - rho * (Delta1 + Delta2)
          Theta <- (-Q + matrixsqrt(Q %*% Q + 8 * rho * diag(p))$w) / (4 * rho + alpha)
          
          Delta1 <- mcp(vec(Theta + Gamma1 / rho), lamvec[i] / rho, ga = 3)
          Delta1 <- matrix(Delta1, p, p) * Const
          Gamma1 <- Gamma1 + rho * (Theta - Delta1)
          
          Delta2 <- positveadj(Theta + Gamma2 / rho, min.eps = mineig)
          Gamma2 <- Gamma2 + rho * (Theta - Delta2)
        }
        
        iter  <- iter + 1
        error <- norm(Theta - Theta1, "f") / sqrt(p)
      }
      
      df <- (sum(Delta1 != 0) - p) / 2
      subvecerror[i, j] <- entropyloss(S2, Delta1) + (log(p * (p - 1) / 2) + log(m)) / m * df
      Thetalist[i, j, , ] <- Delta1
    }
  }
  
  istar <- which.min(rowMedian(subvecerror))
  Thetalist_best <- Thetalist[istar, , , , drop = FALSE]
  Thetalist_best <- array(Thetalist_best, dim = c(subtime, p, p))
  
  K <- Thetalist_best[1, , ] * 0
  for (ii in 1:subtime) {
    K <- K + (Thetalist_best[ii, , ] != 0) / subtime
  }
  
  S <- spearmancov(BETA) * nrow(BETA)
  M <- Rnoise * nrow(BETA)
  S1 <- cov2cor1(S - M)
  
  Theta0 <- glasso::glasso(S1, 0.1)$wi
  Theta  <- Theta0
  Theta1 <- Theta0 * 0
  Delta1 <- Theta
  Gamma1 <- Delta1 * 0
  Delta2 <- Theta
  Gamma2 <- Delta2 * 0
  
  error <- norm(Theta - Theta1, "f") / sqrt(p)
  iter  <- 0
  
  while (error > max.eps && iter < (2 * max.iter)) {
    Theta1 <- Theta
    
    Q <- S1 + Gamma1 - rho * Delta1
    Theta <- (-Q + matrixsqrt(Q %*% Q + 4 * rho * diag(p))$w) / (2 * rho + alpha)
    
    Delta1 <- mcp(vec(Theta + Gamma1 / rho), lamvec[istar] / rho, ga = 3)
    Delta1 <- matrix(Delta1, p, p) * (K > subthres)
    Gamma1 <- Gamma1 + rho * (Theta - Delta1)
    
    if (iter %% 5 == 0 && min(eigen((Theta + t(Theta)) / 2, symmetric = TRUE)$values) < mineig) {
      Q <- S1 + Gamma1 + Gamma2 - rho * (Delta1 + Delta2)
      Theta <- (-Q + matrixsqrt(Q %*% Q + 8 * rho * diag(p))$w) / (4 * rho + alpha)
      
      Delta1 <- mcp(vec(Theta + Gamma1 / rho), lamvec[istar] / rho, ga = 3)
      Delta1 <- matrix(Delta1, p, p)
      Gamma1 <- Gamma1 + rho * (Theta - Delta1)
      
      Delta2 <- positveadj(Theta + Gamma2 / rho, min.eps = mineig)
      Gamma2 <- Gamma2 + rho * (Theta - Delta2)
    }
    
    iter  <- iter + 1
    error <- norm(Theta - Theta1, "f") / sqrt(p)
  }
  
  Theta_hat <- (Delta1 + t(Delta1)) / 2
  list(Theta = Theta_hat, K = K, cv.error = subvecerror, istar = istar)
}

egg33_fit_theta_from_X <- function(
    X, Omega_hat,
    lamvec = seq(0.03, 0.080, length.out = 12),
    rho = 0.05,
    subtime = 40, subfrac = 0.2, subthres = 0.98,
    mineig = 0.01, max.iter = 20, max.eps = 0.005
){
  entropy.mcp.spearman.sampling(
    BETA = X, Rnoise = Omega_hat,
    lamvec = lamvec,
    rho = rho, subtime = subtime, subfrac = subfrac, subthres = subthres,
    mineig = mineig, max.iter = max.iter, max.eps = max.eps
  )
}

############################################################
# 10) hub scores from EGG 3.3 Theta_hat
############################################################
egg_theta_alpha_score <- function(Theta_hat){
  Theta_hat <- (Theta_hat + t(Theta_hat)) / 2
  .alphavals(Theta_hat)
}

egg_theta_degree_score <- function(Theta_hat, theta_perc = 0.7){
  Theta_hat <- (Theta_hat + t(Theta_hat)) / 2
  off <- abs(Theta_hat[upper.tri(Theta_hat, diag = FALSE)])
  lam <- as.numeric(quantile(off, probs = theta_perc, na.rm = TRUE))
  
  Theta_thr <- Theta_hat
  diag_keep <- diag(Theta_thr)
  diag(Theta_thr) <- 0
  Theta_thr[abs(Theta_thr) <= lam] <- 0
  diag(Theta_thr) <- diag_keep
  Theta_thr <- (Theta_thr + t(Theta_thr)) / 2
  
  Theta_off <- Theta_thr
  diag(Theta_off) <- 0
  deg <- apply(Theta_off, 2, function(x) sum(abs(x) > 0))
  
  list(score = deg, lambda = lam)
}

############################################################
# 11) hub selection by mean + 2sd + metrics
############################################################
hub_select_mean2sd <- function(score){
  thr <- mean(score) + 2 * sd(score)
  list(thr = thr, hub_hat = (score > thr), hub_hat_idx = which(score > thr))
}

eval_hubs <- function(score, hub_true_idx){
  p <- length(score)
  truth <- rep(FALSE, p)
  truth[hub_true_idx] <- TRUE
  
  sel <- hub_select_mean2sd(score)
  pred <- sel$hub_hat
  
  tp <- sum(pred & truth)
  fp <- sum(pred & !truth)
  fn <- sum(!pred & truth)
  tn <- sum(!pred & !truth)
  
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1        <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  
  tpr <- recall
  fpr <- if ((fp + tn) == 0) 0 else fp / (fp + tn)
  
  list(
    thr = sel$thr,
    hub_hat_idx = sel$hub_hat_idx,
    n_hat = length(sel$hub_hat_idx),
    tp = tp, fp = fp, fn = fn, tn = tn,
    precision = precision, recall = recall, f1 = f1,
    tpr = tpr, fpr = fpr
  )
}

############################################################
# 12) One complete run (single replicate)
############################################################
run_one <- function(
    # ----- GRID factors -----
    p = 200,
    m = 2000,
    n = 50000,
    T0_prop = 1,
    ph = 0.8,
    UHP = 0,
    
    # ----- fixed hub count -----
    r = 5,
    
    # ----- pm generator baseline -----
    pnh = 0.05,
    pneff = 0.01,
    diagonal_shift = 2,
    type = "unif",
    hmin = 4, hmax = 5,
    nhmin = 4, nhmax = 5,
    neffmin = 4, neffmax = 5,
    
    # ----- Siguu -----
    Siguu = NULL,
    
    # ----- pleiotropy-like injection -----
    UHPratio = 0.10,
    UHPshift = 1,
    UHP_n_traits = 20,
    
    # ----- EGG 3.1/3.2 -----
    w_mult = 10,
    min_eig = 0.001,
    
    # ----- IPC-HD thr only -----
    ipchd_perc = 0.7,
    
    # ----- RawInv-degree sparsify -----
    theta_perc = 0.7,
    
    # ----- EGG 3.3 -----
    egg33_lamvec = seq(0.03, 0.080, length.out = 12),
    egg33_rho = 0.05,
    egg33_subtime = 40,
    egg33_subfrac = 0.2,
    egg33_subthres = 0.98,
    egg33_mineig = 0.01,
    egg33_max_iter = 20,
    egg33_max_eps = 0.005,
    egg_theta_perc = 0.7
){
  # (1) true pm + Sbb
  T0 <- as.integer(p * T0_prop)
  pm_true <- r.sparse.pdhubmat(
    p = p, T0 = T0, r = r, ph = ph, pnh = pnh, pneff = pneff,
    diagonal_shift = diagonal_shift,
    type = type,
    hmin = hmin, hmax = hmax,
    nhmin = nhmin, nhmax = nhmax,
    neffmin = neffmin, neffmax = neffmax
  )
  Sbb <- solve(pm_true)
  hub_true_idx <- 1:r
  
  # (2) Siguu default
  if (is.null(Siguu)) {
    Siguu <- CScov(p, 0.5) * make_Nover(p)
  }
  
  # (3) simulate X
  X <- simulate_X_EGG(
    m = m, p = p, n = n,
    Sbb = Sbb, Siguu = Siguu,
    UHP = UHP, UHPratio = UHPratio, UHPshift = UHPshift,
    UHP_n_traits = UHP_n_traits
  )$X
  
  # (4) EGG 3.1/3.2 -> Sigma_tilde -> Rhat
  est <- estimate_genetic_cov_EGG31_32_w_author(
    X = X, Siguu = Siguu, m_for_w = m, w_mult = w_mult, min_eig = min_eig
  )
  Rhat <- cov2cor(est$Sigma_tilde)
  
  # (5) scores
  score_raw <- rawinv_score(Rhat)
  
  rawdeg <- rawinv_degree_score(Rhat, theta_perc = theta_perc)
  score_raw_deg <- rawdeg$score
  
  score_ip_thr <- ipchd_thr_score(X, Rhat, ipchd_perc = ipchd_perc)
  
  egg33 <- egg33_fit_theta_from_X(
    X = X, Omega_hat = est$Omega_hat,
    lamvec = egg33_lamvec,
    rho = egg33_rho,
    subtime = egg33_subtime,
    subfrac = egg33_subfrac,
    subthres = egg33_subthres,
    mineig = egg33_mineig,
    max.iter = egg33_max_iter,
    max.eps = egg33_max_eps
  )
  Theta_hat <- egg33$Theta
  
  score_egg_alpha <- egg_theta_alpha_score(Theta_hat)
  eggdeg <- egg_theta_degree_score(Theta_hat, theta_perc = egg_theta_perc)
  score_egg_deg <- eggdeg$score
  
  # (6) eval
  list(
    eval = list(
      raw = eval_hubs(score_raw, hub_true_idx),
      raw_deg = eval_hubs(score_raw_deg, hub_true_idx),
      ip_thr = eval_hubs(score_ip_thr, hub_true_idx),
      egg_theta_alpha = eval_hubs(score_egg_alpha, hub_true_idx),
      egg_theta_deg = eval_hubs(score_egg_deg, hub_true_idx)
    )
  )
}

############################################################
# 13) Repeat R times for ONE combo
############################################################
run_reps <- function(R = 20, seed_base = 1, ...){
  methods <- c("raw", "raw_deg", "ip_thr", "egg_theta_alpha", "egg_theta_deg")
  
  out <- vector("list", R)
  for (i in seq_len(R)) {
    set.seed(seed_base + i - 1)
    out[[i]] <- run_one(...)
  }
  
  runs_long <- do.call(rbind, lapply(seq_len(R), function(i){
    z <- out[[i]]
    do.call(rbind, lapply(methods, function(m){
      e <- z$eval[[m]]
      data.frame(
        rep = i,
        seed = seed_base + i - 1,
        method = m,
        thr = e$thr,
        n_hat = e$n_hat,
        tp = e$tp, fp = e$fp, fn = e$fn, tn = e$tn,
        tpr = e$tpr, fpr = e$fpr,
        precision = e$precision, recall = e$recall, f1 = e$f1,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  summarize_method <- function(df_m){
    data.frame(
      tpr_mean = mean(df_m$tpr), tpr_sd = sd(df_m$tpr),
      fpr_mean = mean(df_m$fpr), fpr_sd = sd(df_m$fpr),
      precision_mean = mean(df_m$precision), precision_sd = sd(df_m$precision),
      recall_mean = mean(df_m$recall), recall_sd = sd(df_m$recall),
      f1_mean = mean(df_m$f1), f1_sd = sd(df_m$f1),
      n_hat_mean = mean(df_m$n_hat), n_hat_sd = sd(df_m$n_hat),
      stringsAsFactors = FALSE
    )
  }
  
  summary <- do.call(rbind, lapply(methods, function(m){
    df_m <- runs_long[runs_long$method == m, , drop = FALSE]
    cbind(method = m, summarize_method(df_m))
  }))
  
  list(runs_long = runs_long, summary = summary)
}

############################################################
# 14) GRID definition (216 combos)
############################################################
make_grid <- function(){
  expand.grid(
    T0_prop = c(0.5, 0.75, 1),
    n       = c(30000, 50000, 100000),
    m       = c(500, 1000, 2000),
    UHP     = c(0, 1),
    ph      = c(0.4, 0.8),
    p       = c(100, 200),
    stringsAsFactors = FALSE
  )
}

############################################################
# 15) PARALLEL GRID runner (Linux)
############################################################
run_grid_parallel <- function(
    R = 20,
    seed0 = 123,
    
    # parallel
    ncores = 48,
    
    # output
    out_dir = "results_grid",
    prefix = "hub_egg_grid216_R20",
    
    # fixed params (not in grid)
    ipchd_perc = 0.7,
    theta_perc = 0.7,
    egg_theta_perc = 0.7,
    
    UHPratio = 0.10,
    UHPshift = 1,
    UHP_n_traits = 20,
    
    egg33_subtime = 15,
    egg33_subfrac = 0.2,
    egg33_subthres = 0.98,
    egg33_lamvec = seq(0.03, 0.080, length.out = 12),
    
    diagonal_shift = 2,
    pnh = 0.05,
    pneff = 0.01,
    
    verbose = TRUE
){
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  grid <- make_grid()
  n_combo <- nrow(grid)
  
  if (verbose) {
    cat(sprintf("Grid combos = %d; R = %d; ncores = %d\n", n_combo, R, ncores))
    cat("Output folder: ", normalizePath(out_dir), "\n", sep = "")
  }
  
  run_one_combo <- function(g){
    pars <- grid[g, ]
    
    # deterministic, non-overlapping seeds per combo
    seed_base <- seed0 + (g - 1) * 100000
    
    rr <- run_reps(
      R = R,
      seed_base = seed_base,
      
      # GRID factors
      p = pars$p,
      ph = pars$ph,
      T0_prop = pars$T0_prop,
      n = pars$n,
      m = pars$m,
      UHP = pars$UHP,
      
      # fixed
      diagonal_shift = diagonal_shift,
      pnh = pnh,
      pneff = pneff,
      
      ipchd_perc = ipchd_perc,
      theta_perc = theta_perc,
      egg_theta_perc = egg_theta_perc,
      
      UHPratio = UHPratio,
      UHPshift = UHPshift,
      UHP_n_traits = UHP_n_traits,
      
      egg33_subtime = egg33_subtime,
      egg33_subfrac = egg33_subfrac,
      egg33_subthres = egg33_subthres,
      egg33_lamvec = egg33_lamvec
    )
    
    # annotate runs_long
    rl <- rr$runs_long
    rl$p <- pars$p
    rl$ph <- pars$ph
    rl$T0_prop <- pars$T0_prop
    rl$n <- pars$n
    rl$m <- pars$m
    rl$UHP <- pars$UHP
    rl$combo_id <- g
    rl$R <- R
    
    # annotate summary
    sm <- rr$summary
    sm$p <- pars$p
    sm$ph <- pars$ph
    sm$T0_prop <- pars$T0_prop
    sm$n <- pars$n
    sm$m <- pars$m
    sm$UHP <- pars$UHP
    sm$combo_id <- g
    sm$R <- R
    
    list(runs_long = rl, summary = sm)
  }
  
  # parallel over combos
  res_list <- mclapply(
    X = seq_len(n_combo),
    FUN = function(g){
      if (verbose && (g %% 10 == 0)) {
        cat(sprintf("[worker] starting combo %d/%d\n", g, n_combo))
      }
      run_one_combo(g)
    },
    mc.cores = ncores
  )
  
  runs_long <- do.call(rbind, lapply(res_list, `[[`, "runs_long"))
  summary_grid <- do.call(rbind, lapply(res_list, `[[`, "summary"))
  
  # save CSVs
  f_runs <- file.path(out_dir, paste0(prefix, "_runs_long.csv"))
  f_sum  <- file.path(out_dir, paste0(prefix, "_summary_grid.csv"))
  
  write.csv(runs_long, f_runs, row.names = FALSE)
  write.csv(summary_grid, f_sum, row.names = FALSE)
  
  cat("\nSaved:\n", f_runs, "\n", f_sum, "\n", sep = "")
  
  invisible(list(grid = grid, runs_long = runs_long, summary_grid = summary_grid))
}


# ---- Execute the simulation study ----
if (run_analysis) {
# Full 216-scenario grid. This is a substantial computation, not a quick demo.
ncores <- as.integer(Sys.getenv("HINGE_CORES", "1"))
if (is.na(ncores) || ncores < 1L) stop("HINGE_CORES must be a positive integer.")
if (.Platform$OS.type == "windows" && ncores > 1L) stop("Use one core on Windows; fork parallelism requires Linux/macOS.")
dir.create("results/simulation", recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = "results/simulation/sessionInfo.txt")
run_grid_parallel(R = 100, seed0 = 123, ncores = ncores,
  out_dir = "results/simulation", prefix = "hub_egg_grid216_R100",
  egg33_subtime = 15)

}

# ---- Draw F1, TPR and precision figures ----
if (run_plots) {
# Adapted from the supplied manuscript plotting script; run from repository root.
input_file <- "results/simulation/hub_egg_grid216_R100_summary_grid.csv"
plot_dir <- "results/figures/simulation"
if (!file.exists(input_file)) stop("Run the analysis section of simulation.R first.")
df <- read.csv(input_file)
required <- c("p", "n", "m", "UHP", "ph", "T0_prop", "method", "f1_mean", "tpr_mean", "precision_mean")
if (!all(required %in% names(df))) stop("Missing required simulation summary columns.")
# These publication layouts require the full grid, not the smoke-test summary.
expected <- expand.grid(p=c(100,200), n=c(30000,50000,100000), m=c(500,1000,2000), UHP=c(0,1), ph=c(0.4,0.8), T0_prop=c(0.5,0.75,1), method=c("raw","raw_deg","ip_thr","egg_theta_alpha","egg_theta_deg"), stringsAsFactors=FALSE)
keys <- names(expected)
if (anyDuplicated(df[keys])) stop("Duplicate scenario-method rows in summary.")
if (nrow(merge(expected, df, by=keys)) != nrow(expected)) stop("Plotting requires all 216 scenarios and all five methods.")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(grid)
library(ggh4x)
#F1#############
# display.brewer.all()
# 图形参数设置

method_levels <- c("Raw_alpha", "Raw_deg", "EGG_alpha", "EGG_deg", "HINGE")

col_map <- c(
  Raw_alpha = "#FB8072",  # 深橙红
  Raw_deg   = "#BEBADA",  # 深紫
  EGG_alpha = "#FDB462",  # 金橙
  EGG_deg   = "#B3DE69",  # 深青绿
  HINGE     = "#80B1D3"   # 深蓝
)

lt_map <- c(
  Raw_alpha = "solid",
  Raw_deg   = "dashed",
  EGG_alpha = "dotdash",
  EGG_deg   = "solid",
  HINGE     = "longdash"
)

sh_map <- c(
  Raw_alpha = 16,
  Raw_deg   = 4,
  EGG_alpha = 17,
  EGG_deg   = 15,
  HINGE     = 8
)

# 绘图函数
make_f1_plot <- function(df, p_fix, n_fix, UHP_fix, subtitle_text, show_legend = FALSE) {
  
  df_plot <- df %>%
    filter(
      p == p_fix,
      as.numeric(n) == n_fix,
      UHP == UHP_fix,
      T0_prop %in% c(0.5, 0.75, 1)
    ) %>%
    mutate(
      T = as.integer(T0_prop * p),
      m_num = as.numeric(m)
    ) %>%
    group_by(ph, T, m_num, method) %>%
    summarise(f1_mean = mean(f1_mean), .groups = "drop") %>%
    mutate(
      ph = factor(ph),
      T  = factor(T, levels = sort(unique(T)))
    ) %>%
    mutate(
      method = recode(
        method,
        raw = "Raw_alpha",
        ip_thr = "HINGE",
        egg_theta_deg = "EGG_deg",
        egg_theta_alpha = "EGG_alpha",
        raw_deg = "Raw_deg"
      )
    )
  
  present_levels <- method_levels[method_levels %in% unique(as.character(df_plot$method))]
  df_plot$method <- factor(df_plot$method, levels = present_levels)
  df_plot$method <- droplevels(df_plot$method)
  
  p <- ggplot(
    df_plot,
    aes(
      x = m_num, y = f1_mean,
      color = method, linetype = method, shape = method,
      group = method
    )
  ) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.0, stroke = 0.7) +
    
    facet_grid(
      ph ~ T,
      labeller = labeller(
        T  = function(x) paste0("T = ", x),
        ph = function(x) paste0("Ph = ", x)
      )
    ) +
    
    scale_x_continuous(
      breaks = c(500, 1000, 2000),
      labels = c("500", "1000", "2000"),
      minor_breaks = seq(500, 2000, by = 250)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      minor_breaks = seq(0, 1, by = 0.05)
    ) +
    
    scale_color_manual(values = col_map, breaks = present_levels) +
    scale_linetype_manual(values = lt_map, breaks = present_levels) +
    scale_shape_manual(values = sh_map, breaks = present_levels) +
    
    labs(
      title = "Hub Detection Performance (F1-score)",
      subtitle = subtitle_text,
      x = "Number of SNPs m",
      y = "F1-score",
      color = "Method"
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      
      panel.grid.major = element_line(color = "grey90", linewidth = 0.45),
      panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
      aspect.ratio = 1,
      
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = unit(1.6, "cm"),
      legend.key.height = unit(0.9, "cm")
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(
          linetype  = lt_map[present_levels],
          shape     = sh_map[present_levels],
          linewidth = 0.9,
          size      = 3.0
        )
      ),
      linetype = "none",
      shape = "none"
    )
  
  return(p)
}

# 绘图
# 参数网格
p_vals <- c(100, 200)
n_vals <- c(30000, 50000, 100000)

# 只提取一次图例（选择任意一组参数即可，legend 与具体图形内容无关）
p_legend <- make_f1_plot(
  df = df,
  p_fix = 200,
  n_fix = 30000,
  UHP_fix = 1,
  subtitle_text = "",
  show_legend = TRUE
)
legend <- cowplot::get_legend(p_legend)

# 遍历所有 p 与 n 的组合
for (p_fix in p_vals) {
  for (n_fix in n_vals) {
    
    # 左侧：UHP = 0
    p_left <- make_f1_plot(
      df = df,
      p_fix = p_fix,
      n_fix = n_fix,
      UHP_fix = 0,
      subtitle_text = paste0("Fixed settings: p = ", p_fix, ", n = ", n_fix, ", UHP = 0"),
      show_legend = FALSE
    )
    
    # 右侧：UHP = 1
    p_right <- make_f1_plot(
      df = df,
      p_fix = p_fix,
      n_fix = n_fix,
      UHP_fix = 1,
      subtitle_text = paste0("Fixed settings: p = ", p_fix, ", n = ", n_fix, ", UHP = 1"),
      show_legend = FALSE
    )
    
    # 横向拼接左右两图
    p_main <- cowplot::plot_grid(
      p_left, p_right,
      nrow = 1,
      align = "h",
      axis = "tb",
      rel_widths = c(1, 1)
    )
    
    # 加上右侧图例
    p_final <- cowplot::plot_grid(
      p_main, legend,
      nrow = 1,
      rel_widths = c(1, 0.12)
    )
    
    # 显示（可选）
    # Results are saved below; do not open an extra graphics device.
    
    # 动态文件名保存
    base_name <- paste0("F1_m_fixed_p", p_fix, "_n", n_fix, "_UHP01_side_by_side_final")
    ggsave(
      filename = file.path(plot_dir, paste0(base_name, ".png")),
      plot = p_final,
      width = 24,
      height = 8,
      dpi = 300
    )
    ggsave(
      filename = file.path(plot_dir, paste0(base_name, ".pdf")),
      plot = p_final,
      width = 24,
      height = 8
    )
  }
}

#####TPR#####
# 参数网格
p_vals   <- c(100, 200)
ph_vals  <- c(0.4, 0.8)
T0_vals  <- c(0.5, 0.75, 1)

# 固定不变的部分（方法与配色）
method_levels <- c("raw", "raw_deg", "egg_theta_alpha", "egg_theta_deg", "ip_thr")
method_labels <- c(
  raw             = "Raw_alpha",
  raw_deg         = "Raw_deg",
  egg_theta_alpha = "EGG_alpha",
  egg_theta_deg   = "EGG_deg",
  ip_thr          = "HINGE"
)
fill_map <- c(
  raw             = "#FB8072",
  raw_deg         = "#BEBADA",
  egg_theta_alpha = "#FDB462",
  egg_theta_deg   = "#B3DE69",
  ip_thr          = "#80B1D3"
)

# 遍历所有参数组合
for (fixed_p in p_vals) {
  for (fixed_ph in ph_vals) {
    for (fixed_T0 in T0_vals) {
      
      # 数据筛选 + 预处理（与原来逻辑完全相同）
      dat <- df %>%
        filter(ph == fixed_ph, p == fixed_p, T0_prop == fixed_T0) %>%
        mutate(
          pleio = ifelse(UHP == 1, "10% pleiotropy", "no pleiotropy"),
          n_num = as.numeric(n),
          n_lab = case_when(
            n_num == 30000  ~ "n = 30,000",
            n_num == 50000  ~ "n = 50,000",
            n_num == 100000 ~ "n = 100,000",
            TRUE ~ paste0("n = ", format(n_num, big.mark = ",", scientific = FALSE))
          ),
          tpr_mean_plot = dplyr::coalesce(tpr_mean, recall_mean),
          pleio_f = factor(pleio, levels = c("no pleiotropy", "10% pleiotropy")),
          n_f     = factor(n_lab, levels = c("n = 30,000", "n = 50,000", "n = 100,000")),
          m       = factor(m, levels = sort(unique(m)))
        ) %>%
        mutate(method = factor(method, levels = method_levels))
      
      # 动态 subtitle
      sub_text <- paste0(
        "Fixed settings: ph = ", fixed_ph,
        ", p = ", fixed_p,
        ", T/p = ", fixed_T0
      )
      
      # 作图
      p <- ggplot(dat, aes(x = method, y = tpr_mean_plot, fill = method)) +
        geom_col(width = 0.75, color = "grey35", linewidth = 0.25) +
        coord_flip() +
        facet_nested(
          m ~ pleio_f + n_f,
          labeller = labeller(m = function(x) paste0("m = ", x))
        ) +
        scale_fill_manual(values = fill_map, breaks = method_levels, labels = method_labels) +
        scale_x_discrete(labels = method_labels) +
        scale_y_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.25),
          labels = sprintf("%.2f", seq(0, 1, by = 0.25)),
          expand = expansion(mult = c(0, 0.05))
        ) +
        labs(
          title = "Hub Detection Performance (TPR)",
          subtitle = sub_text,
          x = NULL,
          y = "True Positive Rate (TPR)"
        ) +
        theme_bw(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 13),
          axis.text.x = element_text(size = 10),
          strip.placement = "outside",
          strip.background.x = element_rect(fill = "grey90", color = "grey40", linewidth = 0.6),
          strip.background.y = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
          strip.text.x = element_text(size = 11),
          strip.text.y.right = element_text(angle = -90, face = "bold", size = 10),
          strip.nestline = element_line(color = "grey20", linewidth = 0.6),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey88", linewidth = 0.4),
          panel.spacing.x = unit(0.9, "lines"),
          panel.spacing.y = unit(0.5, "lines"),
          legend.position = "none"
        )
      
      # 动态文件名：将小数转为整数编码，避免小数点
      # ph 保留一位小数，乘以10；T0 保留两位小数，乘以100
      ph_code  <- sprintf("%02d", as.integer(fixed_ph * 10))   # e.g., 04, 08
      T0_code  <- sprintf("%03d", as.integer(fixed_T0 * 100)) # e.g., 050, 075, 100
      base_name <- sprintf("TPR_p%d_ph%s_T0%s", fixed_p, ph_code, T0_code)
      
      # 保存 png 和 pdf
      ggsave(
        filename = file.path(plot_dir, paste0(base_name, ".png")),
        plot = p,
        width = 13,
        height = 8
      )
      ggsave(
        filename = file.path(plot_dir, paste0(base_name, ".pdf")),
        plot = p,
        width = 13,
        height = 8
      )
      
      # 如果需要预览图，可以取消下一行的注释（注意会弹出很多窗口）
      # print(p)
    }
  }
}
#####precision#####
# display.brewer.all()
# 图形参数设置

method_levels <- c("Raw_alpha", "Raw_deg", "EGG_alpha", "EGG_deg", "HINGE")

col_map <- c(
  Raw_alpha = "#FB8072",  # 深橙红
  Raw_deg   = "#BEBADA",  # 深紫
  EGG_alpha = "#FDB462",  # 金橙
  EGG_deg   = "#B3DE69",  # 深青绿
  HINGE     = "#80B1D3"   # 深蓝
)

lt_map <- c(
  Raw_alpha = "solid",
  Raw_deg   = "dashed",
  EGG_alpha = "dotdash",
  EGG_deg   = "solid",
  HINGE     = "longdash"
)

sh_map <- c(
  Raw_alpha = 16,
  Raw_deg   = 4,
  EGG_alpha = 17,
  EGG_deg   = 15,
  HINGE     = 8
)

# 绘图函数
make_precision_plot <- function(df, p_fix, n_fix, UHP_fix, subtitle_text, show_legend = FALSE) {
  
  df_plot <- df %>%
    filter(
      p == p_fix,
      as.numeric(n) == n_fix,
      UHP == UHP_fix,
      T0_prop %in% c(0.5, 0.75, 1)
    ) %>%
    mutate(
      T = as.integer(T0_prop * p),
      m_num = as.numeric(m)
    ) %>%
    group_by(ph, T, m_num, method) %>%
    summarise(precision_mean = mean(precision_mean), .groups = "drop") %>%
    mutate(
      ph = factor(ph),
      T  = factor(T, levels = sort(unique(T)))
    ) %>%
    mutate(
      method = recode(
        method,
        raw = "Raw_alpha",
        ip_thr = "HINGE",
        egg_theta_deg = "EGG_deg",
        egg_theta_alpha = "EGG_alpha",
        raw_deg = "Raw_deg"
      )
    )
  
  present_levels <- method_levels[method_levels %in% unique(as.character(df_plot$method))]
  df_plot$method <- factor(df_plot$method, levels = present_levels)
  df_plot$method <- droplevels(df_plot$method)
  
  p <- ggplot(
    df_plot,
    aes(
      x = m_num, y = precision_mean,
      color = method, linetype = method, shape = method,
      group = method
    )
  ) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 2.0, stroke = 0.7) +
    
    facet_grid(
      ph ~ T,
      labeller = labeller(
        T  = function(x) paste0("T = ", x),
        ph = function(x) paste0("Ph = ", x)
      )
    ) +
    
    scale_x_continuous(
      breaks = c(500, 1000, 2000),
      labels = c("500", "1000", "2000"),
      minor_breaks = seq(500, 2000, by = 250)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      minor_breaks = seq(0, 1, by = 0.05)
    ) +
    
    scale_color_manual(values = col_map, breaks = present_levels) +
    scale_linetype_manual(values = lt_map, breaks = present_levels) +
    scale_shape_manual(values = sh_map, breaks = present_levels) +
    
    labs(
      title = "Hub Detection Performance (Precision)",
      subtitle = subtitle_text,
      x = "Number of SNPs m",
      y = "Precision",
      color = "Method"
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      
      panel.grid.major = element_line(color = "grey90", linewidth = 0.45),
      panel.grid.minor = element_line(color = "grey85", linewidth = 0.25),
      panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "grey85", color = "grey30", linewidth = 0.6),
      aspect.ratio = 1,
      
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 13),
      legend.key.width = unit(1.6, "cm"),
      legend.key.height = unit(0.9, "cm")
    ) +
    
    guides(
      color = guide_legend(
        override.aes = list(
          linetype  = lt_map[present_levels],
          shape     = sh_map[present_levels],
          linewidth = 0.9,
          size      = 3.0
        )
      ),
      linetype = "none",
      shape = "none"
    )
  
  return(p)
}

# 绘图
# 参数网格
p_vals <- c(100, 200)
n_vals <- c(30000, 50000, 100000)

# 只提取一次图例（选择任意一组参数即可，legend 与具体图形内容无关）
p_legend <- make_precision_plot(
  df = df,
  p_fix = 200,
  n_fix = 30000,
  UHP_fix = 1,
  subtitle_text = "",
  show_legend = TRUE
)
legend <- cowplot::get_legend(p_legend)

# 遍历所有 p 与 n 的组合
for (p_fix in p_vals) {
  for (n_fix in n_vals) {
    
    # 左侧：UHP = 0
    p_left <- make_precision_plot(
      df = df,
      p_fix = p_fix,
      n_fix = n_fix,
      UHP_fix = 0,
      subtitle_text = paste0("Fixed settings: p = ", p_fix, ", n = ", n_fix, ", UHP = 0"),
      show_legend = FALSE
    )
    
    # 右侧：UHP = 1
    p_right <- make_precision_plot(
      df = df,
      p_fix = p_fix,
      n_fix = n_fix,
      UHP_fix = 1,
      subtitle_text = paste0("Fixed settings: p = ", p_fix, ", n = ", n_fix, ", UHP = 1"),
      show_legend = FALSE
    )
    
    # 横向拼接左右两图
    p_main <- cowplot::plot_grid(
      p_left, p_right,
      nrow = 1,
      align = "h",
      axis = "tb",
      rel_widths = c(1, 1)
    )
    
    # 加上右侧图例
    p_final <- cowplot::plot_grid(
      p_main, legend,
      nrow = 1,
      rel_widths = c(1, 0.12)
    )
    
    # 显示（可选）
    # Results are saved below; do not open an extra graphics device.
    
    # 动态文件名保存
    base_name <- paste0("Precision_m_fixed_p", p_fix, "_n", n_fix, "_UHP01_side_by_side_final")
    ggsave(
      filename = file.path(plot_dir, paste0(base_name, ".png")),
      plot = p_final,
      width = 24,
      height = 8,
      dpi = 300
    )
    ggsave(
      filename = file.path(plot_dir, paste0(base_name, ".pdf")),
      plot = p_final,
      width = 24,
      height = 8
    )
  }
}

capture.output(sessionInfo(), file = file.path(plot_dir, "sessionInfo.txt"))


}
