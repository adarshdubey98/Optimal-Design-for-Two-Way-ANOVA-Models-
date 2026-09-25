# =============================================================================
# IUT SUCCESSIVE COMPARISON DESIGN (NO2 TIERS x LOCATION TYPE)
#
# Computational structure: 
# 1. Fixed cell means and pooled SD (Table 4)
# 2. Exact two-sided IUT power evaluation via multivariate normal
# 3. Max-min design optimization (1D line search)
# 4. Integer allocation (Hamilton's method)
# 5. Monte Carlo empirical power & sample size search
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages(library(mvtnorm))

# =============================================================================
# 1. FIXED INPUT & SUMMARY STATISTICS
# =============================================================================

M_sampled <- matrix(c(
  121.2296, 161.3533,  84.7333,   # Tier 1: Low NO2
  155.9911, 161.4000, 155.1333,   # Tier 2: Med-Low NO2
  180.2076, 163.7333, 231.7333,   # Tier 3: Med-High NO2
  185.6422, 210.0333, 303.0000    # Tier 4: High NO2
), nrow = 4, ncol = 3, byrow = TRUE)

rownames(M_sampled) <- c("Tier 1: Low NO2", "Tier 2: Med-Low NO2", "Tier 3: Med-High NO2", "Tier 4: High NO2")
colnames(M_sampled) <- c("Residential", "Industrial", "Sensitive")

pooled_sd <- 68.08
K <- nrow(M_sampled); C <- ncol(M_sampled); m <- K - 1;n_sim <- 100000
N_total <- 180; alpha <- 0.05; c_alpha <- qnorm(1 - alpha / 2) # No multiplicity correction for IUT

row_means <- rowMeans(M_sampled)
adj_diffs <- diff(row_means)
min_delta <- min(abs(adj_diffs))

cat("===============================================================================\nFIXED INPUT SUMMARY\n===============================================================================\n")
cat("\nCell-mean matrix (PM10, µg/m³):\n"); print(round(M_sampled, 3))
cat(sprintf("\nRow means            : %s\n", paste(sprintf("%.3f", row_means), collapse = "  ")))
cat(sprintf("Adjacent differences : %s\n", paste(sprintf("%.3f", adj_diffs), collapse = "  ")))
cat(sprintf("Min adjacent gap (δ) : %.3f\nPooled SD (σ̂)        : %.3f\nδ / σ̂                : %.4f\n", min_delta, pooled_sd, min_delta / pooled_sd))

# =============================================================================
# 2. DESIGN CLASS & LFC POWER EVALUATION
# =============================================================================

make_design <- function(w) {
  w <- as.numeric(w[1])
  if (length(w) != 1L || !is.finite(w) || w <= 0 || w >= 0.5) stop("w must lie strictly between 0 and 0.5.")
  p <- c(w, 0.5 - w, 0.5 - w, w)
  stopifnot(all(p > 0), abs(sum(p) - 1) < 1e-12)
  return(p)
}

Sigma_successive <- function(p) {
  Sig <- diag(m)
  for (i in 1:m) {
    for (j in 1:m) {
      if (abs(i - j) == 1) {
        mid <- max(i, j)
        Sig[i, j] <- -(1 / p[mid]) / sqrt((1/p[i] + 1/p[i+1]) * (1/p[j] + 1/p[j+1]))
      }
    }
  }
  return(Sig)
}

power_IUT_successive <- function(w, delta, N_reference = N_total) {
  if (w <= 0 || w >= 0.5) return(0)
  p <- make_design(w); Sigma <- Sigma_successive(p)
  
  noncentrality <- sapply(1:m, function(i) {
    standard_error <- (pooled_sd / sqrt(N_reference)) * sqrt(1 / p[i] + 1 / p[i + 1])
    delta / standard_error
  })
  
  sign_patterns <- as.matrix(expand.grid(rep(list(c(-1, 1)), m)))
  orthant_probabilities <- apply(sign_patterns, 1, function(sv) {
    lower <- ifelse(sv > 0, c_alpha, -Inf); upper <- ifelse(sv > 0, Inf, -c_alpha)
    as.numeric(pmvnorm(lower = lower, upper = upper, mean = noncentrality, sigma = Sigma, algorithm = Miwa(steps = 512)))
  })
  
  min(max(sum(orthant_probabilities), 0), 1)
}

optimize_w_iut <- function(delta, n_coarse = 80) {
  w_grid <- seq(0.01, 0.49, length.out = n_coarse)
  vals <- sapply(w_grid, function(w) power_IUT_successive(w, delta))
  best <- which.max(vals)
  
  lo <- w_grid[max(1, best - 1)]; hi <- w_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(w) power_IUT_successive(w, delta), interval = c(lo, hi), maximum = TRUE)
  
  if (opt$objective >= vals[best]) list(w = opt$maximum, power = opt$objective) else list(w = w_grid[best], power = vals[best])
}

cat("\n===============================================================================\nMAX-MIN DESIGN OPTIMIZATION\n===============================================================================\n")
res <- optimize_w_iut(min_delta)
w_opt <- res$w; p_opt <- make_design(w_opt); p_bal <- rep(1 / K, K)

cat(sprintf("Optimal w (w_opt)       : %.6f\nLFC power at w_opt      : %.6f\nRow proportions p_opt   : %s\n", w_opt, res$power, paste(sprintf("%.6f", p_opt), collapse = "  ")))

# =============================================================================
# 3. ROBUST INTEGER ALLOCATION
# =============================================================================

allocate_integers <- function(p, total_N) {
  # 1. Decimal targets calculate karein
  target <- outer(p, rep(1 / C, C)) * total_N
  
  # 2. Base allocation (sirf integer part lein)
  N_mat <- floor(target)
  
  # 3. Bache hue observations (remainder) nikalen
  rem <- total_N - sum(N_mat)
  
  # 4. Agar kuch bacha hai, toh sabse bade fractions walo ko de dein
  if (rem > 0) {
    frac <- as.vector(target - floor(target))
    idx <- order(frac, decreasing = TRUE, method = "radix")[seq_len(rem)]
    N_mat[idx] <- N_mat[idx] + 1L
  }
  
  dimnames(N_mat) <- dimnames(M_sampled)
  return(N_mat)
}
N_opt <- allocate_integers(p_opt, N_total)
N_bal <- allocate_integers(p_bal, N_total)

cat("\n===============================================================================\nEXACT INTEGER DESIGNS AT N =", N_total, "\n===============================================================================\n")
cat("\nMax-min (cell-wise rounding):\n"); print(N_opt)
cat(sprintf("  Row totals: %s  (sum = %d)\n\nBalanced (cell-wise rounding):\n", paste(rowSums(N_opt), collapse = "  "), sum(N_opt)))
print(N_bal)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = "  "), sum(N_bal)))

# =============================================================================
# 4. EMPIRICAL POWER SIMULATION
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L; row_vars <- rowSums(sd^2 / N_mat) / (C^2)
  diff_se <- sqrt(row_vars[1:(K - 1)] + row_vars[2:K])
  
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    for (i in 1:K) for (j in 1:C) Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
    Z <- diff(rowMeans(Ybar)) / diff_se
    if (all(abs(Z) > c_alpha)) ctr <- ctr + 1L
  }
  return(ctr / n_sim)
}


cat("\n===============================================================================\nEMPIRICAL POWER AT N =", N_total, "\n===============================================================================\n")
cat(sprintf("Critical value : %.4f (no multiplicity correction)\nSimulation size: %d\n\n", c_alpha, n_sim))

power_opt <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim)
power_bal <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim)
se_opt <- sqrt(power_opt * (1 - power_opt) / n_sim); se_bal <- sqrt(power_bal * (1 - power_bal) / n_sim)
gain <- 100 * (power_opt - power_bal) / power_bal

cat(sprintf("Max-min design  : %.4f (SE %.4f)\nBalanced design : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", power_opt, se_opt, power_bal, se_bal, gain))

# =============================================================================
# 5. SAMPLE SIZE SEARCH FOR 90% POWER
# =============================================================================

find_n90 <- function(p, M, sd, c_alpha, target = 0.90, N_min = 20, N_max = 3000) {
  N_val <- N_min; pow <- 0
  cat("  [Coarse step 10] ... ")
  while (pow < target && N_val < N_max) { N_val <- N_val + 10; pow <- simulate_power(allocate_integers(p, N_val), M, sd, c_alpha, 500) }
  cat(sprintf("reached ~%d.  [Fine step 1] ... ", N_val))
  
  N_val <- max(N_min, N_val - 15); pow <- 0
  while (pow < target && N_val < N_max) { N_val <- N_val + 1; pow <- simulate_power(allocate_integers(p, N_val), M, sd, c_alpha, 2000) }
  cat(sprintf("N = %d (power %.2f%%)\n", N_val, 100 * pow))
  list(N = N_val, power = pow)
}

cat("\n===============================================================================\nSAMPLE SIZE FOR 90% POWER\n===============================================================================\n")
cat("> Max-min design:\n"); n90_opt <- find_n90(p_opt, M_sampled, pooled_sd, c_alpha)
cat("> Balanced design:\n"); n90_bal <- find_n90(p_bal, M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 6. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (IUT Successive Comparison) ══════════════════════\n")
cat(sprintf(" Min adjacent gap (delta)  : %.3f\n Pooled SD (sigma)         : %.3f\n", min_delta, pooled_sd))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" Power @ N=%d  Max-min      : %.4f\n Power @ N=%d  Balanced     : %.4f\n Relative gain             : %+.2f%%\n", N_total, power_opt, N_total, power_bal, gain))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" 90%% N  Max-min            : %d (power %.4f)\n 90%% N  Balanced           : %d (power %.4f)\n Sample size savings       : %d (%.2f%% reduction)\n", n90_opt$N, n90_opt$power, n90_bal$N, n90_bal$power, n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N))
cat("═══════════════════════════════════════════════════════════════════════════════════════\n")