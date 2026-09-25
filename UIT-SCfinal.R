# =============================================================================
# MAX-MIN DESIGN FOR THE UIT SUCCESSIVE-COMPARISON TEST (K = 3)
#
# Properties used for optimization:
# (P1) Reflection symmetry of design on Ξ* = { p : p_1 = p_3, p(w) = (w, 1-2w, w) }
# (P2) Mirror LFCs have equal power on Ξ*
# (P3) Symmetric-prior restriction (q1 = q2 = 0.5)
# (P4) H-algorithm collapses to a single 1D optimum-on-the-average problem
# (P5) Deterministic integration (Miwa) for exact numerical precision
# (P6) Robust Hamilton's method for cell-wise integer rounding
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({ library(mvtnorm); library(nloptr) })

# =============================================================================
# 1. FIXED INPUT & GLOBAL SETUP
# =============================================================================

M_sampled <- matrix(c(
  71.600, 142.600, 212.90,    # Tier 1: Low SO2
  168.000, 140.247, 166.80,    # Tier 2: Mid SO2
  168.782, 193.200, 134.95     # Tier 3: High SO2
), nrow = 3, ncol = 3, byrow = TRUE)

rownames(M_sampled) <- c("Tier 1: Low SO2", "Tier 2: Mid SO2", "Tier 3: High SO2")
colnames(M_sampled) <- c("Monsoon", "Summer", "Winter")

pooled_sd <- 65.533
K <- nrow(M_sampled); C <- ncol(M_sampled); N_total <- 90;n_sim <- 100000
alpha <- 0.05; c_alpha <- qnorm(1 - alpha / (2 * (K - 1))) # Bonferroni correction for UIT

row_means <- rowMeans(M_sampled); max_delta <- max(abs(diff(row_means)))

cat("===============================================================================\nFIXED INPUT SUMMARY (UIT SUCCESSIVE K = 3)\n===============================================================================\n")
cat("\nCell mean matrix (PM10, µg/m³):\n"); print(round(M_sampled, 3))
cat(sprintf("\nRow means       : %s\nAdjacent diffs  : %s\nMax gap (delta) : %.3f\nPooled SD       : %.3f\ndelta / sigma   : %.4f\n\n", paste(sprintf("%.3f", row_means), collapse = "  "), paste(sprintf("%.3f", diff(row_means)), collapse = "  "), max_delta, pooled_sd, max_delta / pooled_sd))

# =============================================================================
# 2. DESIGN CLASS & EXACT UIT POWER EVALUATION
# =============================================================================

lfc_list <- lapply(1:(K - 1), function(k) rep(c(rep(max_delta, k), rep(0, K - k)), each = C))
make_design <- function(w) { w1 <- min(max(w[1], 0.01), 0.49); p <- c(w1, 1 - 2 * w1, w1); p / sum(p) }

uit_power <- function(p, z, sd, N_ref = N_total) {
  p <- p / sum(p); if (any(p <= 0)) return(0)
  m <- K - 1; rs <- C^2 / p; z_mat <- matrix(z, nrow = K, ncol = C, byrow = TRUE)
  
  mu <- sapply(1:m, function(i) sqrt(N_ref) * (sum(z_mat[i + 1, ]) - sum(z_mat[i, ])) / (sd * sqrt(rs[i] + rs[i + 1])))
  Sigma <- diag(m)
  if (m > 1) {
    for (i in 1:(m - 1)) {
      Sigma[i, i + 1] <- Sigma[i + 1, i] <- -rs[i + 1] / sqrt((rs[i] + rs[i + 1]) * (rs[i + 1] + rs[i + 2]))
    }
  }
  
  algo <- if (m == 2) Miwa(steps = 512) else GenzBretz(maxpts = 250000, abseps = 1e-8)
  nr_prob <- as.numeric(pmvnorm(lower = rep(-c_alpha, m), upper = rep(c_alpha, m), mean = mu, sigma = Sigma, algorithm = algo)[1])
  min(max(1 - nr_prob, 0), 1)
}

expected_power <- function(q, p, sd, lfc_list, N_ref = N_total) {
  sum(q * sapply(lfc_list, function(z) uit_power(p, z, sd, N_ref)))
}

# =============================================================================
# 3. REDUCED SYMMETRIC H-ALGORITHM & GRID VERIFICATION
# =============================================================================

optimum_average_design <- function(q, sd, lfc_list, N_ref = N_total) {
  res <- nloptr(x0 = 1 / K, eval_f = function(w) -expected_power(q, make_design(w), sd, lfc_list, N_ref), lb = 0.02, ub = 0.48, opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-8, maxeval = 1000))
  list(p = make_design(res$solution), B = -res$objective, w = res$solution)
}

h_algorithm_k3_symmetric <- function(lfc_list, sd, N_ref = N_total) {
  q <- c(0.5, 0.5); fit <- optimum_average_design(q, sd, lfc_list, N_ref); p <- fit$p
  powers <- sapply(lfc_list, function(z) uit_power(p, z, sd, N_ref))
  B_avg <- sum(q * powers); min_power <- min(powers); C1_gap <- B_avg - min_power
  
  cat("===============================================================================\nREDUCED SYMMETRIC H-ALGORITHM (K = 3)\n===============================================================================\n")
  cat(sprintf("Prior q       : %s\nDesign p      : %s\nLFC powers    : %s\nB_avg         : %.6f\nMin power     : %.6f\nC1 gap        : %.3e\nConverged     : %s\n", paste(sprintf("%.6f", q), collapse = " "), paste(sprintf("%.6f", p), collapse = " "), paste(sprintf("%.6f", powers), collapse = " "), B_avg, min_power, C1_gap, isTRUE(C1_gap <= 1e-8)))
  list(p = p, B_avg = B_avg)
}

verify_optimum_k3 <- function(q, sd, lfc_list, N_ref = N_total, n_grid = 2000) {
  w_grid <- sort(seq(0.02, 0.48, length.out = n_grid))
  B_vals <- sapply(w_grid, function(w) expected_power(q, make_design(w), sd, lfc_list, N_ref))
  i <- which.max(B_vals); list(w = w_grid[i], p = make_design(w_grid[i]), B = B_vals[i])
}

fit <- h_algorithm_k3_symmetric(lfc_list, pooled_sd, N_total)
grid_fit <- verify_optimum_k3(c(0.5, 0.5), pooled_sd, lfc_list)
p_opt <- if (abs(grid_fit$w - fit$p[1]) > 1e-4) { cat("\n[!] COBYLA and grid disagree — using grid value.\n"); grid_fit$p } else { fit$p }

cat(sprintf("\nGrid optimum w : %.6f (B = %.6f)\nCOBYLA w       : %.6f (B = %.6f)\nAbsolute diff  : %.3e\n", grid_fit$w, grid_fit$B, fit$p[1], fit$B_avg, abs(grid_fit$w - fit$p[1])))

# =============================================================================
# 4. ROBUST INTEGER ALLOCATION
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
N_bal <- allocate_integers(rep(1 / K, K), N_total)

cat("\n===============================================================================\nEXACT INTEGER DESIGNS AT N =", N_total, "\n===============================================================================\n")
cat("\nMax-min (cell-wise rounding):\n"); print(N_opt)
cat(sprintf("  Row totals: %s  (sum = %d)\n\nBalanced:\n", paste(rowSums(N_opt), collapse = "  "), sum(N_opt)))
print(N_bal)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = "  "), sum(N_bal)))

# =============================================================================
# 5. MONTE CARLO POWER SIMULATION & SAMPLE-SIZE SEARCH
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L; row_vars <- rowSums(sd^2 / N_mat) / (C^2)
  diff_se <- sqrt(row_vars[1:(K - 1)] + row_vars[2:K])
  
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    for (i in 1:K) for (j in 1:C) Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
    if (any(abs(diff(rowMeans(Ybar)) / diff_se) > c_alpha)) ctr <- ctr + 1L
  }
  return(ctr / n_sim)
}


power_opt <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim)
power_bal <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim)
se_opt <- sqrt(power_opt * (1 - power_opt) / n_sim); se_bal <- sqrt(power_bal * (1 - power_bal) / n_sim)
gain <- 100 * (power_opt - power_bal) / power_bal

cat("\n===============================================================================\nEMPIRICAL POWER AT N =", N_total, "\n===============================================================================\n")
cat(sprintf("Bonferroni CV   : %.4f\nSimulation size : %d\nMax-min power   : %.4f (SE %.4f)\nBalanced power  : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", c_alpha, n_sim, power_opt, se_opt, power_bal, se_bal, gain))

find_n90 <- function(p_row, M, sd, c_alpha, target = 0.90, N_min = 20, N_max = 3000) {
  N_val <- N_min; pow <- 0
  cat("  [Coarse step 10] ... ")
  while (pow < target && N_val < N_max) { N_val <- N_val + 10; pow <- simulate_power(allocate_integers(p_row, N_val), M, sd, c_alpha, 5000) }
  cat(sprintf("reached ~%d.  [Fine step 1] ... ", N_val))
  
  N_val <- max(N_min, N_val - 15); pow <- 0
  while (pow < target && N_val < N_max) { N_val <- N_val + 1; pow <- simulate_power(allocate_integers(p_row, N_val), M, sd, c_alpha, 20000) }
  cat(sprintf("N = %d (power %.2f%%)\n", N_val, 100 * pow))
  list(N = N_val, power = pow)
}

cat("\n===============================================================================\nSAMPLE SIZE FOR 90% POWER\n===============================================================================\n")
cat("> Max-min design:\n"); n90_opt <- find_n90(p_opt, M_sampled, pooled_sd, c_alpha)
cat("> Balanced design:\n"); n90_bal <- find_n90(rep(1 / K, K), M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 6. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (UIT Successive Comparison) ══════════════════════\n")
cat(sprintf(" Max adjacent gap (delta) : %.3f\n Pooled SD (sigma)        : %.3f\n delta / sigma            : %.4f\n", max_delta, pooled_sd, max_delta / pooled_sd))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain            : %+.2f%%\n", N_total, power_opt, N_total, power_bal, gain))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction)\n", n90_opt$N, n90_opt$power, n90_bal$N, n90_bal$power, n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(" Max-min integer design:\n"); print(N_opt); cat(sprintf("  Row totals: %s  (sum = %d)\n\n Balanced integer design:\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt))); print(N_bal); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))
cat("═══════════════════════════════════════════════════════════════════════════════════════\n")