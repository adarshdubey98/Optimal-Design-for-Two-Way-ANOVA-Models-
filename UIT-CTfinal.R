# =============================================================================
# UIT CONTROL-VERSUS-TREATMENT DESIGN (KANPUR VS 5 CITIES)
#
# Computational structure:
# 1. Fixed cell means and pooled SD (Table 1)
# 2. Complete two-sided UIT power evaluation (Bonferroni corrected)
# 3. Max-min design optimization (eta parameter, over LFCs)
# 4. Integer allocation (Hamilton's method)
# 5. Monte Carlo empirical power & sample size search (Max-Min vs Dunnett vs Bal)
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages(library(mvtnorm))

# =============================================================================
# 1. FIXED INPUT & SUMMARY STATISTICS
# =============================================================================

M_sampled <- matrix(c(
  200.80, 136.20,   # Kanpur (control)
  202.00, 184.20,   # Lucknow
  206.00, 174.20,   # Mathura
  177.80, 158.80,   # Khurja
  153.94, 137.44,   # Rai Bareilly
  135.80, 137.40    # Unnao
), nrow = 6, ncol = 2, byrow = TRUE)

rownames(M_sampled) <- c("Kanpur (control)", "Lucknow", "Mathura", "Khurja", "Rai Bareilly", "Unnao")
colnames(M_sampled) <- c("Industrial", "Residential/Rural")

pooled_sd <- 39.40
K <- nrow(M_sampled); C <- ncol(M_sampled); m <- K - 1;n_sim <- 100000
N_total <- 60; alpha <- 0.05; c_alpha <- qnorm(1 - alpha / (2 * m)) # Bonferroni correction

row_means <- rowMeans(M_sampled)
ct_gaps <- row_means[-1] - row_means[1]
max_delta <- max(abs(ct_gaps))

cat("===============================================================================\nFIXED INPUT SUMMARY (UIT CONTROL VS TREATMENT)\n===============================================================================\n")
cat("\nCell mean matrix (PM10, µg/m³):\n"); print(round(M_sampled, 2))
cat(sprintf("\nRow means           : %s\nControl-vs-Trt gaps : %s\nMax CT gap (delta)  : %.3f\nPooled SD (sigma)   : %.3f\ndelta / sigma       : %.4f\n", paste(sprintf("%.2f", row_means), collapse = "  "), paste(sprintf("%.2f", ct_gaps), collapse = "  "), max_delta, pooled_sd, max_delta / pooled_sd))

# =============================================================================
# 2. DESIGN CLASS & WORST-CASE UIT POWER EVALUATION
# =============================================================================

design_CT <- function(eta) c(C * eta, rep((1 - C * eta) / (K - 1), K - 1))

kappa_rho_CT <- function(eta) {
  kappa <- sqrt(N_total * C * eta * (1 - C * eta)) / (pooled_sd * sqrt(1 + C * eta * (K - 2)))
  rho <- (1 - C * eta) / (1 + C * eta * (K - 2))
  list(kappa = kappa, rho = rho)
}

power_UIT_CT <- function(Delta, eta) {
  kr <- kappa_rho_CT(eta); Sigma <- matrix(kr$rho, m, m); diag(Sigma) <- 1
  nu <- -kr$kappa * Delta
  1 - as.numeric(pmvnorm(lower = rep(-c_alpha, m), upper = rep(c_alpha, m), mean = nu, sigma = Sigma)[1])
}

Pi_u_eta <- function(u, eta, delta) power_UIT_CT(c(delta, rep(u, m - 1)), eta)

worst_case_eta <- function(eta, delta, n_coarse = 60) {
  u_grid <- seq(0, delta, length.out = n_coarse)
  vals <- sapply(u_grid, function(u) Pi_u_eta(u, eta, delta)); best <- which.min(vals)
  lo <- u_grid[max(1, best - 1)]; hi <- u_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(u) Pi_u_eta(u, eta, delta), interval = c(lo, hi))
  min(vals[best], opt$objective)
}

optimize_eta <- function(delta, n_coarse = 80) {
  eta_grid <- seq(1e-3, 1 / C - 1e-3, length.out = n_coarse)
  vals <- sapply(eta_grid, function(e) worst_case_eta(e, delta)); best <- which.max(vals)
  lo <- eta_grid[max(1, best - 1)]; hi <- eta_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(e) worst_case_eta(e, delta), interval = c(lo, hi), maximum = TRUE)
  if (opt$objective >= vals[best]) list(eta = opt$maximum, power = opt$objective) else list(eta = eta_grid[best], power = vals[best])
}

cat("\n===============================================================================\nMAX-MIN DESIGN OPTIMIZATION\n===============================================================================\n")
res <- optimize_eta(max_delta)
eta_opt <- res$eta; p_opt <- design_CT(eta_opt)
eta_dunnett <- 1 / (C * (1 + sqrt(K - 1))); p_dunnett <- design_CT(eta_dunnett)
p_bal <- rep(1 / K, K)

cat(sprintf("Optimal eta             : %.6f (range: 0 to %.4f)\nWorst-case LFC power    : %.6f\nMax-min row proportions : %s\nDunnett row proportions : %s\n", eta_opt, 1/C, res$power, paste(sprintf("%.6f", p_opt), collapse = "  "), paste(sprintf("%.6f", p_dunnett), collapse = "  ")))

# =============================================================================
# 3. ROBUST INTEGER ALLOCATION
# =============================================================================

allocate_integers <- function(p, total_N) {
  # 1. Decimal targets calculate karein
  target <- outer(p, rep(1 / C, C)) * total_N
  
  # 2. Base allocation 
  N_mat <- floor(target)
  
  # 3. remainder observations
  rem <- total_N - sum(N_mat)
  
  # 4.  remainder observations allocation largest fractions cell
  if (rem > 0) {
    frac <- as.vector(target - floor(target))
    idx <- order(frac, decreasing = TRUE, method = "radix")[seq_len(rem)]
    N_mat[idx] <- N_mat[idx] + 1L
  }
  
  dimnames(N_mat) <- dimnames(M_sampled)
  return(N_mat)
}

N_opt <- allocate_integers(p_opt, N_total)
N_dun <- allocate_integers(p_dunnett, N_total)
N_bal <- allocate_integers(p_bal, N_total)

cat("\n===============================================================================\nEXACT INTEGER DESIGNS AT N =", N_total, "\n===============================================================================\n")
cat("\nMax-min integer design:\n"); print(N_opt); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt)))
cat("\nDunnett integer design:\n"); print(N_dun); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_dun), collapse = " "), sum(N_dun)))
cat("\nBalanced integer design:\n"); print(N_bal); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))

# =============================================================================
# 4. MONTE CARLO POWER SIMULATION & SAMPLE-SIZE SEARCH
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L; row_vars <- rowSums(sd^2 / N_mat) / (C^2)
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    for (i in 1:K) for (j in 1:C) Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
    row_means <- rowMeans(Ybar)
    if (any(abs((row_means[1] - row_means[2:K]) / sqrt(row_vars[1] + row_vars[2:K])) > c_alpha)) ctr <- ctr + 1L
  }
  return(ctr / n_sim)
}


power_opt <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim)
power_dun <- simulate_power(N_dun, M_sampled, pooled_sd, c_alpha, n_sim)
power_bal <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim)

se_opt <- sqrt(power_opt * (1 - power_opt) / n_sim)
se_dun <- sqrt(power_dun * (1 - power_dun) / n_sim)
se_bal <- sqrt(power_bal * (1 - power_bal) / n_sim)
gain <- 100 * (power_opt - power_bal) / power_bal

cat("\n===============================================================================\nEMPIRICAL POWER AT N =", N_total, "\n===============================================================================\n")
cat(sprintf("Bonferroni CV   : %.4f\nSimulation size : %d\nMax-min power   : %.4f (SE %.4f)\nDunnett power   : %.4f (SE %.4f)\nBalanced power  : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", c_alpha, n_sim, power_opt, se_opt, power_dun, se_dun, power_bal, se_bal, gain))

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
cat("> Dunnett design:\n"); n90_dun <- find_n90(p_dunnett, M_sampled, pooled_sd, c_alpha)
cat("> Balanced design:\n"); n90_bal <- find_n90(p_bal, M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 5. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (UIT Control vs Treatment) ══════════════════════\n")
cat(sprintf(" Max CT gap (delta)       : %.3f\n Pooled SD (sigma)        : %.3f\n delta / sigma            : %.4f\n", max_delta, pooled_sd, max_delta / pooled_sd))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Dunnett     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain (vs Bal)   : %+.2f%%\n", N_total, power_opt, N_total, power_dun, N_total, power_bal, gain))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Dunnett           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction vs Bal)\n", n90_opt$N, n90_opt$power, n90_dun$N, n90_dun$power, n90_bal$N, n90_bal$power, n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N))
cat("══════════════════════════════════════════════════════════════════════════════════════\n")
