# =============================================================================
# UIT SUCCESSIVE-COMPARISON DESIGN (K = 3)
#
# Computational structure:
# 1. Data load & feature engineering (3 Exact Tertiles of SO2 × Season)
# 2. Fixed population cell means and pooled SD
# 3. Robust integer allocation (Hamilton's method with tie-break)
# 4. UIT power evaluation (Deterministic Miwa integration)
# 5. Reduced symmetric H-algorithm (Max-min design for successive gaps)
# 6. Monte Carlo empirical power & sample size search (Max-Min vs Balanced)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(mvtnorm)
  library(nloptr)
})

set.seed(27092026)

# =============================================================================
# 1. DATA LOAD & FEATURE ENGINEERING (3 EXACT TERTILES)
# =============================================================================

df <- read.csv("cpcb_dly_aq_uttar_pradesh-2011.csv", check.names = FALSE, na.strings = c("NA", ""))
colnames(df) <- gsub("[ /]", "_", colnames(df))

df$Sampling_Date <- as.Date(df$Sampling_Date, format = "%d/%m/%Y")
df$Month <- as.numeric(format(df$Sampling_Date, "%m"))
df$SO2 <- as.numeric(df$SO2)
df$RSPM_PM10 <- as.numeric(df$RSPM_PM10)

# --- Generate exactly 3 groups using Tertiles (33.33% and 66.67%) ---
quants <- quantile(df$SO2, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

df_clean <- df %>%
  filter(!is.na(RSPM_PM10), !is.na(SO2), !is.na(Month)) %>%
  mutate(
    SO2_Tier = case_when(
      SO2 <= quants[2]                  ~ "Tier 1: Low SO2",
      SO2 >  quants[2] & SO2 <= quants[3] ~ "Tier 2: Mid SO2",
      SO2 >  quants[3]                  ~ "Tier 3: High SO2"
    ),
    Season = case_when(
      Month %in% c(7, 8, 9, 10)  ~ "Season 1: Monsoon",
      Month %in% c(3, 4, 5, 6)   ~ "Season 2: Summer",
      Month %in% c(11, 12, 1, 2) ~ "Season 3: Winter"
    )
  ) %>%
  filter(!is.na(SO2_Tier), !is.na(Season))

factor_a_levels <- c("Tier 1: Low SO2", "Tier 2: Mid SO2", "Tier 3: High SO2")
factor_b_levels <- c("Season 1: Monsoon", "Season 2: Summer", "Season 3: Winter")

# --- Model Parameters ---
K <- length(factor_a_levels)
C <- length(factor_b_levels)
m <- K - 1
N_total <- 198 ##200
alpha <- 0.05
n_sim <- 100000
c_alpha <- qnorm(1 - alpha / (2 * m)) # Bonferroni correction for UIT

cat("===============================================================================\n")
cat("SETUP SUMMARY & CELL COUNTS (3 SO2 TERTILES x SEASON)\n")
cat("===============================================================================\n")
cat(sprintf("K = %d, C = %d, m = %d, N_total = %d, Bonferroni CV = %.4f\n\n", K, C, m, N_total, c_alpha))
print(table(df_clean$SO2_Tier, df_clean$Season))

# =============================================================================
# 2. FIXED POPULATION CELL MEANS & STATISTICS
# =============================================================================

cell_stats <- df_clean %>%
  group_by(SO2_Tier, Season) %>%
  summarise(
    n_obs = n(),
    mean  = mean(RSPM_PM10, na.rm = TRUE),
    sd    = sd(RSPM_PM10, na.rm = TRUE),
    .groups = "drop"
  )

M_sampled <- cell_stats %>%
  select(SO2_Tier, Season, mean) %>%
  pivot_wider(names_from = Season, values_from = mean) %>%
  as.data.frame()

rownames(M_sampled) <- M_sampled$SO2_Tier
M_sampled <- as.matrix(M_sampled[, factor_b_levels, drop = FALSE])
M_sampled <- M_sampled[factor_a_levels, , drop = FALSE]

ss_within <- sum((cell_stats$n_obs - 1) * cell_stats$sd^2, na.rm = TRUE)
df_within <- sum(cell_stats$n_obs - 1, na.rm = TRUE)
pooled_sd <- sqrt(ss_within / df_within)

row_means <- rowMeans(M_sampled)
adj_diffs <- diff(row_means)
max_delta <- max(abs(adj_diffs))

cat("\n===============================================================================\n")
cat("FULL-DATA CELL MEANS (PM10, µg/m³)\n")
cat("===============================================================================\n\n")
print(round(M_sampled, 3))
cat(sprintf("\nσ̂ (pooled) = %.3f | δ = %.3f | δ/σ = %.4f\n\n", pooled_sd, max_delta, max_delta / pooled_sd))

# =============================================================================
# 3. ROBUST INTEGER ALLOCATION
# =============================================================================

allocate_integers <- function(x, total_N, seed = 26092026) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  if (is.null(dim(x))) x <- outer(x, rep(1 / C, C)) * total_N
  
  N_mat <- floor(x)
  rem <- total_N - sum(N_mat)
  if (rem > 0) {
    set.seed(seed)
    frac <- as.vector(x - floor(x))
    uniq_fracs <- sort(unique(round(frac, 12)), decreasing = TRUE)
    idx_selected <- integer(0)
    remaining <- rem
    for (f in uniq_fracs) {
      if (remaining <= 0) break
      group <- which(abs(frac - f) < 1e-12)
      take <- min(length(group), remaining)
      chosen <- if (take == length(group)) group else sample(group, take)
      idx_selected <- c(idx_selected, chosen)
      remaining <- remaining - take
    }
    N_mat[idx_selected] <- N_mat[idx_selected] + 1L
  }
  if (!is.null(dim(N_mat))) dimnames(N_mat) <- dimnames(M_sampled)
  return(N_mat)
}

# =============================================================================
# 4. DETERMINISTIC UIT POWER & MAX-MIN DESIGN (H-Algorithm)
# =============================================================================

uit_power <- function(p, z, sd, N_ref = N_total) {
  p <- p / sum(p); if (any(p <= 0)) stop("All row proportions must be positive.")
  z_mat <- matrix(z, nrow = K, ncol = C, byrow = TRUE); rs <- C^2 / p
  
  mu <- numeric(m)
  for (i in seq_len(m)) {
    se <- sqrt(rs[i] + rs[i + 1])
    mu[i] <- sqrt(N_ref) * (sum(z_mat[i + 1, ]) - sum(z_mat[i, ])) / (sd * se)
  }
  
  Sigma <- diag(m)
  if (m > 1) {
    for (i in seq_len(m - 1)) {
      se_i <- sqrt(rs[i] + rs[i + 1]); se_next <- sqrt(rs[i + 1] + rs[i + 2])
      Sigma[i, i + 1] <- -rs[i + 1] / (se_i * se_next); Sigma[i + 1, i] <- Sigma[i, i + 1]
    }
  }
  
  nr_prob <- as.numeric(pmvnorm(lower = rep(-c_alpha, m), upper = rep(c_alpha, m), mean = mu, sigma = Sigma, algorithm = Miwa(steps = 512)))
  min(max(1 - nr_prob, 0), 1)
}

make_design <- function(w) { w1 <- w[1]; p <- c(w1, 1 - 2 * w1, w1); p / sum(p) }
expected_power <- function(q, p, sd, lfc_list, N_ref = N_total) { sum(q * sapply(lfc_list, function(z) uit_power(p, z, sd, N_ref))) }

optimum_average_design <- function(q, sd, lfc_list, N_ref = N_total) {
  res <- nloptr(x0 = 1 / K, eval_f = function(w) -expected_power(q, make_design(w), sd, lfc_list, N_ref), lb = 0.02, ub = 0.48, opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-8, maxeval = 1000))
  list(p = make_design(res$solution), B = -res$objective, w = res$solution)
}

h_algorithm_symmetric_k3 <- function(lfc_list, sd, N_ref = N_total, tol = 1e-8, verbose = TRUE) {
  q <- c(0.5, 0.5)
  fit <- optimum_average_design(q, sd, lfc_list, N_ref)
  p_opt <- fit$p; powers <- vapply(lfc_list, function(z) uit_power(p_opt, z, sd, N_ref), numeric(1))
  B_avg <- sum(q * powers); min_power <- min(powers); gap <- B_avg - min_power
  if (verbose) {
    cat(sprintf("  p (design) = %s\n  LFC powers = %s\n  B_avg = %.6f | min = %.6f | gap = %.3e\n", paste(sprintf("%.6f", p_opt), collapse = " "), paste(sprintf("%.6f", powers), collapse = " "), B_avg, min_power, gap))
  }
  list(p = p_opt, q = q, B_avg = B_avg, min_power = min_power, pow_each = powers, gap = gap, converged = (gap <= tol))
}

lfc_list <- lapply(1:m, function(k) { rep(c(rep(max_delta, k), rep(0, K - k)), each = C) })

cat("\n===============================================================================\nREDUCED SYMMETRIC H-ALGORITHM (K = 3) — FULL POPULATION\n===============================================================================\n")
fit <- h_algorithm_symmetric_k3(lfc_list, pooled_sd, N_ref = N_total, tol = 1e-8, verbose = TRUE)
p_opt <- fit$p; p_bal <- rep(1 / K, K)

# =============================================================================
# 5. INTEGER DESIGNS & EMPIRICAL POWER SIMULATION
# =============================================================================

N_opt <- allocate_integers(p_opt, N_total)
N_bal <- allocate_integers(p_bal, N_total)

cat("\n===============================================================================\nEXACT INTEGER DESIGNS AT N =", N_total, "\n===============================================================================\n")
cat("\nMax-min integer design:\n"); print(N_opt); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt)))
cat("\nBalanced integer design:\n"); print(N_bal); cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    for (i in 1:K) for (j in 1:C) Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
    row_means <- rowMeans(Ybar); row_vars <- rowSums(sd^2 / N_mat) / (C^2)
    diff_se <- sqrt(row_vars[1:m] + row_vars[2:K])
    if (any(abs(diff(row_means) / diff_se) > c_alpha)) ctr <- ctr + 1L
  }
  return(ctr / n_sim)
}

power_bal <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim)
power_opt <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim)
se_bal <- sqrt(power_bal * (1 - power_bal) / n_sim); se_opt <- sqrt(power_opt * (1 - power_opt) / n_sim)
gain   <- 100 * (power_opt - power_bal) / power_bal

cat("\n===============================================================================\nEMPIRICAL POWER AT N =", N_total, "\n===============================================================================\n")
cat(sprintf("Simulation size : %d\nMax-min power   : %.4f (SE %.4f)\nBalanced power  : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", n_sim, power_opt, se_opt, power_bal, se_bal, gain))

# =============================================================================
# 6. SAMPLE SIZE SEARCH FOR 90% POWER
# =============================================================================

find_n90 <- function(p_row, M, sd, c_alpha, target = 0.90, N_min = 30, N_max = 3000) {
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
cat("> Balanced design:\n"); n90_bal <- find_n90(p_bal, M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 7. FINAL SUMMARY
# =============================================================================

cat("\n════════════════════ FINAL SUMMARY (UIT Successive Comparison) ════════════════════\n")
cat(sprintf(" Max gap (delta)          : %.3f\n Pooled SD (sigma)        : %.3f\n delta / sigma            : %.4f\n Row proportions p        : %s\n B_avg / min LFC power    : %.6f / %.6f\n", max_delta, pooled_sd, max_delta / pooled_sd, paste(sprintf("%.6f", p_opt), collapse = " "), fit$B_avg, fit$min_power))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain (vs Bal)   : %+.2f%%\n", N_total, power_opt, N_total, power_bal, gain))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(" 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction vs Bal)\n", n90_opt$N, n90_opt$power, n90_bal$N, n90_bal$power, n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N))
cat("══════════════════════════════════════════════════════════════════════════════════════\n")
