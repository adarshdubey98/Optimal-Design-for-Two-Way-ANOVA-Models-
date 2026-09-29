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

suppressPackageStartupMessages({
  library(mvtnorm)
  library(dplyr)
  library(tidyr)
})

set.seed(27092026)

# =============================================================================
# 1. FIXED INPUT & SUMMARY STATISTICS
# =============================================================================

# --- 1. Data Load ---
df <- read.csv("cpcb_dly_aq_uttar_pradesh-2011.csv", check.names = FALSE)
df$`RSPM/PM10` <- as.numeric(as.character(df$`RSPM/PM10`))

# NOTE: change this to whatever your CSV calls the city column
# (commonly "City/Town/Village/Area" in the CPCB open-data extracts).
city_col <- "City/Town/Village/Area"

df_clean <- df %>%
  filter(!is.na(`RSPM/PM10`), !is.na(`Type of Location`), !is.na(.data[[city_col]]))

# --- 2. Factor Levels Setup ---

factor_a_levels <- c("Kanpur", "Lucknow", "Mathura", "Khurja", "Rai Bareilly", "Unnao")
factor_b_levels <- c("Industrial Area", "Residential, Rural and other Areas")

df_clean <- df_clean %>%
  filter(.data[[city_col]] %in% factor_a_levels) %>%
  mutate(
    Factor_A = factor(.data[[city_col]], levels = factor_a_levels),
    Factor_B = if_else(grepl("Industrial", `Type of Location`),
                       "Industrial Area", "Residential, Rural and other Areas"),
    Factor_B = factor(Factor_B, levels = factor_b_levels)
  )

# --- Cell-wise full-data means ---
cell_means_full <- df_clean %>%
  group_by(Factor_A, Factor_B) %>%
  summarise(
    n_obs = n(),
    mean  = mean(`RSPM/PM10`, na.rm = TRUE),
    sd    = sd(`RSPM/PM10`, na.rm = TRUE),
    .groups = "drop"
  )

# --- Wide matrix (base R) ---
M_full <- matrix(
  NA,
  nrow = length(factor_a_levels),
  ncol = length(factor_b_levels),
  dimnames = list(factor_a_levels, factor_b_levels)
)

N_full <- matrix(
  NA_integer_,
  nrow = length(factor_a_levels),
  ncol = length(factor_b_levels),
  dimnames = list(factor_a_levels, factor_b_levels)
)

for (i in seq_along(factor_a_levels)) {
  for (j in seq_along(factor_b_levels)) {
    vals <- df_clean$`RSPM/PM10`[
      df_clean$Factor_A == factor_a_levels[i] &
      df_clean$Factor_B == factor_b_levels[j]
    ]
    M_full[i, j] <- if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA
    N_full[i, j] <- length(vals)
  }
}

# --- Print Summaries ---
cat("=========================================================\n")
cat("FULL-DATA CELL MEANS (RSPM/PM10, µg/m³)\n")
cat("=========================================================\n\n")
print(round(M_full, 2))

cat("\n--- Cell sample sizes ---\n")
print(N_full)

# --- Pooled within-cell σ ---
ss_within <- sum((cell_means_full$n_obs - 1) * cell_means_full$sd^2, na.rm = TRUE)
df_within <- sum(cell_means_full$n_obs - 1, na.rm = TRUE)
sigma_full <- sqrt(ss_within / df_within)
cat(sprintf("\nσ̂ (pooled) = %.3f\n", sigma_full))

# --- Row means, gaps, δ/σ ---
row_means_full <- rowMeans(M_full, na.rm = TRUE)
ct_gaps_full   <- row_means_full[-1] - row_means_full[1]
delta_full     <- max(abs(ct_gaps_full))

cat("\n--- Row means ---\n")
print(round(row_means_full, 3))
cat("\n--- Control vs treatment gaps ---\n")
print(round(ct_gaps_full, 3))
cat(sprintf("\nδ = %.3f\nσ = %.3f\nδ/σ = %.4f\n", delta_full, sigma_full, delta_full / sigma_full))

# --- Model Parameters ---
M_sampled <- M_full 
pooled_sd <- sigma_full 
K <- nrow(M_sampled)
C <- ncol(M_sampled)
m <- K - 1
n_sim <- 100000
N_total <- 78
alpha <- 0.05
c_alpha <- qnorm(1 - alpha / (2 * m)) # Bonferroni correction

row_means <- rowMeans(M_sampled)
ct_gaps <- row_means[-1] - row_means[1]
max_delta <- max(abs(ct_gaps))

cat("\n===============================================================================\n")
cat("FIXED INPUT SUMMARY (UIT CONTROL VS TREATMENT)\n")
cat("===============================================================================\n")
cat("\nCell mean matrix (PM10, µg/m³):\n")
print(round(M_sampled, 2))
cat(sprintf(
  "\nRow means           : %s\nControl-vs-Trt gaps : %s\nMax CT gap (delta)  : %.3f\nPooled SD (sigma)   : %.3f\ndelta / sigma       : %.4f\n", 
  paste(sprintf("%.2f", row_means), collapse = "  "), 
  paste(sprintf("%.2f", ct_gaps), collapse = "  "), 
  max_delta, 
  pooled_sd, 
  max_delta / pooled_sd
))

# =============================================================================
# 2. DESIGN CLASS & WORST-CASE UIT POWER EVALUATION
# =============================================================================

design_CT <- function(eta) {
  c(C * eta, rep((1 - C * eta) / (K - 1), K - 1))
}

kappa_rho_CT <- function(eta) {
  kappa <- sqrt(N_total * C * eta * (1 - C * eta)) / (pooled_sd * sqrt(1 + C * eta * (K - 2)))
  rho <- (1 - C * eta) / (1 + C * eta * (K - 2))
  list(kappa = kappa, rho = rho)
}

power_UIT_CT <- function(Delta, eta) {
  kr <- kappa_rho_CT(eta)
  Sigma <- matrix(kr$rho, m, m)
  diag(Sigma) <- 1
  nu <- -kr$kappa * Delta
  
  1 - as.numeric(pmvnorm(lower = rep(-c_alpha, m), upper = rep(c_alpha, m), mean = nu, sigma = Sigma)[1])
}

Pi_u_eta <- function(u, eta, delta) {
  power_UIT_CT(c(delta, rep(u, m - 1)), eta)
}

worst_case_eta <- function(eta, delta, n_coarse = 60) {
  u_grid <- seq(0, delta, length.out = n_coarse)
  vals <- sapply(u_grid, function(u) Pi_u_eta(u, eta, delta))
  best <- which.min(vals)
  
  lo <- u_grid[max(1, best - 1)]
  hi <- u_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(u) Pi_u_eta(u, eta, delta), interval = c(lo, hi))
  
  min(vals[best], opt$objective)
}

optimize_eta <- function(delta, n_coarse = 80) {
  eta_grid <- seq(1e-3, 1 / C - 1e-3, length.out = n_coarse)
  vals <- sapply(eta_grid, function(e) worst_case_eta(e, delta))
  best <- which.max(vals)
  
  lo <- eta_grid[max(1, best - 1)]
  hi <- eta_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(e) worst_case_eta(e, delta), interval = c(lo, hi), maximum = TRUE)
  
  if (opt$objective >= vals[best]) {
    list(eta = opt$maximum, power = opt$objective)
  } else {
    list(eta = eta_grid[best], power = vals[best])
  }
}

cat("\n===============================================================================\n")
cat("MAX-MIN DESIGN OPTIMIZATION\n")
cat("===============================================================================\n")

res <- optimize_eta(max_delta)
eta_opt <- res$eta
p_opt <- design_CT(eta_opt)

eta_dunnett <- 1 / (C * (1 + sqrt(K - 1)))
p_dunnett <- design_CT(eta_dunnett)

p_bal <- rep(1 / K, K)

cat(sprintf(
  "Optimal eta             : %.6f (range: 0 to %.4f)\nWorst-case LFC power    : %.6f\nMax-min row proportions : %s\nDunnett row proportions : %s\n", 
  eta_opt, 1/C, res$power, 
  paste(sprintf("%.6f", p_opt), collapse = "  "), 
  paste(sprintf("%.6f", p_dunnett), collapse = "  ")
))

# =============================================================================
# 3. ROBUST INTEGER ALLOCATION
# =============================================================================

allocate_integers <- function(x, total_N, seed = 26092026) {
  # --- RNG state save/restore (to avoid disturbing the global stream) ---
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv),
      add = TRUE
    )
  }
  
  # --- If x is a vector, expand it into a cell matrix ---
  if (is.null(dim(x))) {
    x <- outer(x, rep(1 / C, C)) * total_N
  }
  
  # --- Base allocation ---
  N_mat <- floor(x)
  rem <- total_N - sum(N_mat)
  
  if (rem > 0) {
    set.seed(seed)   # temporary local stream
    frac <- as.vector(x - floor(x))
    
    # Step 1: unique fractional values, descending
    uniq_fracs <- sort(unique(round(frac, 12)), decreasing = TRUE)
    
    # Step 2: group-wise increments
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
  
  # --- Apply dimnames only if it's a matrix ---
  if (!is.null(dim(N_mat))) {
    dimnames(N_mat) <- dimnames(M_sampled)
  }
  
  return(N_mat)
}

N_opt <- allocate_integers(p_opt, N_total)
N_dun <- allocate_integers(p_dunnett, N_total)
N_bal <- allocate_integers(p_bal, N_total)

cat("\n===============================================================================\n")
cat("EXACT INTEGER DESIGNS AT N =", N_total, "\n")
cat("===============================================================================\n")

cat("\nMax-min integer design:\n")
print(N_opt)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt)))

cat("\nDunnett integer design:\n")
print(N_dun)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_dun), collapse = " "), sum(N_dun)))

cat("\nBalanced integer design:\n")
print(N_bal)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))

# =============================================================================
# 4. MONTE CARLO POWER SIMULATION & SAMPLE-SIZE SEARCH
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L
  row_vars <- rowSums(sd^2 / N_mat) / (C^2)
  
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    
    for (i in 1:K) {
      for (j in 1:C) {
        Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
      }
    }
    
    row_means <- rowMeans(Ybar)
    
    if (any(abs((row_means[1] - row_means[2:K]) / sqrt(row_vars[1] + row_vars[2:K])) > c_alpha)) {
      ctr <- ctr + 1L
    }
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

cat("\n===============================================================================\n")
cat("EMPIRICAL POWER AT N =", N_total, "\n")
cat("===============================================================================\n")
cat(sprintf(
  "Bonferroni CV   : %.4f\nSimulation size : %d\nMax-min power   : %.4f (SE %.4f)\nDunnett power   : %.4f (SE %.4f)\nBalanced power  : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", 
  c_alpha, n_sim, power_opt, se_opt, power_dun, se_dun, power_bal, se_bal, gain
))

find_n90 <- function(p_row, M, sd, c_alpha, target = 0.90, N_min = 20, N_max = 3000) {
  N_val <- N_min
  pow <- 0
  cat("  [Coarse step 10] ... ")
  
  while (pow < target && N_val < N_max) { 
    N_val <- N_val + 10
    pow <- simulate_power(allocate_integers(p_row, N_val), M, sd, c_alpha, 5000) 
  }
  cat(sprintf("reached ~%d.  [Fine step 1] ... ", N_val))
  
  N_val <- max(N_min, N_val - 15)
  pow <- 0
  
  while (pow < target && N_val < N_max) { 
    N_val <- N_val + 1
    pow <- simulate_power(allocate_integers(p_row, N_val), M, sd, c_alpha, 20000) 
  }
  cat(sprintf("N = %d (power %.2f%%)\n", N_val, 100 * pow))
  
  list(N = N_val, power = pow)
}

cat("\n===============================================================================\n")
cat("SAMPLE SIZE FOR 90% POWER\n")
cat("===============================================================================\n")
cat("> Max-min design:\n")
n90_opt <- find_n90(p_opt, M_sampled, pooled_sd, c_alpha)
cat("> Dunnett design:\n")
n90_dun <- find_n90(p_dunnett, M_sampled, pooled_sd, c_alpha)
cat("> Balanced design:\n")
n90_bal <- find_n90(p_bal, M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 5. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (UIT Control vs Treatment) ══════════════════════\n")
cat(sprintf(
  " Max CT gap (delta)       : %.3f\n Pooled SD (sigma)        : %.3f\n delta / sigma            : %.4f\n", 
  max_delta, pooled_sd, max_delta / pooled_sd
))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Dunnett     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain (vs Bal)   : %+.2f%%\n", 
  N_total, power_opt, N_total, power_dun, N_total, power_bal, gain
))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Dunnett           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction vs Bal)\n", 
  n90_opt$N, n90_opt$power, n90_dun$N, n90_dun$power, n90_bal$N, n90_bal$power, 
  n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N
))
cat("══════════════════════════════════════════════════════════════════════════════════════\n")
