# =============================================================================
# IUT SUCCESSIVE COMPARISON DESIGN (NO2 TIERS x LOCATION TYPE)
#
# Computational structure: 
# 1. Data load & feature engineering (3 Exact Quantiles/Tertiles of NO2)
# 2. Fixed cell means and pooled SD (Table 4)
# 3. Exact two-sided IUT power evaluation via multivariate normal
# 4. Max-min design optimization (1D line search)
# 5. Robust integer allocation (Hamilton's method)
# 6. Monte Carlo empirical power & sample size search
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
colnames(df) <- gsub(" ", "_", colnames(df))
colnames(df) <- gsub("/", "_", colnames(df))

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
K <- length(factor_a_levels)   # 3
C <- length(factor_b_levels)   # 3
m <- K - 1                     # 2
N_total <- 65
alpha   <- 0.05
c_alpha <- qnorm(1 - alpha / 2) # IUT critical value
n_sim   <- 100000

cat("===============================================================================\n")
cat("SETUP & CELL COUNTS (3 NO2 TERTILES)\n")
cat("===============================================================================\n")
print(table(df_clean$NO2_Tier, df_clean$Location))

# =============================================================================
# 2. POPULATION PLANNING VALUES & STATISTICS
# =============================================================================

cell_stats <- df_clean %>%
  group_by(NO2_Tier, Location) %>%
  summarise(
    n_obs = n(),
    mean  = mean(RSPM_PM10, na.rm = TRUE),
    sd    = sd(RSPM_PM10, na.rm = TRUE),
    .groups = "drop"
  )

# --- Wide mean matrix ---
M_sampled <- cell_stats %>%
  select(NO2_Tier, Location, mean) %>%
  pivot_wider(names_from = Location, values_from = mean) %>%
  as.data.frame()

rownames(M_sampled) <- M_sampled$NO2_Tier
M_sampled <- as.matrix(M_sampled[, factor_b_levels, drop = FALSE])
M_sampled <- M_sampled[factor_a_levels, , drop = FALSE]

# --- Pooled within-cell SD ---
ss_within <- sum((cell_stats$n_obs - 1) * cell_stats$sd^2, na.rm = TRUE)
df_within <- sum(cell_stats$n_obs - 1, na.rm = TRUE)
pooled_sd <- sqrt(ss_within / df_within)

# --- Row means, gaps, δ/σ ---
row_means <- rowMeans(M_sampled)
adj_diffs <- diff(row_means)
min_delta <- min(abs(adj_diffs))

cat("\n===============================================================================\n")
cat("FULL-DATA CELL MEANS (PM10, µg/m³)\n")
cat("===============================================================================\n\n")
print(round(M_sampled, 3))
cat(sprintf(
  "\nRow means            : %s\nAdjacent differences : %s\nMin adj gap (delta)  : %.3f\nPooled SD (sigma)    : %.3f\nDelta / sigma        : %.4f\nIUT critical value   : %.6f\n", 
  paste(sprintf("%.3f", row_means), collapse = "  "), 
  paste(sprintf("%.3f", adj_diffs), collapse = "  "), 
  min_delta, pooled_sd, min_delta / pooled_sd, c_alpha
))

# =============================================================================
# 3. DESIGN CLASS & LFC POWER EVALUATION
# =============================================================================

make_design <- function(w) {
  w <- as.numeric(w[1])
  if (length(w) != 1L || !is.finite(w) || w <= 0 || w >= 0.5) {
    stop("w must lie strictly between 0 and 0.5.")
  }
  
  # Both ends get 'w', the remaining (1 - 2w) is split equally among the middle tiers.
  p_mid <- (1 - 2 * w) / (K - 2)
  p <- c(w, rep(p_mid, K - 2), w)
  
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
  
  p <- make_design(w)
  Sigma <- Sigma_successive(p)
  
  noncentrality <- sapply(1:m, function(i) {
    standard_error <- (pooled_sd / sqrt(N_reference)) * sqrt(1 / p[i] + 1 / p[i + 1])
    delta / standard_error
  })
  
  sign_patterns <- as.matrix(expand.grid(rep(list(c(-1, 1)), m)))
  
  orthant_probabilities <- apply(sign_patterns, 1, function(sv) {
    lower <- ifelse(sv > 0, c_alpha, -Inf)
    upper <- ifelse(sv > 0, Inf, -c_alpha)
    as.numeric(pmvnorm(lower = lower, upper = upper, mean = noncentrality, sigma = Sigma, algorithm = Miwa(steps = 512)))
  })
  
  min(max(sum(orthant_probabilities), 0), 1)
}

optimize_w_iut <- function(delta, n_coarse = 80) {
  w_grid <- seq(0.01, 0.49, length.out = n_coarse)
  vals <- sapply(w_grid, function(w) power_IUT_successive(w, delta))
  best <- which.max(vals)
  
  lo <- w_grid[max(1, best - 1)]
  hi <- w_grid[min(n_coarse, best + 1)]
  opt <- optimize(function(w) power_IUT_successive(w, delta), interval = c(lo, hi), maximum = TRUE)
  
  if (opt$objective >= vals[best]) {
    list(w = opt$maximum, power = opt$objective)
  } else {
    list(w = w_grid[best], power = vals[best])
  }
}

cat("\n===============================================================================\n")
cat("MAX-MIN DESIGN OPTIMIZATION\n")
cat("===============================================================================\n")

res <- optimize_w_iut(min_delta)
w_opt <- res$w
p_opt <- make_design(w_opt)
p_bal <- rep(1 / K, K)

cat(sprintf(
  "Optimal w (w_opt)       : %.6f\nLFC power at w_opt      : %.6f\nRow proportions p_opt   : %s\nRow proportions p_bal   : %s\n", 
  w_opt, res$power, 
  paste(sprintf("%.6f", p_opt), collapse = "  "),
  paste(sprintf("%.6f", p_bal), collapse = "  ")
))

# =============================================================================
# 4. ROBUST INTEGER ALLOCATION
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
N_bal <- allocate_integers(p_bal, N_total)

cat("\n===============================================================================\n")
cat("EXACT INTEGER DESIGNS AT N =", N_total, "\n")
cat("===============================================================================\n")

cat("\nMax-min (cell-wise rounding):\n")
print(N_opt)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_opt), collapse = "  "), sum(N_opt)))

cat("\nBalanced (cell-wise rounding):\n")
print(N_bal)
cat(sprintf("  Row totals: %s  (sum = %d)\n", paste(rowSums(N_bal), collapse = "  "), sum(N_bal)))

# =============================================================================
# 5. EMPIRICAL POWER SIMULATION
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  ctr <- 0L
  row_vars <- rowSums(sd^2 / N_mat) / (C^2)
  diff_se <- sqrt(row_vars[1:m] + row_vars[2:K])
  
  for (b in 1:n_sim) {
    Ybar <- matrix(NA_real_, nrow = K, ncol = C)
    for (i in 1:K) {
      for (j in 1:C) {
        Ybar[i, j] <- mean(rnorm(N_mat[i, j], mean = M[i, j], sd = sd))
      }
    }
    
    Z <- diff(rowMeans(Ybar)) / diff_se
    
    # IUT Rejection logic: ALL Z-statistics must exceed critical value
    if (all(abs(Z) > c_alpha)) {
      ctr <- ctr + 1L
    }
  }
  return(ctr / n_sim)
}

power_opt <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim)
power_bal <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim)

se_opt <- sqrt(power_opt * (1 - power_opt) / n_sim)
se_bal <- sqrt(power_bal * (1 - power_bal) / n_sim)
gain <- 100 * (power_opt - power_bal) / power_bal

cat("\n===============================================================================\n")
cat("EMPIRICAL POWER AT N =", N_total, "\n")
cat("===============================================================================\n")
cat(sprintf(
  "Simulation size : %d\nMax-min design  : %.4f (SE %.4f)\nBalanced design : %.4f (SE %.4f)\nRelative gain   : %+.2f%%\n", 
  n_sim, power_opt, se_opt, power_bal, se_bal, gain
))

# =============================================================================
# 6. SAMPLE SIZE SEARCH FOR 90% POWER
# =============================================================================

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
    pow <- simulate_power(allocate_integers(p_row, N_val), M, sd, c_alpha, 100000) 
  }
  cat(sprintf("N = %d (power %.2f%%)\n", N_val, 100 * pow))
  
  list(N = N_val, power = pow)
}

cat("\n===============================================================================\n")
cat("SAMPLE SIZE FOR 90% POWER\n")
cat("===============================================================================\n")
cat("> Max-min design:\n")
n90_opt <- find_n90(p_opt, M_sampled, pooled_sd, c_alpha)
cat("> Balanced design:\n")
n90_bal <- find_n90(p_bal, M_sampled, pooled_sd, c_alpha)

# =============================================================================
# 7. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (IUT Successive Comparison) ══════════════════════\n")
cat(sprintf(
  " Min adjacent gap (delta) : %.3f\n Pooled SD (sigma)        : %.3f\n Delta / sigma            : %.4f\n", 
  min_delta, pooled_sd, min_delta / pooled_sd
))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " Optimal w (w_opt)        : %.6f\n Max-min row proportions  : %s\n", 
  w_opt, paste(sprintf("%.6f", p_opt), collapse = " ")
))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain            : %+.2f%%\n", 
  N_total, power_opt, N_total, power_bal, gain
))
cat("───────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction vs Bal)\n", 
  n90_opt$N, n90_opt$power, n90_bal$N, n90_bal$power, 
  n90_bal$N - n90_opt$N, 100 * (n90_bal$N - n90_opt$N) / n90_bal$N
))
cat("═══════════════════════════════════════════════════════════════════════════════════════\n")
