# =============================================================================
# IUT CONTROL-VERSUS-TREATMENT DESIGN
#
# Fixed planning values: Agra is the control row, five cities are treatment rows,
# two land-use categories form the columns.
# Computational structure: 
# 1. Fixed cell means and pooled SD 
# 2. Complete two-sided IUT max-min design
# 3. Integer allocations 
# 4. Monte Carlo power comparison 
# 5. Sample size required for 90% power
# IUT rejection rule: min_i |Z_{1i}| > Phi^{-1}(1 - alpha/2).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(27092026)

# =============================================================================
# 1. FIXED INPUT & SETTINGS
# =============================================================================

df <- read.csv("cpcb_dly_aq_uttar_pradesh-2011.csv", check.names = FALSE, na.strings = c("NA", ""))
colnames(df) <- gsub("[ /]", "_", colnames(df))

cat("Columns after rename:\n")
print(colnames(df))

df$RSPM_PM10 <- as.numeric(df$RSPM_PM10)
city_col <- "City_Town_Village_Area"  
df_clean <- df %>%
  filter(!is.na(RSPM_PM10), !is.na(Type_of_Location), !is.na(.data[[city_col]]))

# --- City & location factor levels ---
factor_a_levels <- c("Agra", "Gajraula", "Noida", "Mathura")
factor_b_levels <- c("Industrial Area", "Residential, Rural and other Areas")

df_clean <- df_clean %>%
  filter(.data[[city_col]] %in% factor_a_levels) %>%
  mutate(
    Factor_A = factor(.data[[city_col]], levels = factor_a_levels),
    Factor_B = if_else(grepl("Industrial", Type_of_Location),
                       "Industrial Area", "Residential, Rural and other Areas"),
    Factor_B = factor(Factor_B, levels = factor_b_levels)
  )

# --- Cell-wise statistics from data ---
cell_stats <- df_clean %>%
  group_by(Factor_A, Factor_B) %>%
  summarise(
    n_obs = n(),
    mean  = mean(RSPM_PM10, na.rm = TRUE),
    sd    = sd(RSPM_PM10, na.rm = TRUE),
    .groups = "drop"
  )

# --- Population mean matrix ---
M_sampled <- cell_stats %>%
  select(Factor_A, Factor_B, mean) %>%
  pivot_wider(names_from = Factor_B, values_from = mean) %>%
  as.data.frame()

rownames(M_sampled) <- M_sampled$Factor_A
M_sampled <- as.matrix(M_sampled[, factor_b_levels, drop = FALSE])
M_sampled <- M_sampled[factor_a_levels, , drop = FALSE]

# --- Pooled within-cell SD (population) ---
ss_within <- sum((cell_stats$n_obs - 1) * cell_stats$sd^2, na.rm = TRUE)
df_within <- sum(cell_stats$n_obs - 1, na.rm = TRUE)
pooled_sd <- sqrt(ss_within / df_within)

# --- Model Parameters ---
K <- nrow(M_sampled)
C <- ncol(M_sampled)
m <- K - 1
N_total <- 190
alpha <- 0.05
c_alpha <- qnorm(1 - alpha / 2) # IUT two-sided critical value without Bonferroni correction
minimum_cell_count <- 2
target_power <- 0.90
maximum_N <- 4000
n_sim_power <- 100000
n_sim_coarse <- 5000
n_sim_fine <- 100000
n_sim_validation <- 100000

# =============================================================================
# 2. PLANNING STATISTICS
# =============================================================================

row_means <- rowMeans(M_sampled)
ct_gaps <- row_means[-1] - row_means[1]
min_delta <- min(abs(ct_gaps))

cat("\n===============================================================================\n")
cat("FIXED INPUT SUMMARY (IUT CONTROL VS TREATMENT)\n")
cat("===============================================================================\n")
cat("\nCell-mean matrix (PM10, µg/m³):\n")
print(round(M_sampled, 3))
cat(sprintf(
  "\nRow means              : %s\nTreatment-control gaps : %s\nMinimum absolute gap   : %.3f\nPooled SD (sigma)      : %.3f\nDelta / sigma          : %.4f\nIUT critical value     : %.6f\n", 
  paste(sprintf("%.3f", row_means), collapse = "  "), 
  paste(sprintf("%.3f", ct_gaps), collapse = "  "), 
  min_delta, pooled_sd, min_delta / pooled_sd, c_alpha
))

# =============================================================================
# 3. DESIGN 
# =============================================================================

design_CT_cells <- function(eta) {
  eta <- as.numeric(eta[1])
  if (!is.finite(eta) || eta <= 0 || eta >= 1 / C) stop("eta must lie strictly between 0 and 1/C.")
  
  treatment_cell <- (1 - C * eta) / (C * (K - 1))
  design <- matrix(c(rep(eta, C), rep(treatment_cell, (K - 1) * C)), 
                   nrow = K, ncol = C, byrow = TRUE, dimnames = dimnames(M_sampled))
  
  stopifnot(all(design > 0), abs(sum(design) - 1) < 1e-12)
  return(design)
}

design_CT_rows <- function(eta) {
  rowSums(design_CT_cells(eta))
}

# =============================================================================
# 4. IUT POWER 
# =============================================================================

power_IUT_CT_sign_class <- function(eta, delta, k_plus, N_reference = N_total, rel_tol = 1e-8, abs_tol = 1e-10) {
  eta <- as.numeric(eta[1])
  if (!is.finite(eta) || eta <= 0 || eta >= 1 / C) return(NA_real_)
  if (!is.finite(delta) || delta <= 0) stop("delta must be positive and finite.")
  if (k_plus < 0 || k_plus > m || k_plus != as.integer(k_plus)) stop("k_plus must be an integer between 0 and m.")
  
  p_row <- design_CT_rows(eta)
  p_control <- p_row[1]
  p_treatment <- p_row[2]
  
  standard_error <- (pooled_sd / sqrt(N_reference)) * sqrt(1 / p_control + 1 / p_treatment)
  rho <- (1 / p_control) / (1 / p_control + 1 / p_treatment)
  
  if (!is.finite(rho) || rho <= 0 || rho >= 1) return(NA_real_)
  
  sign_vector <- c(rep(1, k_plus), rep(-1, m - k_plus))
  noncentrality <- (delta / standard_error) * sign_vector
  sqrt_rho <- sqrt(rho)
  conditional_sd <- sqrt(1 - rho)
  
  integrand <- function(z) {
    vapply(z, function(z_value) {
      conditional_means <- noncentrality + sqrt_rho * z_value
      upper_tail <- pnorm(q = c_alpha, mean = conditional_means, sd = conditional_sd, lower.tail = FALSE)
      lower_tail <- pnorm(q = -c_alpha, mean = conditional_means, sd = conditional_sd, lower.tail = TRUE)
      component_probabilities <- pmin(pmax(upper_tail + lower_tail, 0), 1)
      prod(component_probabilities) * dnorm(z_value)
    }, numeric(1))
  }
  
  integration_result <- integrate(
    f = integrand, lower = -Inf, upper = Inf, 
    rel.tol = rel_tol, abs.tol = abs_tol, subdivisions = 1000L, stop.on.error = TRUE
  )
  min(max(integration_result$value, 0), 1)
}

# =============================================================================
# 5. Power at LFC
# =============================================================================

worst_case_IUT_CT <- function(eta, delta, N_reference = N_total) {
  sign_classes <- 0:m
  powers <- vapply(sign_classes, function(k) power_IUT_CT_sign_class(eta, delta, k, N_reference), numeric(1))
  
  if (any(!is.finite(powers))) {
    return(list(worst_power = NA_real_, worst_k_plus = NA_integer_, sign_class_powers = powers))
  }
  
  worst_index <- which.min(powers)
  list(worst_power = powers[worst_index], worst_k_plus = sign_classes[worst_index], sign_class_powers = powers)
}

# =============================================================================
# 6. CONSTRAINED MAX-MIN DESIGN
# =============================================================================

optimize_CT_IUT_design <- function(delta, N_reference = N_total, min_n = minimum_cell_count, n_grid = 301) {
  eta_lower <- min_n / N_reference
  eta_upper <- (1 - min_n * C * (K - 1) / N_reference) / C
  
  if (!is.finite(eta_lower) || !is.finite(eta_upper) || eta_upper <= eta_lower) {
    stop("No eta satisfies the minimum-cell allocation requirement.")
  }
  
  eta_grid <- seq(eta_lower, eta_upper, length.out = n_grid)
  grid_results <- lapply(eta_grid, function(e) worst_case_IUT_CT(e, delta, N_reference))
  grid_power <- vapply(grid_results, function(r) r$worst_power, numeric(1))
  
  if (any(!is.finite(grid_power))) stop("Non-finite worst-case powers occurred during optimization.")
  
  best_index <- which.max(grid_power)
  lower_index <- max(1L, best_index - 1L)
  upper_index <- min(n_grid, best_index + 1L)
  local_lower <- eta_grid[lower_index]
  local_upper <- eta_grid[upper_index]
  
  if (local_lower < local_upper) {
    refined_result <- optimize(
      f = function(e) worst_case_IUT_CT(e, delta, N_reference)$worst_power, 
      interval = c(local_lower, local_upper), maximum = TRUE, tol = 1e-7
    )
    if (is.finite(refined_result$objective) && refined_result$objective >= grid_power[best_index]) {
      eta_opt <- refined_result$maximum
      max_min_power <- refined_result$objective
    } else {
      eta_opt <- eta_grid[best_index]
      max_min_power <- grid_power[best_index]
    }
  } else {
    eta_opt <- eta_grid[best_index]
    max_min_power <- grid_power[best_index]
  }
  
  worst_result <- worst_case_IUT_CT(eta_opt, delta, N_reference)
  boundary_solution <- (abs(eta_opt - eta_lower) <= 1e-5 || abs(eta_opt - eta_upper) <= 1e-5)
  cell_design <- design_CT_cells(eta_opt)
  
  list(
    eta_opt = eta_opt, max_min_power = max_min_power, optimal_cell_design = cell_design, 
    optimal_row_design = rowSums(cell_design), worst_k_plus = worst_result$worst_k_plus, 
    sign_class_powers = worst_result$sign_class_powers, eta_lower = eta_lower, 
    eta_upper = eta_upper, boundary_solution = boundary_solution, 
    eta_grid = eta_grid, power_grid = grid_power
  )
}

cat("\n===============================================================================\n")
cat("MAX-MIN DESIGN OPTIMIZATION\n")
cat("===============================================================================\n")

design_result <- optimize_CT_IUT_design(delta = min_delta, N_reference = N_total, min_n = minimum_cell_count, n_grid = 301)
eta_opt <- design_result$eta_opt
xi_opt <- design_result$optimal_cell_design
p_opt <- design_result$optimal_row_design

xi_bal <- matrix(1 / (K * C), nrow = K, ncol = C, dimnames = dimnames(M_sampled))
p_bal <- rowSums(xi_bal)

eta_dunnett <- 1 / (C * (1 + sqrt(K - 1)))
xi_dunnett <- design_CT_cells(eta_dunnett)
p_dunnett <- rowSums(xi_dunnett)

cat(sprintf(
  "Admissible eta interval : [%.8f, %.8f]\nOptimal eta             : %.8f\nOptimal row proportions : %s\nSum of row proportions  : %.12f\nWorst-case IUT power    : %.8f\nWorst sign class k_plus : %d of %d\nSign-class powers       : %s\nBoundary solution       : %s\nDunnett limiting eta    : %.8f\nDunnett row proportions : %s\n", 
  design_result$eta_lower, design_result$eta_upper, eta_opt, 
  paste(sprintf("%.6f", p_opt), collapse = " "), sum(p_opt), 
  design_result$max_min_power, design_result$worst_k_plus, m, 
  paste(sprintf("%d:%.6f", 0:m, design_result$sign_class_powers), collapse = "  "), 
  design_result$boundary_solution, eta_dunnett, 
  paste(sprintf("%.6f", p_dunnett), collapse = " ")
))

# =============================================================================
# 7. INTEGER ALLOCATION
# =============================================================================

allocate_integers <- function(target_matrix, total_N, seed = 26092026) {
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
  
  # --- Base allocation ---
  N_mat <- floor(target_matrix)
  rem <- total_N - sum(N_mat)
  
  if (rem > 0) {
    set.seed(seed)   # temporary local stream
    frac <- as.vector(target_matrix - floor(target_matrix))
    
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
  
  dimnames(N_mat) <- dimnames(target_matrix)
  return(N_mat)
}

N_opt <- allocate_integers(N_total * xi_opt, N_total)
N_dun <- allocate_integers(N_total * xi_dunnett, N_total)
N_bal <- allocate_integers(N_total * xi_bal, N_total)

cat("\n===============================================================================\n")
cat("INTEGER ALLOCATIONS AT N =", N_total, "\n")
cat("===============================================================================\n")

cat("\nMax-min target counts:\n")
print(round(N_total * xi_opt, 6))
cat("\nMax-min integer allocation:\n")
print(N_opt)
cat(sprintf("Row totals: %s; total = %d\n\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt)))

cat("Dunnett integer allocation:\n")
print(N_dun)
cat(sprintf("Row totals: %s; total = %d\n\n", paste(rowSums(N_dun), collapse = " "), sum(N_dun)))

cat("Balanced integer allocation:\n")
print(N_bal)
cat(sprintf("Row totals: %s; total = %d\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))

# =============================================================================
# 8. EMPIRICAL POWER
# =============================================================================

simulate_power <- function(N_mat, M, sd, c_alpha, n_sim) {
  K_local <- nrow(N_mat)
  C_local <- ncol(N_mat)
  
  if (!all(dim(N_mat) == dim(M))) stop("N_mat and M must have identical dimensions.")
  if (any(N_mat <= 0)) stop("All cell counts must be positive.")
  
  row_variances <- rowSums(sd^2 / N_mat) / C_local^2
  contrast_standard_errors <- sqrt(row_variances[1] + row_variances[2:K_local])
  rejection_count <- 0L
  
  for (simulation_index in seq_len(n_sim)) {
    simulated_cell_means <- matrix(NA_real_, nrow = K_local, ncol = C_local)
    for (i in seq_len(K_local)) {
      for (j in seq_len(C_local)) {
        simulated_cell_means[i, j] <- rnorm(1, mean = M[i, j], sd = sd / sqrt(N_mat[i, j]))
      }
    }
    
    simulated_row_means <- rowMeans(simulated_cell_means)
    Z_statistics <- (simulated_row_means[1] - simulated_row_means[2:K_local]) / contrast_standard_errors
    
    # IUT Rejection logic: ALL Z-statistics must exceed critical value
    if (all(abs(Z_statistics) > c_alpha)) rejection_count <- rejection_count + 1L
  }
  
  power <- rejection_count / n_sim
  list(power = power, monte_carlo_se = sqrt(power * (1 - power) / n_sim))
}

power_opt_result <- simulate_power(N_opt, M_sampled, pooled_sd, c_alpha, n_sim_power)
power_dun_result <- simulate_power(N_dun, M_sampled, pooled_sd, c_alpha, n_sim_power)
power_bal_result <- simulate_power(N_bal, M_sampled, pooled_sd, c_alpha, n_sim_power)

power_opt <- power_opt_result$power
power_dun <- power_dun_result$power
power_bal <- power_bal_result$power
relative_gain <- 100 * (power_opt - power_bal) / power_bal

cat("\n===============================================================================\n")
cat("EMPIRICAL POWER AT N =", N_total, "\n")
cat("===============================================================================\n")
cat(sprintf(
  "Max-min design  : %.4f (MC SE %.5f)\nDunnett design  : %.4f (MC SE %.5f)\nBalanced design : %.4f (MC SE %.5f)\nRelative gain (max-min vs balanced): %.2f%%\n", 
  power_opt, power_opt_result$monte_carlo_se, 
  power_dun, power_dun_result$monte_carlo_se, 
  power_bal, power_bal_result$monte_carlo_se, 
  relative_gain
))

# =============================================================================
# 9. SAMPLE SIZE REQUIRED FOR 90% POWER
# =============================================================================

evaluate_power_at_N <- function(N_candidate, cell_proportions, simulations) {
  allocation <- allocate_integers(N_candidate * cell_proportions, N_candidate)
  result <- simulate_power(allocation, M_sampled, pooled_sd, c_alpha, simulations)
  list(N = N_candidate, allocation = allocation, power = result$power, monte_carlo_se = result$monte_carlo_se)
}

find_n90 <- function(cell_proportions, design_label, target = target_power, N_start = N_total, N_max = maximum_N, coarse_step = 10) {
  cat(sprintf("\nTesting %s:\n", design_label))
  N_current <- N_start
  
  cat("  [Coarse step 10] ... ")
  coarse_result <- evaluate_power_at_N(N_current, cell_proportions, n_sim_coarse)
  
  while (coarse_result$power < target) {
    N_current <- N_current + coarse_step
    if (N_current > N_max) stop(paste("Target power was not reached by N =", N_max))
    coarse_result <- evaluate_power_at_N(N_current, cell_proportions, n_sim_coarse)
  }
  cat(sprintf("reached ~%d.\n  [Fine step 1] ... ", N_current))
  
  fine_lower <- max(N_start, N_current - 2 * coarse_step)
  fine_results <- lapply(seq(fine_lower, N_current, by = 1), function(n) evaluate_power_at_N(n, cell_proportions, n_sim_fine))
  fine_powers <- vapply(fine_results, function(r) r$power, numeric(1))
  
  feasible_indices <- which(fine_powers >= target)
  if (length(feasible_indices) == 0L) stop("No candidate in the fine-search interval attained target power.")
  selected_result <- fine_results[[min(feasible_indices)]]
  
  validation_result <- evaluate_power_at_N(selected_result$N, cell_proportions, n_sim_validation)
  
  while (validation_result$power < target) {
    next_N <- validation_result$N + 1L
    if (next_N > N_max) stop("Target power was not reached during validation.")
    validation_result <- evaluate_power_at_N(next_N, cell_proportions, n_sim_validation)
  }
  
  repeat {
    previous_N <- validation_result$N - 1L
    if (previous_N < N_start) break
    previous_result <- evaluate_power_at_N(previous_N, cell_proportions, n_sim_validation)
    if (previous_result$power >= target) {
      validation_result <- previous_result 
    } else {
      break
    }
  }
  
  cat(sprintf("Validated N = %d (power %.4f, MC SE = %.5f)\n", validation_result$N, validation_result$power, validation_result$monte_carlo_se))
  return(validation_result)
}

cat("\n===============================================================================\n")
cat("SAMPLE SIZE REQUIRED FOR 90% POWER\n")
cat("===============================================================================\n")
n90_opt <- find_n90(xi_opt, "Max-min design")
n90_dun <- find_n90(xi_dunnett, "Dunnett design")
n90_bal <- find_n90(xi_bal, "Balanced design")

sample_savings <- n90_bal$N - n90_opt$N
sample_reduction_percent <- 100 * sample_savings / n90_bal$N

# =============================================================================
# 10. FINAL SUMMARY
# =============================================================================

cat("\n══════════════════════ FINAL SUMMARY (IUT Control vs Treatment) ══════════════════════\n")
cat(sprintf(
  " Minimum CT gap (delta)   : %.3f\n Pooled SD (sigma)        : %.3f\n Delta / sigma            : %.4f\n", 
  min_delta, pooled_sd, min_delta / pooled_sd
))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " Optimal eta              : %.6f\n Constrained max-min prop : %s\n Worst sign class k_plus  : %d\n", 
  eta_opt, paste(sprintf("%.6f", p_opt), collapse = " "), design_result$worst_k_plus
))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " Power @ N=%d  Max-min     : %.4f\n Power @ N=%d  Dunnett     : %.4f\n Power @ N=%d  Balanced    : %.4f\n Relative gain (vs Bal)   : %+.2f%%\n", 
  N_total, power_opt, N_total, power_dun, N_total, power_bal, relative_gain
))
cat("──────────────────────────────────────────────────────────────────────────────────────\n")
cat(sprintf(
  " 90%% N  Max-min           : %d (power %.4f)\n 90%% N  Dunnett           : %d (power %.4f)\n 90%% N  Balanced          : %d (power %.4f)\n Sample size savings      : %d (%.2f%% reduction vs Bal)\n", 
  n90_opt$N, n90_opt$power, n90_dun$N, n90_dun$power, n90_bal$N, n90_bal$power, 
  sample_savings, sample_reduction_percent
))
cat("══════════════════════════════════════════════════════════════════════════════════════\n")
