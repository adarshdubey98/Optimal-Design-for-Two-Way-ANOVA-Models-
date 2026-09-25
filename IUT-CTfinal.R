# =============================================================================
# IUT CONTROL-VERSUS-TREATMENT DESIGN
#
# Fixed planning values: Agra is the control row, five cities are treatment rows,
# two land-use categories form the columns.
# Computational structure: 1. Fixed cell means and pooled SD, 2. Complete two-sided 
# IUT max-min design, 3. Integer allocations, 4. Monte Carlo power comparison, 
# 5. Sample size required for 90% power
# IUT rejection rule: min_i |Z_{1i}| > Phi^{-1}(1 - alpha/2).
# =============================================================================

rm(list = ls())

# =============================================================================
# 1. FIXED INPUT & 2. SETTINGS
# =============================================================================

M_sampled <- matrix(c(
  191.125, 173.71429, # Agra (control)
  92.500,  86.28571,  # Gajraula
  69.250,  62.42857,  # Gorakhpur
  172.500, 151.14286, # Khurja
  160.150, 120.64286, # Rai Bareilly
  137.375, 135.28571  # Unnao
), nrow = 6, ncol = 2, byrow = TRUE)

rownames(M_sampled) <- c("Agra (control)", "Gajraula", "Gorakhpur", "Khurja", "Rai Bareilly", "Unnao")
colnames(M_sampled) <- c("Industrial Area", "Residential, Rural and other Areas")
pooled_sd <- 34.7

K <- nrow(M_sampled); C <- ncol(M_sampled); m <- K - 1
N_total <- 90; alpha <- 0.05; c_alpha <- qnorm(1 - alpha / 2)
minimum_cell_count <- 2; target_power <- 0.90; maximum_N <- 4000
n_sim_power <- 100000; n_sim_coarse <- 5000; n_sim_fine <- 20000; n_sim_validation <- 100000

# =============================================================================
# 3. PLANNING STATISTICS
# =============================================================================

row_means <- rowMeans(M_sampled)
ct_gaps <- row_means[-1] - row_means[1]
min_delta <- min(abs(ct_gaps))

cat("============================================================\nFIXED INPUT SUMMARY\n============================================================\n")
cat("\nCell-mean matrix:\n"); print(round(M_sampled, 3))
cat(sprintf("\nRow means              : %s\n", paste(sprintf("%.3f", row_means), collapse = "  ")))
cat(sprintf("Treatment-control gaps: %s\n", paste(sprintf("%.3f", ct_gaps), collapse = "  ")))
cat(sprintf("Minimum absolute gap  : %.3f\n", min_delta))
cat(sprintf("Pooled SD              : %.3f\n", pooled_sd))
cat(sprintf("Delta / sigma          : %.4f\n", min_delta / pooled_sd))
cat(sprintf("IUT critical value     : %.6f\n", c_alpha))

# =============================================================================
# 4. CONTROL-VERSUS-TREATMENT DESIGN CLASS
# =============================================================================

design_CT_cells <- function(eta) {
  eta <- as.numeric(eta[1])
  if (!is.finite(eta) || eta <= 0 || eta >= 1 / C) stop("eta must lie strictly between 0 and 1/C.")
  treatment_cell <- (1 - C * eta) / (C * (K - 1))
  design <- matrix(c(rep(eta, C), rep(treatment_cell, (K - 1) * C)), nrow = K, ncol = C, byrow = TRUE, dimnames = dimnames(M_sampled))
  stopifnot(all(design > 0), abs(sum(design) - 1) < 1e-12)
  return(design)
}

design_CT_rows <- function(eta) rowSums(design_CT_cells(eta))

# =============================================================================
# 5. COMPLETE TWO-SIDED IUT POWER AT A SIGN CLASS
# =============================================================================

power_IUT_CT_sign_class <- function(eta, delta, k_plus, N_reference = N_total, relative_tolerance = 1e-8, absolute_tolerance = 1e-10) {
  eta <- as.numeric(eta[1])
  if (!is.finite(eta) || eta <= 0 || eta >= 1 / C) return(NA_real_)
  if (!is.finite(delta) || delta <= 0) stop("delta must be positive and finite.")
  if (k_plus < 0 || k_plus > m || k_plus != as.integer(k_plus)) stop("k_plus must be an integer between 0 and m.")
  
  p_row <- design_CT_rows(eta); p_control <- p_row[1]; p_treatment <- p_row[2]
  standard_error <- (pooled_sd / sqrt(N_reference)) * sqrt(1 / p_control + 1 / p_treatment)
  rho <- (1 / p_control) / (1 / p_control + 1 / p_treatment)
  if (!is.finite(rho) || rho <= 0 || rho >= 1) return(NA_real_)
  
  sign_vector <- c(rep(1, k_plus), rep(-1, m - k_plus))
  noncentrality <- (delta / standard_error) * sign_vector
  sqrt_rho <- sqrt(rho); conditional_sd <- sqrt(1 - rho)
  
  integrand <- function(z) {
    vapply(z, function(z_value) {
      conditional_means <- noncentrality + sqrt_rho * z_value
      upper_tail <- pnorm(q = c_alpha, mean = conditional_means, sd = conditional_sd, lower.tail = FALSE)
      lower_tail <- pnorm(q = -c_alpha, mean = conditional_means, sd = conditional_sd, lower.tail = TRUE)
      component_probabilities <- pmin(pmax(upper_tail + lower_tail, 0), 1)
      prod(component_probabilities) * dnorm(z_value)
    }, numeric(1))
  }
  
  integration_result <- integrate(f = integrand, lower = -Inf, upper = Inf, rel.tol = relative_tolerance, abs.tol = absolute_tolerance, subdivisions = 1000L, stop.on.error = TRUE)
  min(max(integration_result$value, 0), 1)
}

# =============================================================================
# 6. WORST-CASE POWER OVER SIGN CLASSES
# =============================================================================

worst_case_IUT_CT <- function(eta, delta, N_reference = N_total) {
  sign_classes <- 0:m
  powers <- vapply(sign_classes, function(k) power_IUT_CT_sign_class(eta, delta, k, N_reference), numeric(1))
  if (any(!is.finite(powers))) return(list(worst_power = NA_real_, worst_k_plus = NA_integer_, sign_class_powers = powers))
  worst_index <- which.min(powers)
  list(worst_power = powers[worst_index], worst_k_plus = sign_classes[worst_index], sign_class_powers = powers)
}

# =============================================================================
# 7. CONSTRAINED MAX-MIN DESIGN
# =============================================================================

optimize_CT_IUT_design <- function(delta, N_reference = N_total, min_n = minimum_cell_count, n_grid = 301) {
  eta_lower <- min_n / N_reference
  eta_upper <- (1 - min_n * C * (K - 1) / N_reference) / C
  if (!is.finite(eta_lower) || !is.finite(eta_upper) || eta_upper <= eta_lower) stop("No eta satisfies the minimum-cell allocation requirement.")
  
  eta_grid <- seq(eta_lower, eta_upper, length.out = n_grid)
  grid_results <- lapply(eta_grid, function(e) worst_case_IUT_CT(e, delta, N_reference))
  grid_power <- vapply(grid_results, function(r) r$worst_power, numeric(1))
  if (any(!is.finite(grid_power))) stop("Non-finite worst-case powers occurred during optimization.")
  
  best_index <- which.max(grid_power)
  lower_index <- max(1L, best_index - 1L); upper_index <- min(n_grid, best_index + 1L)
  local_lower <- eta_grid[lower_index]; local_upper <- eta_grid[upper_index]
  
  if (local_lower < local_upper) {
    refined_result <- optimize(f = function(e) worst_case_IUT_CT(e, delta, N_reference)$worst_power, interval = c(local_lower, local_upper), maximum = TRUE, tol = 1e-7)
    if (is.finite(refined_result$objective) && refined_result$objective >= grid_power[best_index]) {
      eta_opt <- refined_result$maximum; max_min_power <- refined_result$objective
    } else {
      eta_opt <- eta_grid[best_index]; max_min_power <- grid_power[best_index]
    }
  } else {
    eta_opt <- eta_grid[best_index]; max_min_power <- grid_power[best_index]
  }
  
  worst_result <- worst_case_IUT_CT(eta_opt, delta, N_reference)
  boundary_solution <- (abs(eta_opt - eta_lower) <= 1e-5 || abs(eta_opt - eta_upper) <= 1e-5)
  cell_design <- design_CT_cells(eta_opt)
  
  list(eta_opt = eta_opt, max_min_power = max_min_power, optimal_cell_design = cell_design, optimal_row_design = rowSums(cell_design), worst_k_plus = worst_result$worst_k_plus, sign_class_powers = worst_result$sign_class_powers, eta_lower = eta_lower, eta_upper = eta_upper, boundary_solution = boundary_solution, eta_grid = eta_grid, power_grid = grid_power)
}

cat("\n============================================================\nMAX-MIN DESIGN OPTIMIZATION\n============================================================\n")
design_result <- optimize_CT_IUT_design(delta = min_delta, N_reference = N_total, min_n = minimum_cell_count, n_grid = 301)
eta_opt <- design_result$eta_opt; xi_opt <- design_result$optimal_cell_design; p_opt <- design_result$optimal_row_design

xi_bal <- matrix(1 / (K * C), nrow = K, ncol = C, dimnames = dimnames(M_sampled)); p_bal <- rowSums(xi_bal)
eta_dunnett <- 1 / (C * (1 + sqrt(K - 1)))
xi_dunnett <- design_CT_cells(eta_dunnett); p_dunnett <- rowSums(xi_dunnett)

cat(sprintf("Admissible eta interval : [%.8f, %.8f]\n", design_result$eta_lower, design_result$eta_upper))
cat(sprintf("Optimal eta             : %.8f\n", eta_opt))
cat("Optimal row proportions :", paste(sprintf("%.6f", p_opt), collapse = " "), "\n")
cat(sprintf("Sum of row proportions  : %.12f\nWorst-case IUT power    : %.8f\n", sum(p_opt), design_result$max_min_power))
cat(sprintf("Worst sign class k_plus : %d of %d\n", design_result$worst_k_plus, m))
cat("Sign-class powers       :", paste(sprintf("%d:%.6f", 0:m, design_result$sign_class_powers), collapse = "  "), "\n")
cat(sprintf("Boundary solution       : %s\n\nDunnett limiting eta    : %.8f\n", design_result$boundary_solution, eta_dunnett))
cat("Dunnett row proportions :", paste(sprintf("%.6f", p_dunnett), collapse = " "), "\n")

# =============================================================================
# 8. INTEGER ALLOCATION
# =============================================================================

# allocate_integers <- function(target_matrix, total_N) {
#   target_matrix <- as.matrix(target_matrix)
#   
#   if (any(!is.finite(target_matrix)) || any(target_matrix < 0)) {
#     stop("target_matrix must contain finite nonnegative values.")
#   }
#   
#   if (!is.finite(total_N) || total_N < 0 || total_N != as.integer(total_N)) {
#     stop("total_N must be a nonnegative integer.")
#   }
#   
#   target_sum <- sum(target_matrix)
#   if (target_sum <= 0) {
#     stop("The sum of target_matrix must be positive.")
#   }
#   
#   # Rescale the targets to ensure that they sum exactly to total_N.
#   target_counts <- target_matrix * total_N / target_sum
#   
#   # Initial integer allocation
#   allocation <- floor(target_counts)
#   
#   # Number of observations still to be allocated
#   remaining <- total_N - sum(allocation)
#   
#   if (remaining > 0L) {
#     fractional_parts <- as.vector(target_counts - floor(target_counts))
#     selected_cells <- order(fractional_parts, decreasing = TRUE, method = "radix")[seq_len(remaining)]
#     allocation[selected_cells] <- allocation[selected_cells] + 1L
#   }
#   
#   allocation <- matrix(allocation, nrow = nrow(target_matrix), ncol = ncol(target_matrix), dimnames = dimnames(target_matrix))
#   
#   stopifnot(sum(allocation) == total_N, all(allocation >= 0))
#   
#   return(allocation)
# }
allocate_integers <- function(target_matrix, total_N) {
  # 1. Base allocation (sirf integer part)
  N_mat <- floor(target_matrix)
  
  # 2. Bachi hui observations (remainder)
  rem <- total_N - sum(N_mat)
  
  # 3. Agar remainder bacha hai, toh sabse bade fraction walo ko de do
  if (rem > 0) {
    frac <- as.vector(target_matrix - floor(target_matrix))
    idx <- order(frac, decreasing = TRUE, method = "radix")[seq_len(rem)]
    N_mat[idx] <- N_mat[idx] + 1L
  }
  
  dimnames(N_mat) <- dimnames(target_matrix)
  return(N_mat)
}

N_opt <- allocate_integers(N_total * xi_opt, N_total)
N_dun <- allocate_integers(N_total * xi_dunnett, N_total)
N_bal <- allocate_integers(N_total * xi_bal, N_total)

cat("\n============================================================\nINTEGER ALLOCATIONS AT N = 90\n============================================================\n")
cat("\nMax-min target counts:\n"); print(round(N_total * xi_opt, 6))
cat("\nMax-min integer allocation:\n"); print(N_opt)
cat(sprintf("Row totals: %s; total = %d\n\nDunnett integer allocation:\n", paste(rowSums(N_opt), collapse = " "), sum(N_opt)))
print(N_dun)
cat(sprintf("Row totals: %s; total = %d\n\nBalanced integer allocation:\n", paste(rowSums(N_dun), collapse = " "), sum(N_dun)))
print(N_bal)
cat(sprintf("Row totals: %s; total = %d\n", paste(rowSums(N_bal), collapse = " "), sum(N_bal)))

# =============================================================================
# 9. EMPIRICAL POWER
# =============================================================================

simulate_power <- function(allocation, mean_matrix, sigma, simulations) {
  K_local <- nrow(allocation); C_local <- ncol(allocation)
  if (!all(dim(allocation) == dim(mean_matrix))) stop("allocation and mean_matrix must have identical dimensions.")
  if (any(allocation <= 0)) stop("All cell counts must be positive.")
  
  row_variances <- rowSums(sigma^2 / allocation) / C_local^2
  contrast_standard_errors <- sqrt(row_variances[1] + row_variances[2:K_local])
  rejection_count <- 0L
  
  for (simulation_index in seq_len(simulations)) {
    simulated_cell_means <- matrix(NA_real_, nrow = K_local, ncol = C_local)
    for (i in seq_len(K_local)) {
      for (j in seq_len(C_local)) {
        simulated_cell_means[i, j] <- rnorm(1, mean = mean_matrix[i, j], sd = sigma / sqrt(allocation[i, j]))
      }
    }
    simulated_row_means <- rowMeans(simulated_cell_means)
    Z_statistics <- (simulated_row_means[1] - simulated_row_means[2:K_local]) / contrast_standard_errors
    if (all(abs(Z_statistics) > c_alpha)) rejection_count <- rejection_count + 1L
  }
  power <- rejection_count / simulations
  list(power = power, monte_carlo_se = sqrt(power * (1 - power) / simulations))
}

power_opt_result <- simulate_power(N_opt, M_sampled, pooled_sd, n_sim_power)
power_dun_result <- simulate_power(N_dun, M_sampled, pooled_sd, n_sim_power)
power_bal_result <- simulate_power(N_bal, M_sampled, pooled_sd, n_sim_power)

power_opt <- power_opt_result$power; power_dun <- power_dun_result$power; power_bal <- power_bal_result$power
relative_gain <- 100 * (power_opt - power_bal) / power_bal

cat("\n============================================================\nEMPIRICAL POWER AT N = 90\n============================================================\n")
cat(sprintf("Max-min design  : %.4f (MC SE %.5f)\n", power_opt, power_opt_result$monte_carlo_se))
cat(sprintf("Dunnett design  : %.4f (MC SE %.5f)\n", power_dun, power_dun_result$monte_carlo_se))
cat(sprintf("Balanced design : %.4f (MC SE %.5f)\n", power_bal, power_bal_result$monte_carlo_se))
cat(sprintf("Relative gain, max-min vs balanced: %.2f%%\n", relative_gain))

# =============================================================================
# 10. SAMPLE SIZE REQUIRED FOR 90% POWER
# =============================================================================

evaluate_power_at_N <- function(N_candidate, cell_proportions, simulations) {
  allocation <- allocate_integers(N_candidate * cell_proportions, N_candidate)
  result <- simulate_power(allocation, M_sampled, pooled_sd, simulations)
  list(N = N_candidate, allocation = allocation, power = result$power, monte_carlo_se = result$monte_carlo_se)
}

find_sample_size <- function(cell_proportions, design_label, target = target_power, N_start = N_total, N_max = maximum_N, coarse_step = 10) {
  cat(sprintf("\nTesting %s:\n", design_label))
  N_current <- N_start
  coarse_result <- evaluate_power_at_N(N_current, cell_proportions, n_sim_coarse)
  
  while (coarse_result$power < target) {
    N_current <- N_current + coarse_step
    if (N_current > N_max) stop(paste("Target power was not reached by N =", N_max))
    coarse_result <- evaluate_power_at_N(N_current, cell_proportions, n_sim_coarse)
  }
  cat(sprintf("  Coarse search crossed the target near N = %d.\n", N_current))
  
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
    if (previous_result$power >= target) validation_result <- previous_result else break
  }
  
  cat(sprintf("  Validated N = %d, power = %.4f, MC SE = %.5f.\n", validation_result$N, validation_result$power, validation_result$monte_carlo_se))
  validation_result
}

cat("\n============================================================\nSAMPLE SIZE REQUIRED FOR 90% POWER\n============================================================\n")
n90_opt <- find_sample_size(xi_opt, "fixed minimum-cell-constrained max-min design")
n90_dun <- find_sample_size(xi_dunnett, "fixed Dunnett limiting design")
n90_bal <- find_sample_size(xi_bal, "balanced design")

sample_savings <- n90_bal$N - n90_opt$N
sample_reduction_percent <- 100 * sample_savings / n90_bal$N

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n===============================================================================\nFINAL SUMMARY: IUT CONTROL VERSUS TREATMENT\n===============================================================================\n")
cat(sprintf("Minimum CT gap                 : %.3f\nPooled SD                      : %.3f\nDelta / sigma                  : %.4f\n", min_delta, pooled_sd, min_delta / pooled_sd))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("Optimal eta                    : %.6f\nConstrained max-min row prop   : %s\nBoundary solution              : %s\nWorst sign class, k_plus       : %d\n", eta_opt, paste(sprintf("%.6f", p_opt), collapse = " "), design_result$boundary_solution, design_result$worst_k_plus))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("Power at N=%d, max-min         : %.4f\nPower at N=%d, Dunnett         : %.4f\nPower at N=%d, balanced        : %.4f\nRelative gain vs balanced      : %.2f%%\n", N_total, power_opt, N_total, power_dun, N_total, power_bal, relative_gain))
cat("-------------------------------------------------------------------------------\n")
cat(sprintf("90%% power N, max-min           : %d (power %.4f)\n90%% power N, Dunnett           : %d (power %.4f)\n90%% power N, balanced          : %d (power %.4f)\nSample-size saving vs balanced : %d (%.2f%% reduction)\n", n90_opt$N, n90_opt$power, n90_dun$N, n90_dun$power, n90_bal$N, n90_bal$power, sample_savings, sample_reduction_percent))
cat("-------------------------------------------------------------------------------\n\nMax-min integer allocation at N = 90:\n")
print(N_opt)
cat("\nDunnett integer allocation at N = 90:\n")
print(N_dun)
cat("\nBalanced integer allocation at N = 90:\n")
print(N_bal)
cat("===============================================================================\n")