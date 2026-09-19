# =========================================================================
# MULTI-TOPOLOGY & MULTI-PARAMETRIC DENSITY OUTLIER SIMULATION
# Grid Expansion: Sample Sizes n in {100, 500} | Contamination r_0 in {1, 3, 5, 10}
# Output: Exactly 16 Superimposed Plots (2 per Topology: Separability & ROC)
# =========================================================================

library(ggplot2)
library(stats)
library(pROC)
library(gridExtra)

# --- 1. CORE EXPERIMENTAL GRID CONFIGURATION ---
set.seed(2026)

sample_sizes   <- c(100, 500)
r_true_grid    <- c(1, 3, 5, 10)
replications   <- 200
bootstrap_reps <- 250

topologies <- c("Normal", "Exponential", "Cauchy", "Beta", 
                "Gamma", "Weibull", "Student_t", "Uniform")

main_plot_dir <- "manuscript_plots"
if (!dir.exists(main_plot_dir)) dir.create(main_plot_dir)

summary_results <- data.frame(
  Topology        = character(),
  n_Sample        = integer(),
  r_True          = integer(),
  FPR             = numeric(),
  TPR             = numeric(),
  AUC             = numeric(),
  Mean_r_hat_H0   = numeric(),
  Mean_r_hat_H1   = numeric(),
  r_hat_Accuracy  = numeric(),
  stringsAsFactors = FALSE
)

# --- 2. ALGORITHMIC FUNCTIONS ---

get_robust_bandwidth <- function(x) {
  n <- length(x)
  sigma_robust <- median(abs(x - median(x))) * 1.4826
  if (sigma_robust == 0) sigma_robust <- sd(x)
  h <- 0.9 * sigma_robust * (n^(-1/5))
  return(max(h, 1e-5))
}

estimate_r_hat <- function(densities, n, h, r_max = NULL) {
  if (is.null(r_max)) r_max <- max(1, floor(5 * log(n)))
  
  sorted_densities <- sort(densities, decreasing = FALSE)
  sorted_densities[sorted_densities <= 0] <- .Machine$double.eps
  log_dens <- log(sorted_densities)
  
  C_r <- 1.0
  tau_n <- C_r * sqrt(log(n) / (n * h))
  
  r_hat <- r_max
  for (k in 1:(r_max - 1)) {
    gap <- log_dens[k + 1] - log_dens[k]
    if (gap < tau_n) {
      r_hat <- k
      break
    }
  }
  return(r_hat)
}

compute_t_statistic <- function(x) {
  n <- length(x)
  h <- get_robust_bandwidth(x)
  
  densities <- numeric(n)
  for (i in 1:n) {
    densities[i] <- mean(dnorm((x[i] - x) / h)) / h
  }
  
  r_hat <- estimate_r_hat(densities, n, h)
  
  sorted_densities <- sort(densities, decreasing = FALSE)
  isolated_points <- sorted_densities[1:r_hat]
  isolated_points[isolated_points <= 0] <- .Machine$double.eps
  
  stat_value <- -mean(log(isolated_points))
  return(list(stat = stat_value, r_hat = r_hat))
}

# --- 3. EXECUTION LOOP WITH PLOT DATA ACCUMULATION ---

for (topo in topologies) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf(" Processing Topology: [%s]\n", topo))
  cat(sprintf("==================================================\n"))
  
  # Data accumulators for superimposing all grid settings of this topology
  topo_density_data <- data.frame()
  topo_roc_data     <- data.frame()
  
  for (n_curr in sample_sizes) {
    for (r_curr in r_true_grid) {
      
      if (r_curr >= n_curr / 2) next 
      
      cat(sprintf(" Grid Config -> Topology: %-11s | n = %3d | r_0 = %2d ... ", 
                  topo, n_curr, r_curr))
      
      topo_results <- data.frame(Scenario = character(), Stat_Value = numeric(), 
                                 r_hat = numeric(), stringsAsFactors = FALSE)
      all_p_values <- numeric(replications * 2)
      ground_truth <- c(rep(0, replications), rep(1, replications))
      
      r_hat_h0_vec <- numeric(replications)
      r_hat_h1_vec <- numeric(replications)
      
      # Null Calibration
      t_bootstrap <- numeric(bootstrap_reps)
      for (b in 1:bootstrap_reps) {
        x_ref <- switch(topo,
                        "Normal"      = rnorm(n_curr, mean = 0, sd = 1),
                        "Exponential" = rexp(n_curr, rate = 1),
                        "Cauchy"      = rcauchy(n_curr, location = 0, scale = 1),
                        "Beta"        = rbeta(n_curr, shape1 = 2, shape2 = 5),
                        "Gamma"       = rgamma(n_curr, shape = 2, scale = 2),
                        "Weibull"     = rweibull(n_curr, shape = 1.5, scale = 1),
                        "Student_t"   = rt(n_curr, df = 3),
                        "Uniform"     = runif(n_curr, min = 0, max = 10))
        boot_res <- compute_t_statistic(x_ref)
        t_bootstrap[b] <- boot_res$stat
      }
      null_vector_sorted <- sort(t_bootstrap)
      t_threshold <- null_vector_sorted[ceiling(0.95 * bootstrap_reps)]
      
      # H0 Evaluation
      for (i in 1:replications) {
        x_null <- switch(topo,
                         "Normal"      = rnorm(n_curr, mean = 0, sd = 1),
                         "Exponential" = rexp(n_curr, rate = 1),
                         "Cauchy"      = rcauchy(n_curr, location = 0, scale = 1),
                         "Beta"        = rbeta(n_curr, shape1 = 2, shape2 = 5),
                         "Gamma"       = rgamma(n_curr, shape = 2, scale = 2),
                         "Weibull"     = rweibull(n_curr, shape = 1.5, scale = 1),
                         "Student_t"   = rt(n_curr, df = 3),
                         "Uniform"     = runif(n_curr, min = 0, max = 10))
        
        res_null <- compute_t_statistic(x_null)
        r_hat_h0_vec[i] <- res_null$r_hat
        topo_results <- rbind(topo_results, data.frame(Scenario = "Null_H0", 
                                                       Stat_Value = res_null$stat, 
                                                       r_hat = res_null$r_hat))
        all_p_values[i] <- mean(null_vector_sorted >= res_null$stat)
      }
      
      # H1 Evaluation
      for (i in 1:replications) {
        x_base <- switch(topo,
                         "Normal"      = rnorm(n_curr - r_curr, mean = 0, sd = 1),
                         "Exponential" = rexp(n_curr - r_curr, rate = 1),
                         "Cauchy"      = rcauchy(n_curr - r_curr, location = 0, scale = 1),
                         "Beta"        = rbeta(n_curr - r_curr, shape1 = 2, shape2 = 5),
                         "Gamma"       = rgamma(n_curr - r_curr, shape = 2, scale = 2),
                         "Weibull"     = rweibull(n_curr - r_curr, shape = 1.5, scale = 1),
                         "Student_t"   = rt(n_curr - r_curr, df = 3),
                         "Uniform"     = runif(n_curr - r_curr, min = 0, max = 10))
        
        x_outliers <- switch(topo,
                             "Normal"      = rnorm(r_curr, mean = 4, sd = 1),
                             "Exponential" = rexp(r_curr, rate = 0.08),
                             "Cauchy"      = rcauchy(r_curr, location = 12, scale = 1),
                             "Beta"        = rbeta(r_curr, shape1 = 18, shape2 = 2),
                             "Gamma"       = rgamma(r_curr, shape = 12, scale = 2),
                             "Weibull"     = rweibull(r_curr, shape = 1.5, scale = 8),
                             "Student_t"   = rt(r_curr, df = 3) + 8,
                             "Uniform"     = runif(r_curr, min = 15, max = 18))
        
        x_contaminated <- c(x_base, x_outliers)
        res_h1 <- compute_t_statistic(x_contaminated)
        r_hat_h1_vec[i] <- res_h1$r_hat
        topo_results <- rbind(topo_results, data.frame(Scenario = "Contaminated_H1", 
                                                       Stat_Value = res_h1$stat, 
                                                       r_hat = res_h1$r_hat))
        all_p_values[replications + i] <- mean(null_vector_sorted >= res_h1$stat)
      }
      
      topo_results$Decision <- ifelse(topo_results$Stat_Value > t_threshold, 1, 0)
      fpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Null_H0"])
      tpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Contaminated_H1"])
      
      roc_curve <- roc(ground_truth, 1 - all_p_values, quiet = TRUE)
      auc_val   <- as.numeric(auc(roc_curve))
      
      mean_r_h0 <- mean(r_hat_h0_vec)
      mean_r_h1 <- mean(r_hat_h1_vec)
      r_hat_acc <- mean(r_hat_h1_vec == r_curr)
      
      summary_results <- rbind(summary_results, data.frame(
        Topology        = topo,
        n_Sample        = n_curr,
        r_True          = r_curr,
        FPR             = fpr_val,
        TPR             = tpr_val,
        AUC             = auc_val,
        Mean_r_hat_H0   = mean_r_h0,
        Mean_r_hat_H1   = mean_r_h1,
        r_hat_Accuracy  = r_hat_acc
      ))
      
      # Aggregate plot data
      topo_results$n_Sample <- factor(sprintf("n = %d", n_curr))
      topo_results$r_True   <- factor(sprintf("r0 = %d", r_curr), levels = c("r0 = 1", "r0 = 3", "r0 = 5", "r0 = 10"))
      topo_density_data     <- rbind(topo_density_data, topo_results)
      
      roc_df <- data.frame(
        FPR = 1 - roc_curve$specificities,
        TPR = roc_curve$sensitivities,
        n_Sample = factor(sprintf("n = %d", n_curr)),
        r_True   = factor(sprintf("r0 = %d", r_curr), levels = c("r0 = 1", "r0 = 3", "r0 = 5", "r0 = 10"))
      )
      topo_roc_data <- rbind(topo_roc_data, roc_df)
      
      cat("Done.\n")
    }
  }
  
  # --- 4. EXPORT SUPERIMPOSED PLOTS FOR TOPOLOGY (2 PLOTS PER TOPOLOGY = 16 TOTAL) ---
  
  # 1. Superimposed Separability Density Profile across all r_0 and n settings
  p_sep <- ggplot(topo_density_data, aes(x = Stat_Value, color = r_True, linetype = Scenario)) +
    geom_density(linewidth = 0.8) +
    facet_wrap(~ n_Sample, scales = "free_y") +
    labs(title = sprintf("Separability Density Profile — %s Topology", topo),
         subtitle = expression("Comparing Null (" * H[0] * ") vs. Contaminated (" * H[1] * ") across Contamination Levels (" * r[0] * ")"),
         x = expression("Test Statistic Value " * T(X * ";" * hat(r)[n])),
         y = "Density",
         color = "Contamination (r0)",
         linetype = "Hypothesis") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  ggsave(file.path(main_plot_dir, sprintf("separability_combined_%s.png", tolower(topo))),
         plot = p_sep, width = 8.5, height = 5)
  
  # 2. Superimposed ROC Curves across all r_0 and n settings
  p_roc <- ggplot(topo_roc_data, aes(x = FPR, y = TPR, color = r_True)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
    facet_wrap(~ n_Sample) +
    labs(title = sprintf("Superimposed ROC Curves — %s Topology", topo),
         subtitle = expression("ROC performance across sample sizes (" * n * ") and contamination levels (" * r[0] * ")"),
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)",
         color = "Contamination (r0)") +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(main_plot_dir, sprintf("roc_combined_%s.png", tolower(topo))),
         plot = p_roc, width = 8.5, height = 5)
}

# --- 5. EXPORT SUMMARY MATRIX ---
write.csv(summary_results, "expanded_simulation_results.csv", row.names = FALSE)
