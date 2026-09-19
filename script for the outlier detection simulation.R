# =========================================================================
# COMPREHENSIVE EXTENDED MULTI-TOPOLOGY DENSITY OUTLIER SIMULATION STUDY
# Updated: Adaptive Sequential Density-Gap Estimator (\hat{r}_n) Integration
# =========================================================================

library(ggplot2)
library(stats)
library(pROC)
library(gridExtra)

# --- 1. CORE MATHEMATICAL ATTRIBUTES & DEFAULTS ---
set.seed(2026) # Reproducible execution state
n_sample <- 100
r_true   <- 3    # True contamination parameter r_0 under H1
replications <- 200
bootstrap_reps <- 250

# Extended to 8 total distributions
topologies <- c("Normal", "Exponential", "Cauchy", "Beta", 
                "Gamma", "Weibull", "Student_t", "Uniform")

# Tracking dataframe for final comparative matrix (including r_hat accuracy)
summary_results <- data.frame(
  Topology        = character(),
  FPR             = numeric(),
  TPR             = numeric(),
  AUC             = numeric(),
  Mean_r_hat_H0   = numeric(),
  Mean_r_hat_H1   = numeric(),
  r_hat_Accuracy  = numeric(), # Proportion of times \hat{r}_n == r_0 under H1
  stringsAsFactors = FALSE
)

# --- 2. ALGORITHMIC ARCHITECTURE IMPLEMENTATION ---

#' Calculate Robust Bandwidth (Silverman's rule with M-estimation Scale)
get_robust_bandwidth <- function(x) {
  n <- length(x)
  sigma_robust <- median(abs(x - median(x))) * 1.4826
  if (sigma_robust == 0) sigma_robust <- sd(x)
  h <- 0.9 * sigma_robust * (n^(-1/5))
  return(max(h, 1e-5))
}

#' Estimate r dynamically via Sequential Log-Density Gap Thresholding (Def 2.1)
estimate_r_hat <- function(densities, n, h, r_max = NULL) {
  if (is.null(r_max)) {
    r_max <- max(1, floor(5 * log(n)))
  }
  
  sorted_densities <- sort(densities, decreasing = FALSE)
  sorted_densities[sorted_densities <= 0] <- .Machine$double.eps
  
  # Log-density sequence
  log_dens <- log(sorted_densities)
  
  # Adaptive threshold: tau_n = C_r * sqrt(ln(n) / (n * h_n))
  C_r <- 1.0
  tau_n <- C_r * sqrt(log(n) / (n * h))
  
  # Sequential log-density gap rule: min { k : ln f(x_{k+1}) - ln f(x_k) < tau_n }
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

#' Compute Test Statistic T(X; \hat{r}_n) using Adaptive r_hat Estimation
compute_t_statistic <- function(x) {
  n <- length(x)
  h <- get_robust_bandwidth(x)
  
  densities <- numeric(n)
  for (i in 1:n) {
    densities[i] <- mean(dnorm((x[i] - x) / h)) / h
  }
  
  # Estimate r_hat dynamically per draw
  r_hat <- estimate_r_hat(densities, n, h)
  
  sorted_densities <- sort(densities, decreasing = FALSE)
  isolated_points <- sorted_densities[1:r_hat]
  isolated_points[isolated_points <= 0] <- .Machine$double.eps
  
  stat_value <- -mean(log(isolated_points))
  
  return(list(stat = stat_value, r_hat = r_hat))
}

# --- 3. EXECUTION PROCESSOR THROUGH ALL ENVIRONMENTS ---

for (topo in topologies) {
  cat(sprintf("\nProcessing Simulation Infrastructure for Topology: [%s]\n", topo))
  
  topo_results <- data.frame(Scenario = character(), Stat_Value = numeric(), r_hat = numeric(), stringsAsFactors = FALSE)
  all_p_values <- numeric(replications * 2)
  ground_truth <- c(rep(0, replications), rep(1, replications))
  
  # Track \hat{r}_n diagnostics
  r_hat_h0_vec <- numeric(replications)
  r_hat_h1_vec <- numeric(replications)
  
  # A. Phase 1: Reference Null Space Distribution Calibration (Bootstrap Loop)
  cat("  -> Calibrating structural null distribution space via dynamic bootstrap...\n")
  t_bootstrap <- numeric(bootstrap_reps)
  for (b in 1:bootstrap_reps) {
    x_ref <- switch(topo,
                    "Normal"      = rnorm(n_sample, mean = 0, sd = 1),
                    "Exponential" = rexp(n_sample, rate = 1),
                    "Cauchy"      = rcauchy(n_sample, location = 0, scale = 1),
                    "Beta"        = rbeta(n_sample, shape1 = 2, shape2 = 5),
                    "Gamma"       = rgamma(n_sample, shape = 2, scale = 2),
                    "Weibull"     = rweibull(n_sample, shape = 1.5, scale = 1),
                    "Student_t"   = rt(n_sample, df = 3),
                    "Uniform"     = runif(n_sample, min = 0, max = 10))
    
    boot_res <- compute_t_statistic(x_ref)
    t_bootstrap[b] <- boot_res$stat
  }
  null_vector_sorted <- sort(t_bootstrap)
  t_threshold <- null_vector_sorted[ceiling(0.95 * bootstrap_reps)]
  
  # B. Phase 2: H0 Simulation Loop (Empirical Type I Error Analysis)
  for (i in 1:replications) {
    x_null <- switch(topo,
                     "Normal"      = rnorm(n_sample, mean = 0, sd = 1),
                     "Exponential" = rexp(n_sample, rate = 1),
                     "Cauchy"      = rcauchy(n_sample, location = 0, scale = 1),
                     "Beta"        = rbeta(n_sample, shape1 = 2, shape2 = 5),
                     "Gamma"       = rgamma(n_sample, shape = 2, scale = 2),
                     "Weibull"     = rweibull(n_sample, shape = 1.5, scale = 1),
                     "Student_t"   = rt(n_sample, df = 3),
                     "Uniform"     = runif(n_sample, min = 0, max = 10))
    
    res_null <- compute_t_statistic(x_null)
    t_calc <- res_null$stat
    r_hat_h0_vec[i] <- res_null$r_hat
    
    topo_results <- rbind(topo_results, data.frame(Scenario = "Null_H0", Stat_Value = t_calc, r_hat = res_null$r_hat))
    all_p_values[i] <- mean(null_vector_sorted >= t_calc)
  }
  
  # C. Phase 3: H1 Simulation Loop (Statistical Power Assessment)
  for (i in 1:replications) {
    # Generate baseline uncontaminated support space
    x_base <- switch(topo,
                     "Normal"      = rnorm(n_sample - r_true, mean = 0, sd = 1),
                     "Exponential" = rexp(n_sample - r_true, rate = 1),
                     "Cauchy"      = rcauchy(n_sample - r_true, location = 0, scale = 1),
                     "Beta"        = rbeta(n_sample - r_true, shape1 = 2, shape2 = 5),
                     "Gamma"       = rgamma(n_sample - r_true, shape = 2, scale = 2),
                     "Weibull"     = rweibull(n_sample - r_true, shape = 1.5, scale = 1),
                     "Student_t"   = rt(n_sample - r_true, df = 3),
                     "Uniform"     = runif(n_sample - r_true, min = 0, max = 10))
    
    # Introduce slippage outliers configured to localized domains
    x_outliers <- switch(topo,
                         "Normal"      = rnorm(r_true, mean = 4, sd = 1),
                         "Exponential" = rexp(r_true, rate = 0.08),
                         "Cauchy"      = rcauchy(r_true, location = 12, scale = 1),
                         "Beta"        = rbeta(r_true, shape1 = 18, shape2 = 2),
                         "Gamma"       = rgamma(r_true, shape = 12, scale = 2),
                         "Weibull"     = rweibull(r_true, shape = 1.5, scale = 8),
                         "Student_t"   = rt(r_true, df = 3) + 8,
                         "Uniform"     = runif(r_true, min = 15, max = 18))
    
    x_contaminated <- c(x_base, x_outliers)
    res_h1 <- compute_t_statistic(x_contaminated)
    t_calc <- res_h1$stat
    r_hat_h1_vec[i] <- res_h1$r_hat
    
    topo_results <- rbind(topo_results, data.frame(Scenario = "Contaminated_H1", Stat_Value = t_calc, r_hat = res_h1$r_hat))
    all_p_values[replications + i] <- mean(null_vector_sorted >= t_calc)
  }
  
  # --- 4. QUANTITATIVE METRICS EXTRACTION ---
  topo_results$Decision <- ifelse(topo_results$Stat_Value > t_threshold, 1, 0)
  fpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Null_H0"])
  tpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Contaminated_H1"])
  
  roc_curve <- roc(ground_truth, 1 - all_p_values, quiet = TRUE)
  auc_val   <- as.numeric(auc(roc_curve))
  
  # Estimator Convergence Diagnostics (Lemma 2.2 verification)
  mean_r_h0 <- mean(r_hat_h0_vec)
  mean_r_h1 <- mean(r_hat_h1_vec)
  r_hat_acc <- mean(r_hat_h1_vec == r_true)
  
  summary_results <- rbind(summary_results, data.frame(
    Topology       = topo,
    FPR            = fpr_val,
    TPR            = tpr_val,
    AUC            = auc_val,
    Mean_r_hat_H0  = mean_r_h0,
    Mean_r_hat_H1  = mean_r_h1,
    r_hat_Accuracy = r_hat_acc
  ))
  
  # --- 5. GRAPHICAL EXPORT COMPILING PER TOPOLOGY ---
  p1 <- ggplot(topo_results, aes(x = Stat_Value, fill = Scenario)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = t_threshold, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("label", x = min(topo_results$Stat_Value) + (max(topo_results$Stat_Value)-min(topo_results$Stat_Value))*0.30, 
             y = 0.4, 
             label = sprintf("FPR: %.4f | TPR: %.4f\nAvg r_hat (H1): %.2f (Acc: %.1f%%)", 
                             fpr_val, tpr_val, mean_r_h1, r_hat_acc * 100), 
             fill = "white", fontface = "bold", size = 3.2) +
    labs(title = sprintf("Separability Profile (%s Topology)", topo),
         subtitle = expression("Evaluated under Adaptive Estimator " * hat(r)[n] * " Thresholding"),
         x = expression("Calculated Test Statistic Value " * T(X * ";" * hat(r)[n])), 
         y = "Density") +
    scale_fill_manual(values = c("Null_H0" = "#3498db", "Contaminated_H1" = "#e74c3c")) +
    theme_minimal()
  
  ggsave(sprintf("separability_%s.png", tolower(topo)), plot = p1, width = 6.5, height = 4.5)
  
  roc_data <- data.frame(Specificity = roc_curve$specificities, Sensitivity = roc_curve$sensitivities)
  p2 <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "#2c3e50", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
    annotate("label", x = 0.75, y = 0.25, 
             label = sprintf("AUC-ROC = %.4f\nPower = %.2f%%", auc_val, tpr_val * 100), 
             fill = "#fcf8e3", size = 4, fontface = "bold") +
    labs(title = sprintf("ROC Curve (%s Topology)", topo),
         subtitle = expression("Framework performance under adaptive " * hat(r)[n] * " selection"),
         x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    theme_minimal()
  
  ggsave(sprintf("roc_%s.png", tolower(topo)), plot = p2, width = 6.5, height = 4.5)
}

# --- 6. GLOBAL EVALUATION MATRIX SUMMARY OUTPUT ---
cat("\n=========================================================================================\n")
cat("          FINAL SIMULATION METRICS & ESTIMATOR CONVERGENCE ACROSS ALL TOPOLOGIES          \n")
cat("=========================================================================================\n")
print(summary_results, row.names = FALSE)
cat("=========================================================================================\n")
