# =========================================================================
# COMPREHENSIVE EXTENDED MULTI-TOPOLOGY DENSITY OUTLIER SIMULATION STUDY
# Target: Validation across eight continuous distributional topologies
# =========================================================================

library(ggplot2)
library(stats)
library(pROC)
library(gridExtra)

# --- 1. CORE MATHEMATICAL ATTRIBUTES & DEFAULTS ---
set.seed(2026) # Reproducible execution state
n_sample <- 100
r_suspected <- 3
replications <- 200
bootstrap_reps <- 250

# Extended to 8 total distributions
topologies <- c("Normal", "Exponential", "Cauchy", "Beta", 
                "Gamma", "Weibull", "Student_t", "Uniform")

# Tracking dataframe for final comparative matrix
summary_results <- data.frame(
  Topology = character(),
  FPR      = numeric(),
  TPR      = numeric(),
  AUC      = numeric(),
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

#' Compute the Generalized Non-Parametric Test Statistic T(X; r)
compute_t_statistic <- function(x, r) {
  n <- length(x)
  h <- get_robust_bandwidth(x)
  
  densities <- numeric(n)
  for (i in 1:n) {
    densities[i] <- mean(dnorm((x[i] - x) / h)) / h
  }
  
  sorted_densities <- sort(densities, decreasing = FALSE)
  isolated_points <- sorted_densities[1:r]
  isolated_points[isolated_points <= 0] <- .Machine$double.eps
  
  return(-mean(log(isolated_points)))
}

# --- 3. EXECUTION PROCESSOR THROUGH ALL ENVIRONMENTS ---

for (topo in topologies) {
  cat(sprintf("\nProcessing Simulation Infrastructure for Topology: [%s]\n", topo))
  
  topo_results <- data.frame(Scenario = character(), Stat_Value = numeric(), stringsAsFactors = FALSE)
  all_p_values <- numeric(replications * 2)
  ground_truth <- c(rep(0, replications), rep(1, replications))
  
  # A. Phase 1: Reference Null Space Distribution Calibration
  cat("  -> Calibrating structural null distribution space via bootstrap...\n")
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
    t_bootstrap[b] <- compute_t_statistic(x_ref, r_suspected)
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
    t_calc <- compute_t_statistic(x_null, r_suspected)
    
    topo_results <- rbind(topo_results, data.frame(Scenario = "Null_H0", Stat_Value = t_calc))
    all_p_values[i] <- mean(null_vector_sorted >= t_calc)
  }
  
  # C. Phase 3: H1 Simulation Loop (Statistical Power Assessment)
  for (i in 1:replications) {
    # Generate baseline uncontaminated support space
    x_base <- switch(topo,
                     "Normal"      = rnorm(n_sample - r_suspected, mean = 0, sd = 1),
                     "Exponential" = rexp(n_sample - r_suspected, rate = 1),
                     "Cauchy"      = rcauchy(n_sample - r_suspected, location = 0, scale = 1),
                     "Beta"        = rbeta(n_sample - r_suspected, shape1 = 2, shape2 = 5),
                     "Gamma"       = rgamma(n_sample - r_suspected, shape = 2, scale = 2),
                     "Weibull"     = rweibull(n_sample - r_suspected, shape = 1.5, scale = 1),
                     "Student_t"   = rt(n_sample - r_suspected, df = 3),
                     "Uniform"     = runif(n_sample - r_suspected, min = 0, max = 10))
    
    # Introduce slippage outliers configured to localized domains
    x_outliers <- switch(topo,
                         "Normal"      = rnorm(r_suspected, mean = 4, sd = 1),
                         "Exponential" = rexp(r_suspected, rate = 0.08),
                         "Cauchy"      = rcauchy(r_suspected, location = 12, scale = 1),
                         "Beta"        = rbeta(r_suspected, shape1 = 18, shape2 = 2),
                         "Gamma"       = rgamma(r_suspected, shape = 12, scale = 2),   # Right-tail contamination
                         "Weibull"     = rweibull(r_suspected, shape = 1.5, scale = 8),# Scale enlargement
                         "Student_t"   = rt(r_suspected, df = 3) + 8,                 # Location shift
                         "Uniform"     = runif(r_suspected, min = 15, max = 18))       # Out-of-bounds shift
    
    x_contaminated <- c(x_base, x_outliers)
    t_calc <- compute_t_statistic(x_contaminated, r_suspected)
    
    topo_results <- rbind(topo_results, data.frame(Scenario = "Contaminated_H1", Stat_Value = t_calc))
    all_p_values[replications + i] <- mean(null_vector_sorted >= t_calc)
  }
  
  # --- 4. QUANTITATIVE METRICS EXTRACTION ---
  topo_results$Decision <- ifelse(topo_results$Stat_Value > t_threshold, 1, 0)
  fpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Null_H0"])
  tpr_val <- mean(topo_results$Decision[topo_results$Scenario == "Contaminated_H1"])
  
  roc_curve <- roc(ground_truth, 1 - all_p_values, quiet = TRUE)
  auc_val   <- as.numeric(auc(roc_curve))
  
  summary_results <- rbind(summary_results, data.frame(
    Topology = topo,
    FPR      = fpr_val,
    TPR      = tpr_val,
    AUC      = auc_val
  ))
  
  # --- 5. GRAPHICAL EXPORT COMPILING PER TOPOLOGY ---
  p1 <- ggplot(topo_results, aes(x = Stat_Value, fill = Scenario)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = t_threshold, linetype = "dashed", color = "red", size = 1) +
    annotate("label", x = min(topo_results$Stat_Value) + (max(topo_results$Stat_Value)-min(topo_results$Stat_Value))*0.25, 
             y = 0.4, 
             label = sprintf("Empirical FPR: %.4f\nEmpirical TPR (Power): %.4f", fpr_val, tpr_val), 
             fill = "white", fontface = "bold", size = 3.5) +
    labs(title = sprintf("Separability Profile (%s Topology)", topo),
         subtitle = "Red dashed line signifies the calibrated validation boundary",
         x = "Calculated Test Statistic Value", y = "Density") +
    scale_fill_manual(values = c("Null_H0" = "#3498db", "Contaminated_H1" = "#e74c3c")) +
    theme_minimal()
  
  ggsave(sprintf("separability_%s.png", tolower(topo)), plot = p1, width = 6.5, height = 4.5)
  
  roc_data <- data.frame(Specificity = roc_curve$specificities, Sensitivity = roc_curve$sensitivities)
  p2 <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_line(color = "#2c3e50", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
    annotate("label", x = 0.75, y = 0.25, 
             label = sprintf("AUC-ROC = %.4f\nPower = %.2f%%", auc_val, tpr_val * 100), 
             fill = "#fcf8e3", size = 4, fontface = "bold") +
    labs(title = sprintf("ROC Curve (%s Topology)", topo),
         subtitle = "Framework performance under distributional limits",
         x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)") +
    theme_minimal()
  
  ggsave(sprintf("roc_%s.png", tolower(topo)), plot = p2, width = 6.5, height = 4.5)
}

# --- 6. GLOBAL EVALUATION MATRIX SUMMARY OUTPUT ---
cat("\n====================================================\n")
cat("   FINAL SIMULATION METRICS ACROSS ALL TOPOLOGIES     \n")
cat("====================================================\n")
print(summary_results, row.names = FALSE)
cat("====================================================\n")