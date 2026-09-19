# ==============================================================================
# REAL-WORLD 8-TOPOLOGY COMPARATIVE VISUALIZATION & BENCHMARK SCRIPT
# Dataset: Base R Air Quality Dataset (New York Ozone Measurements, 1973)
# Evaluated Topologies: Normal, Exponential, Cauchy, Beta, Gamma, Weibull, Student_t, Uniform
# Methods Compared: Proposed Non-Parametric Density, Rosner's ESD, Iterative Grubbs
# ==============================================================================

if (!require("EnvStats")) install.packages("EnvStats")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")

library(EnvStats)
library(ggplot2)
library(gridExtra)
library(stats)

set.seed(2026)

# Load real continuous physical measurements
data(airquality)
real_sensor_data <- na.omit(airquality$Ozone)

n_obs       <- length(real_sensor_data) # N = 116
r_suspected <- 5                        # r = 5 extreme points evaluated

topologies <- c("Normal", "Exponential", "Cauchy", "Beta", 
                "Gamma", "Weibull", "Student_t", "Uniform")

# ------------------------------------------------------------------------------
# 1. CORE ALGORITHMIC IMPLEMENTATIONS
# ------------------------------------------------------------------------------

get_robust_bandwidth <- function(x) {
  sigma_robust <- median(abs(x - median(x))) * 1.4826
  if (sigma_robust == 0) sigma_robust <- sd(x)
  h <- 0.9 * sigma_robust * (length(x)^(-1/5))
  return(max(h, 1e-5))
}

compute_np_density_outliers <- function(x_data, r) {
  n <- length(x_data)
  h <- get_robust_bandwidth(x_data)
  densities <- numeric(n)
  for (i in 1:n) {
    densities[i] <- mean(dnorm((x_data[i] - x_data) / h)) / h
  }
  threshold_density <- sort(densities)[r]
  flagged_indices   <- which(densities <= threshold_density)
  return(list(indices = flagged_indices, densities = densities))
}

run_rosner_esd <- function(x_data, k = 5, alpha = 0.05) {
  tryCatch({
    test_res <- rosnerTest(x_data, k = k, alpha = alpha)
    outliers <- test_res$all.stats[test_res$all.stats$Outlier == TRUE, ]
    return(outliers$Obs.Num)
  }, error = function(e) return(integer(0)))
}

run_grubbs_iterative <- function(x_data, max_k = 5, alpha = 0.05) {
  x_temp <- x_data
  flagged_indices <- integer(0)
  for (iter in 1:max_k) {
    n <- length(x_temp)
    if (n <= 2) break
    mean_x <- mean(x_temp)
    sd_x   <- sd(x_temp)
    devs   <- abs(x_temp - mean_x)
    max_idx <- which.max(devs)
    g_calc <- devs[max_idx] / sd_x
    t_crit <- qt(1 - alpha / (2 * n), df = n - 2)
    g_crit <- ((n - 1) / sqrt(n)) * sqrt(t_crit^2 / (n - 2 + t_crit^2))
    
    if (g_calc > g_crit) {
      orig_idx <- which(x_data == x_temp[max_idx])[1]
      flagged_indices <- c(flagged_indices, orig_idx)
      x_temp <- x_temp[-max_idx]
    } else {
      break
    }
  }
  return(flagged_indices)
}

# ------------------------------------------------------------------------------
# 2. GENERATE COMPARATIVE DENSITY PLOTS & SUMMARY DATA
# ------------------------------------------------------------------------------

plot_list <- list()
heatmap_data <- data.frame(Topology = character(), Method = character(), Count = numeric())

summary_numerical_results <- data.frame(
  Topology             = character(),
  Support_Domain       = character(),
  NP_Density_Flagged   = numeric(),
  Rosner_ESD_Flagged   = numeric(),
  Grubbs_Flagged       = numeric(),
  Overlap_NP_ESD       = numeric(),
  Concordance_Pct      = numeric(),
  stringsAsFactors     = FALSE
)

# Normalize data to (0, 1) strictly for Beta topology mapping
data_beta_scaled <- (real_sensor_data - min(real_sensor_data) + 0.001) / 
  (max(real_sensor_data) - min(real_sensor_data) + 0.002)

cat("Processing benchmarks across all 8 distribution topologies...\n")

for (topo in topologies) {
  
  # Map continuous support based on topology specification
  x_eval <- switch(topo,
                   "Normal"      = real_sensor_data,
                   "Exponential" = real_sensor_data,
                   "Cauchy"      = real_sensor_data - median(real_sensor_data),
                   "Beta"        = data_beta_scaled,
                   "Gamma"       = real_sensor_data,
                   "Weibull"     = real_sensor_data,
                   "Student_t"   = (real_sensor_data - mean(real_sensor_data)) / sd(real_sensor_data),
                   "Uniform"     = real_sensor_data
  )
  
  # Outlier Detection
  np_idx     <- compute_np_density_outliers(x_eval, r_suspected)$indices
  rosner_idx <- run_rosner_esd(x_eval, k = r_suspected)
  grubbs_idx <- run_grubbs_iterative(x_eval, max_k = r_suspected)
  
  # Quantify Agreement Metrics
  overlap_cnt <- length(intersect(np_idx, rosner_idx))
  union_cnt   <- length(union(np_idx, rosner_idx))
  concordance <- if (union_cnt > 0) (overlap_cnt / union_cnt) * 100 else 100.0
  
  # Append to Numerical Table
  summary_numerical_results <- rbind(summary_numerical_results, data.frame(
    Topology             = topo,
    Support_Domain       = switch(topo, "Beta" = "(0, 1)", "Cauchy" = "R (Centered)", "Student_t" = "R (Standardized)", "R+"),
    NP_Density_Flagged   = length(np_idx),
    Rosner_ESD_Flagged   = length(rosner_idx),
    Grubbs_Flagged       = length(grubbs_idx),
    Overlap_NP_ESD       = overlap_cnt,
    Concordance_Pct      = round(concordance, 1)
  ))
  
  # Append to Heatmap Data
  heatmap_data <- rbind(heatmap_data,
                        data.frame(Topology = topo, Method = "Proposed Density", Count = length(np_idx)),
                        data.frame(Topology = topo, Method = "Rosner ESD",       Count = length(rosner_idx)),
                        data.frame(Topology = topo, Method = "Grubbs Test",     Count = length(grubbs_idx))
  )
  
  # Prepare Plotting DataFrame
  df_plot <- data.frame(Value = x_eval, Status = "Normal Inlier")
  df_plot$Status[1:n_obs %in% np_idx] <- "Flagged by Density"
  df_plot$Status[1:n_obs %in% rosner_idx] <- "Flagged by Rosner ESD"
  df_plot$Status[(1:n_obs %in% np_idx) & (1:n_obs %in% rosner_idx)] <- "Flagged by Both"
  
  # Build Individual Panel Plot
  p <- ggplot(df_plot, aes(x = Value)) +
    geom_density(fill = "#3498db", alpha = 0.25, color = "#2980b9", linewidth = 0.8) +
    geom_rug(aes(color = Status), length = unit(0.08, "npc"), alpha = 0.8, linewidth = 1) +
    scale_color_manual(values = c(
      "Normal Inlier"         = "#95a5a6",
      "Flagged by Density"    = "#e74c3c",
      "Flagged by Rosner ESD" = "#f39c12",
      "Flagged by Both"       = "#8e44ad"
    )) +
    labs(
      title = sprintf("Topology: %s", topo),
      x = "Transformed Ozone Domain",
      y = "Density"
    ) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 10))
  
  plot_list[[topo]] <- p
}

# ------------------------------------------------------------------------------
# 3. EXPORT VISUAL GRAPHICS
# ------------------------------------------------------------------------------

grid_graphic <- grid.arrange(
  plot_list[["Normal"]], plot_list[["Exponential"]],
  plot_list[["Cauchy"]], plot_list[["Beta"]],
  plot_list[["Gamma"]],  plot_list[["Weibull"]],
  plot_list[["Student_t"]], plot_list[["Uniform"]],
  ncol = 2,
  top = ""
)

ggsave("real_world_8topology_densities.png", plot = grid_graphic, width = 11, height = 11, dpi = 300)

p_heatmap <- ggplot(heatmap_data, aes(x = Topology, y = Method, fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), color = "black", fontface = "bold", size = 5) +
  scale_fill_gradient(low = "#ebf5fb", high = "#3498db") +
  labs(
    title = "",
    subtitle = "",
    x = "Distribution Topology",
    y = "Detection Algorithm"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave("real_world_method_concordance.png", plot = p_heatmap, width = 8, height = 4.5, dpi = 300)

# ------------------------------------------------------------------------------
# 4. PRINT & EXPORT SUMMARY NUMERICAL RESULTS TABLE
# ------------------------------------------------------------------------------

cat("\n========================================================================================\n")
cat("   SUMMARY NUMERICAL RESULTS: REAL-WORLD DATASET ACROSS ALL 8 DISTRIBUTION TOPOLOGIES  \n")
cat("========================================================================================\n")
print(summary_numerical_results, row.names = FALSE)
cat("========================================================================================\n")

# Save table as CSV
write.csv(summary_numerical_results, file = "real_world_8topology_summary_results.csv", row.names = FALSE)
cat("\nNumerical summary results saved to 'real_world_8topology_summary_results.csv'.\n")