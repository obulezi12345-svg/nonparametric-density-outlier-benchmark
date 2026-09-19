# nonparametric-density-outlier-benchmark
A non-parametric kernel density estimation framework for outlier detection in univariate continuous data across 8 distribution topologies, benchmarked against Rosner's ESD and Grubbs' tests.
# Non-Parametric Density Estimation Outlier Detection & 8-Topology Benchmark

An implementation and comparative evaluation framework for generalized non-parametric kernel density estimation (KDE) outlier detection in continuous univariate data. 

This repository provides simulation and real-world evaluation workflows comparing the proposed non-parametric test statistic $\mathcal{T}(X; r)$ against standard parametric techniques (Rosner's Generalized ESD test and Iterative Grubbs' test) across **eight distinct distributional topologies**:
- **Normal** ($\mathbb{R}$)
- **Exponential** ($\mathbb{R}^+$)
- **Cauchy** (Heavy-tailed)
- **Beta** (Bounded support $(0, 1)$)
- **Gamma** (Skewed $\mathbb{R}^+$)
- **Weibull** (Asymmetric scale/shape)
- **Student's $t$** (Fat-tailed)
- **Uniform** (Bounded step support)

---

## 📌 Features

1. **Non-Parametric Density Outlier Statistic $\mathcal{T}(X; r)$**:
   - Uses robust bandwidth estimation via Median Absolute Deviation (MAD) scaled Silvermans' rule.
   - Evaluates logarithmic density infima over $r$ suspected anomalous observations.
2. **Bootstrap Null Calibration**:
   - Non-parametric critical value determination under standard non-contaminated reference null distributions ($H_0$).
3. **Multi-Topology Benchmark Suite**:
   - Monte Carlo simulation evaluating False Positive Rate (FPR), True Positive Rate / Power (TPR), and Area Under the ROC Curve (AUC).
4. **Real-World Empirical Validation**:
   - Application to the historic New York Ozone Air Quality dataset (`airquality$Ozone`).
   - Comparative method concordance heatmaps and density overlap visualizations.

---

## 📁 Repository Structure

```text
.
├── airquality data application.R        # Real-world benchmark on NYC Ozone dataset
├── script for the outlier detection simulation.R  # Monte Carlo multi-topology simulation study
├── real_world_8topology_densities.png   # Output grid plot (Generated upon run)
├── real_world_method_concordance.png    # Output heatmap plot (Generated upon run)
├── real_world_8topology_summary_results.csv # Empirical result output
└── README.md                            # Project documentation
