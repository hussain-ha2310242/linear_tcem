# Linear Profile Monitoring under Two-Component Measurement Error (TCEM)

Research replication repository for:
**"Effect of Two-Component Error Model on Linear Profile Monitoring"**

This repository contains reproducible R code, simulation pipelines, and publication-grade results investigating the performance of Phase II linear profile monitoring schemes under the **Two-Component Error Model (TCEM)** (Rocke & Lorenzato, 1995).

---

## 1. Research Overview

Standard Statistical Process Control (SPC) techniques for profile monitoring typically assume that measurement error is purely additive with constant variance ($\epsilon \sim N(0, \sigma_\epsilon^2)$). However, in analytical chemistry, spectrometry, and precision manufacturing, measurement errors frequently scale with the level of the response.

We analyze profile monitoring under the **Two-Component Error Model (TCEM)**:

$$Y_{ij} = A_0 + A_1 X_i e^{\eta_{ij}} + \epsilon_{ij}$$

where:
* $A_0, A_1$: Linear profile intercept and slope.
* $X_i$: Independent design/concentration points.
* $\eta_{ij} \sim N(0, \sigma_\eta^2)$: Multiplicative measurement error.
* $\epsilon_{ij} \sim N(0, \sigma_\epsilon^2)$: Additive measurement error.

We evaluate and compare two Phase II control charting techniques:
1. **EWMA-R Residual Chart** (Kang & Albin, 2000): An EWMA chart on average profile residuals to detect mean shifts, paired with a Range ($R$) chart on residuals to detect dispersion changes.
2. **Hotelling's $T^2$ Chart**: Joint bivariate monitoring of least-squares parameter estimates $(\hat{A}_0, \hat{A}_1)^T$.

---

## 2. Repository Layout

Each simulation study is fully self-contained with its own code, calibration routines, and dedicated `results/` folder for tables and figures:

```
linear_tcem/
├── README.md                      # Repository overview, theoretical model, and replication guide
│
├── simulation1/                   # Simulation Study I: Controlled MFC Setting (Kang & Albin, 2000)
│   ├── README.md                  # Comprehensive guide for Simulation Study I
│   ├── tcem_helpers.R             # Core profile generator, run length functions, constants
│   ├── calibration.R              # Binary search calibrating L and UCL_T2 to ARL0 = 200
│   ├── shifts_sim.R               # Parallelized Monte Carlo simulation across 1D & 2D shifts
│   ├── summary.R                  # Generates summary tables and publication-grade plots
│   └── results/                   # Self-contained outputs for Study I
│       ├── tables/                # 30 CSV tables (calibrated limits, 1D shifts, 2D grids)
│       └── plots/                 # 21 high-res PNG plots (ARL curves, multi-panel heatmaps)
│
└── simulation2/                   # Simulation Study II: Empirical Nickel ICP-MS Case Study
    ├── README.md                  # Comprehensive guide for Simulation Study II
    ├── est.R                      # Standalone calibration functions (target ARL0 = 200)
    ├── main.R                     # Phase I estimation, 4 shift types, dual-metric reporting
    └── results/                   # Self-contained outputs for Study II
        ├── tables/                # results_summary.csv (dual metrics across all shifts)
        └── plots/                 # 24 plots (12 PNG + 12 vector Cairo PDF diagnostic stacks)
```

---

## 3. Simulation Studies Summary

| Feature | Simulation Study I | Simulation Study II |
| :--- | :--- | :--- |
| **Application** | Mass Flow Controller (MFC) chamber pressure | Nickel ICP-MS chemical calibration |
| **Reference** | Kang & Albin (2000) | Wilson et al. (2004) / Rocke & Lorenzato (1995) |
| **Design Points ($X$)** | $X = (2, 4, 6, 8)$ ($n = 4$) | $X = (0, 10, 20, 50, 100)$ ppt ($n = 5$) |
| **In-Control Parameters** | $A_0 = 3, \; A_1 = 2, \; \sigma_\epsilon = 1$ | $A_0 = 670.0, \; A_1 = 10.3, \; \sigma_\epsilon = 43.3, \; \sigma_\eta = 0.226$ |
| **Multiplicative Levels** | $\sigma_\eta \in \{0, 0.02, 0.05, 0.08, 0.10\}$ | $\sigma_\eta = 0.226$ (fixed baseline) |
| **Target $\text{ARL}_0$** | $200$ (recalibrated per $\sigma_\eta$ level) | $200$ ($L = 6.534, \; \alpha = 0.003697$) |
| **Phase I Sample** | Theoretical parameters | $m = 30$ in-control profiles |
| **Evaluated Shifts** | Intercept ($\lambda$), slope ($\beta$), $\sigma_\epsilon$ inflation, combined $\delta$, joint 2D grid | Slope ($\beta$), intercept ($\alpha$), $\sigma_\epsilon$ inflation, $\sigma_\eta$ inflation |
| **Evaluation Metrics** | Sequential run length (Average Run Length $\text{ARL}$) | Dual metrics: (1) % OOC in 30 runs, and (2) Sequential $\text{ARL} \pm \text{SE}$ |

---

## 4. Replication Instructions

### Prerequisites

R (version $\ge 4.2$) with required packages:

```r
install.packages(c("tidyverse", "patchwork", "viridis", "parallel"))
```

### Reproducing Simulation Study I

Navigate to `simulation1/`:

```bash
cd simulation1

# 1. (Optional) Re-run control limit calibration
Rscript calibration.R

# 2. Run parallel shift simulations
Rscript shifts_sim.R

# 3. Generate summary tables and plots
Rscript summary.R
```

*All tables and plots are written to `simulation1/results/tables/` and `simulation1/results/plots/`.*

### Reproducing Simulation Study II

Navigate to `simulation2/`:

```bash
cd simulation2

# Run Phase I baseline estimation, Phase II shift simulation, and plot export
Rscript main.R
```

*All tables and figures are written to `simulation2/results/tables/` and `simulation2/results/plots/`.*
