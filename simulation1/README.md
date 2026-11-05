# Simulation Study I: Controlled Assessment of Multiplicative Measurement Error

This folder contains the complete simulation and analysis code for **Simulation Study I**, evaluating the Phase II performance of linear profile monitoring schemes under the Two-Component Error Model (TCEM) using the mass-flow-controller (MFC) experimental setting of Kang & Albin (2000).

All generated output tables and publication figures are stored directly inside **`simulation1/results/`**.

---

## 1. Directory Structure

```
simulation1/
├── README.md              # Detailed study and reproducibility guide
├── tcem_helpers.R         # Core mathematical helpers (data generator, run lengths, constants)
├── calibration.R          # Control-limit calibration (bisection search for target ARL0 = 200)
├── shifts_sim.R           # Parallel Monte Carlo shift evaluations (1D & 2D parameter shifts)
├── summary.R              # Summary table generation and publication-grade plot rendering
└── results/
    ├── tables/            # 30 CSV tables (calibrated limits, 1D shifts, 2D grids, validations)
    └── plots/             # 21 high-res PNG plots (ARL curves, multi-panel heatmaps, line charts)
```

---

## 2. Model & Experimental Parameters

* **Profile Equation**:
  $$Y_{ij} = A_0 + A_1 X_i e^{\eta_{ij}} + \epsilon_{ij}$$
* **Design Points**: $X = (2, 4, 6, 8)$ ($n = 4$)
* **In-Control Baseline**: $A_0 = 3, \; A_1 = 2, \; \sigma_\epsilon = 1$
* **Error Components**:
  * Multiplicative: $\eta_{ij} \sim N(0, \sigma_\eta^2), \quad \sigma_\eta \in \{0, 0.02, 0.05, 0.08, 0.10\}$
  * Additive: $\epsilon_{ij} \sim N(0, \sigma_\epsilon^2), \quad \sigma_\epsilon = 1$
* **Target In-Control Performance**: $\text{ARL}_0 = 200$
* **EWMA Smoothing Parameter**: $\theta = 0.20$

---

## 3. Control Charts Monitored

1. **EWMA-R Residual Chart** (Kang & Albin):
   * Monitors mean residual shifts via EWMA:
     $$z_j = \theta \bar{e}_j + (1 - \theta) z_{j-1}, \quad \text{UCL/LCL} = \pm L \sigma_\epsilon \sqrt{\frac{\theta}{(2 - \theta) n}}$$
   * Monitors residual dispersion via Range ($R$):
     $$\text{UCL}_R = \sigma_\epsilon (d_2 + L d_3), \quad \text{LCL}_R = \max(0, \sigma_\epsilon (d_2 - L d_3))$$
     where $d_2 = 2.059$ and $d_3 = 0.880$ for $n = 4$.
2. **Hotelling's $T^2$ Chart**:
   * Jointly monitors bivariate parameter estimates $Z_j = (\hat{A}_{0j}, \hat{A}_{1j})^T$:
     $$T_j^2 = (Z_j - U)^T \Sigma^{-1} (Z_j - U)$$

---

## 4. Plot & Table Reproducibility Guide

All figures and tables in `results/` are completely reproducible from the scripts:

| Script | Outputs Produced in `results/` |
| :--- | :--- |
| `calibration.R` | • `tables/calibrated_L.csv`<br>• `tables/calibrated_UCL_T2.csv` |
| `shifts_sim.R` | • `tables/shift_arl_intercept.csv`<br>• `tables/shift_arl_slope.csv`<br>• `tables/shift_arl_sigma_eps.csv`<br>• `tables/shift_arl_combined.csv`<br>• `tables/shift2d_arl_eta*.csv`<br>• `plots/shift_plot_intercept.png`<br>• `plots/shift_plot_slope.png`<br>• `plots/shift_plot_sigma_eps.png`<br>• `plots/shift_plot_combined.png`<br>• `plots/shift2d_plot_eta*.png` |
| `summary.R` | • `tables/summary_calibrated_limits.csv`<br>• `tables/summary_incontrol_eta_fixed_limits.csv`<br>• `tables/summary_simulation_check.csv`<br>• `tables/summary_slope_intercept_arl_long.csv`<br>• `tables/summary_slope_intercept_arl_eta*.csv`<br>• `plots/summary_incontrol_eta_fixed_limits.png`<br>• `plots/summary_calibrated_limits.png`<br>• `plots/summary_slope_intercept_arl_grid.png`<br>• `plots/summary_slope_intercept_arl_lines.png`<br>• `plots/summary_slope_intercept_arl_eta*_1x2.png` |

---

## 5. Execution Instructions

From the terminal within the `simulation1/` folder:

```bash
# Step 1 (Optional): Run control-limit calibration (ARL0 = 200)
Rscript calibration.R

# Step 2: Run high-throughput shift simulations (multi-core mclapply)
Rscript shifts_sim.R

# Step 3: Generate summary tables and publication figures
Rscript summary.R
```
