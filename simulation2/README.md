# Simulation Study II: Nickel ICP-MS Calibration Application

This folder contains the code for **Simulation Study II**, which evaluates Phase II profile monitoring under the Two-Component Error Model (TCEM) using nickel calibration parameters from Inductively Coupled Plasma Mass Spectrometry (ICP-MS) (Wilson et al., 2004; Rocke & Lorenzato, 1995).

---

## 1. Directory Structure

```
simulation2/
├── README.md              # This file
├── est.R                  # Complete calibration of L and alpha (target ARL0 = 200)
├── main.R                 # Phase II simulation, out-of-control percentages and plots
└── results/
    ├── tables/            # results_summary.csv (out-of-control percentages)
    └── plots/             # 12 PNG + 12 PDF three-panel diagnostic stacks
```

---

## 2. Experimental Parameters & Baseline Design

* **TCEM model**:
  $$Y_{ij} = A_0 + A_1 X_i e^{\eta_{ij}} + \epsilon_{ij}$$
* **Parameters (Rocke & Lorenzato, 1995 / Wilson et al., 2004)**:
  * Intercept: $A_0 = 670.0$
  * Slope: $A_1 = 10.3$
  * Additive error SD: $\sigma_\epsilon = 43.3$
  * Multiplicative error SD: $\sigma_\eta = 0.226$
* **Design**: $X = (0, 10, 20, 50, 100)$ ppt, $n = 5$, one observation per level.
* **Phase I**: $m = 30$ in-control profiles (seed 123) estimate $\bar{A}_0$, $\bar{A}_1$, $\hat{\sigma}$, the mean of the parameter vectors and $S^{-1}$.
* **Calibrated limits (ARL$_0$ = 200, conditional on the Phase I sample above)**:
  * EWMA-R multiplier: $L = 8.244736$
  * Hotelling $T^2$ significance level: $\alpha = 0.002737542$, giving UCL$_{T^2}$ = 15.202
  * Verified ARL$_0$ with 10,000 replications: EWMA-R 201.1 (SE 2.0), $T^2$ 203.2 (SE 2.1)

---

## 3. Evaluated Shift Scenarios

Shifts start at Phase II profile 31 and are applied one at a time (no combined shifts), using the single-$\sigma_\epsilon$ convention (Kang & Albin, 2000):

1. **Slope shift ($\beta$)**: $A_1^* = A_1 + \beta \sigma_\epsilon$, with $\beta \in \{0, 0.10, 0.25\}$
2. **Intercept shift ($\lambda$)**: $A_0^* = A_0 + \lambda \sigma_\epsilon$, with $\lambda \in \{0, 0.80, 2.00\}$
3. **Additive error inflation**: $\sigma_\epsilon^* = \gamma \sigma_\epsilon$, with $\gamma \in \{1.0, 2.0, 3.0\}$
4. **Multiplicative error inflation**: $\sigma_\eta^* = \gamma_\eta \sigma_\eta$, with $\gamma_\eta \in \{1.0, 2.0, 3.0\}$

---

## 4. Metric

For each scenario, `main.R` reports the **percentage of out-of-control profiles among the 30 Phase II profiles**, averaged over $N = 500$ Monte Carlo replications. Columns:

* `avg_pct_EWMA`: EWMA component alone
* `avg_pct_R`: range component alone
* `avg_pct_EWMA_R`: EWMA-R union chart (EWMA or range signal)
* `avg_pct_T2`: Hotelling $T^2$ chart

In-control scenarios (shift = 0 or the baseline value) give false-alarm rates near 0.5% per profile for an ARL$_0$ of 200.

---

## 5. Reproducing the Results

Run both scripts from the `simulation2/` folder:

```bash
Rscript est.R     # calibration: writes results/tables/calibration_limits.csv
Rscript main.R    # simulation: writes results/tables/results_summary.csv and plots
```

`main.R` uses the constants hard-coded in its section 2c. If the calibration changes, update those values to match the `est.R` output.

| Output file | Description |
| :--- | :--- |
| `results/tables/calibration_limits.csv` | Calibrated L, α, UCL$_{T^2}$ and verified ARL$_0$ (from `est.R`) |
| `results/tables/results_summary.csv` | Out-of-control percentages for all 12 scenarios |
| `results/plots/beta_*.png / .pdf` | Three-panel diagnostic stacks (EWMA, range, $T^2$) for slope shifts |
| `results/plots/alpha_*.png / .pdf` | Same for intercept shifts |
| `results/plots/sigma_eps_*.png / .pdf` | Same for additive error inflation |
| `results/plots/sigma_eta_*.png / .pdf` | Same for multiplicative error inflation |

Runtime: `est.R` takes several minutes. `main.R` with 500 replications takes a few minutes. Required packages (installed automatically by `main.R` via `pacman`): tidyverse, patchwork, viridis.
