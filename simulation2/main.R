# ==============================================================================
# Simulation Study II: Nickel ICP-MS Case Study
# Profile Monitoring under Two-Component Error Model (TCEM)
#
# Generates Phase I reference data (Wilson et al., 2004 nickel calibration)
# Simulates Phase II shifts in slope (beta), intercept (alpha), additive noise
# (sigma_eps), and multiplicative noise (sigma_eta).
#
# Metric Reported:
#   Percentage of out-of-control profiles across a fixed 30-run Phase II horizon,
#   averaged over n_sims replications.
# ==============================================================================

# 1. INITIAL SETUP
set.seed(123)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load("tidyverse", "patchwork", "viridis")

# Output directory setup (contained within simulation2)
out_dir_plots  <- "results/plots"
out_dir_tables <- "results/tables"

dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)


# 2. STUDY DESIGN & PARAMETERS

## 2a. Ground Truth (Wilson et al., 2004 / Rocke & Lorenzato, 1995 Nickel calibration)
alpha_true     <- 670.0    # Intercept
beta_true      <- 10.3     # Slope
sigma_eps_true <- 43.3     # Additive standard deviation
sigma_eta_true <- 0.226    # Multiplicative standard deviation

## 2b. Replication level per concentration level
reps <- 1   # Baseline online profile monitoring: 1 observation per level (n = 5)

## 2c. Calibrated constants for ARL0 = 200, conditional on this Phase I sample (from est.R)
calib <- list(
  "1" = list(L = 8.244736, alpha = 0.002737542, d2 = 2.326, d3 = 0.864)
)
cc <- calib[[as.character(reps)]]

## 2d. Profile design points and simulation settings
levels_x    <- c(0, 10, 20, 50, 100)     # 5-level concentration design (ppt)
X_design    <- rep(levels_x, each = reps)
n_x         <- length(X_design)
X_bar       <- mean(X_design)
S_xx        <- sum((X_design - X_bar)^2)

p1_runs     <- 30                        # Phase I baseline sample size
p2_runs     <- 30                        # Phase II fixed evaluation window
run_range2  <- (p1_runs + 1):(p1_runs + p2_runs)
n_sims      <- 500                       # Monte Carlo replications
max_rl      <- 5000                      # Censoring cap for sequential run-length

theta       <- 0.2                       # EWMA smoothing parameter
p_dim       <- 2                         # Hotelling T2 dimension (intercept, slope)

## 2e. Shift design (Kang & Albin single-sigma convention)
# All parameter shifts are expressed as multiples of residual SD sigma_eps_true = 43.3
beta_mults      <- c(0, 0.10, 0.25)     # Null, moderate, large slope shift
alpha_mults     <- c(0, 0.80, 2.00)     # Null, moderate, large intercept shift
sigma_eps_mults <- c(1.0, 2.0, 3.0)     # Additive error inflation
sigma_eta_mults <- c(1.0, 2.0, 3.0)     # Multiplicative error inflation

## 2f. Two-Component Error Model (TCEM) profile generator
tcem_gen <- function(a0, a1, eta_sd, eps_sd, X) {
  a0 + a1 * X * exp(rnorm(length(X), 0, eta_sd)) + rnorm(length(X), 0, eps_sd)
}

## Fast OLS helper
ols_fast <- function(X, Y) {
  xb <- mean(X); yb <- mean(Y)
  b1 <- sum((X - xb) * (Y - yb)) / sum((X - xb)^2)
  c(a0 = yb - b1 * xb, a1 = b1)
}

## 2g. Shift values in physical measurement units
beta_shifts      <- beta_true      + beta_mults      * sigma_eps_true
alpha_shifts     <- alpha_true     + alpha_mults     * sigma_eps_true
sigma_eps_shifts <- sigma_eps_true * sigma_eps_mults
sigma_eta_shifts <- sigma_eta_true * sigma_eta_mults


# 3. PHASE I REFERENCE ESTIMATION & CONTROL LIMITS

## 3a. Generate Phase I historical profiles
sim_p1   <- expand.grid(X = X_design, run = 1:p1_runs)
sim_p1$Y <- tcem_gen(alpha_true, beta_true, sigma_eta_true, sigma_eps_true, sim_p1$X)

## 3b. Per-profile OLS parameter estimation
p1_ols <- sim_p1 %>%
  group_by(run) %>%
  summarise(
    a0  = coef(lm(Y ~ X))[1],
    a1  = coef(lm(Y ~ X))[2],
    mse = sum(resid(lm(Y ~ X))^2) / (n() - 2),
    .groups = "drop"
  )

a0_ref    <- mean(p1_ols$a0)
a1_ref    <- mean(p1_ols$a1)
sigma_ref <- sqrt(mean(p1_ols$mse))

## 3c. Phase I residual statistics and covariance matrix
p1_resid <- sim_p1 %>%
  mutate(e = Y - (a0_ref + a1_ref * X)) %>%
  group_by(run) %>%
  summarise(ebar = mean(e), R = max(e) - min(e), .groups = "drop") %>%
  mutate(ewma = accumulate(ebar, ~ theta * .y + (1 - theta) * .x))

coef_mat  <- as.matrix(p1_ols[, c("a0", "a1")])
z_bar     <- colMeans(coef_mat)
S_inv     <- solve(cov(coef_mat))

p1_ols$T2 <- apply(coef_mat, 1, function(z)
  as.numeric(t(z - z_bar) %*% S_inv %*% (z - z_bar)))
p1_resid$T2 <- p1_ols$T2

## 3d. Control Limits calibrated to target ARL0 = 200
UCL_EWMA <-  cc$L * sigma_ref * sqrt(theta / ((2 - theta) * n_x))
LCL_EWMA <- -UCL_EWMA

UCL_R <- sigma_ref * (cc$d2 + cc$L * cc$d3)
LCL_R <- max(0, sigma_ref * (cc$d2 - cc$L * cc$d3))
CL_R  <- sigma_ref * cc$d2

ucl_T2 <- (p_dim * (p1_runs - 1) / (p1_runs - p_dim)) *
  qf(1 - cc$alpha, df1 = p_dim, df2 = p1_runs - p_dim)


# 4. PHASE II SIMULATION WITH DUAL METRICS

simulate_scenario <- function(shift_vals, make_params, param_name) {

  # ----------------------------------------------------------------------------
  # Metric 1: Fixed 30-Run Phase II Horizon (% of Out-of-Control Points)
  # ----------------------------------------------------------------------------
  raw <- map_dfr(shift_vals, function(v) {
    p <- make_params(v)
    map_dfr(seq_len(n_sims), function(sim_id) {
      df <- expand.grid(X = X_design, run = run_range2)
      df$Y      <- tcem_gen(p$a0, p$a1, p$eta_sd, p$eps_sd, df$X)
      df$sim_id <- sim_id
      df$shift  <- v
      df
    })
  })

  p2_resid <- raw %>%
    mutate(e = Y - (a0_ref + a1_ref * X)) %>%
    group_by(shift, sim_id, run) %>%
    summarise(ebar = mean(e), R = max(e) - min(e), .groups = "drop")

  p2_ols <- raw %>%
    group_by(shift, sim_id, run) %>%
    summarise(a0 = coef(lm(Y ~ X))[1], a1 = coef(lm(Y ~ X))[2], .groups = "drop") %>%
    mutate(T2 = apply(cbind(a0, a1), 1, function(z)
      as.numeric(t(z - z_bar) %*% S_inv %*% (z - z_bar))))

  # Prepend Phase I for smooth EWMA warm-up
  grid_sims <- distinct(p2_resid, shift, sim_id)
  p1_expand <- tidyr::crossing(grid_sims, p1_resid %>% select(run, ebar, R, T2))

  charts <- bind_rows(
    p1_expand,
    p2_resid %>% left_join(p2_ols %>% select(shift, sim_id, run, T2),
                           by = c("shift", "sim_id", "run"))
  ) %>%
    arrange(shift, sim_id, run) %>%
    group_by(shift, sim_id) %>%
    mutate(ewma = accumulate(ebar, ~ theta * .y + (1 - theta) * .x)) %>%
    ungroup()

  per_run_pct <- charts %>%
    filter(run > p1_runs) %>%
    mutate(
      ooc_ewma   = ewma > UCL_EWMA | ewma < LCL_EWMA,
      ooc_r      = R    > UCL_R    | R    < LCL_R,
      ooc_t2     = T2   > ucl_T2,
      ooc_ewma_r = ooc_ewma | ooc_r
    ) %>%
    group_by(shift, sim_id) %>%
    summarise(
      pct_ewma   = 100 * mean(ooc_ewma),
      pct_r      = 100 * mean(ooc_r),
      pct_t2     = 100 * mean(ooc_t2),
      pct_ewma_r = 100 * mean(ooc_ewma_r),
      .groups    = "drop"
    )

  pct_table <- per_run_pct %>%
    group_by(shift) %>%
    summarise(
      avg_pct_EWMA   = round(mean(pct_ewma), 2),
      avg_pct_R      = round(mean(pct_r), 2),
      avg_pct_EWMA_R = round(mean(pct_ewma_r), 2),
      avg_pct_T2     = round(mean(pct_t2), 2),
      .groups        = "drop"
    )

  # Scenario table: out-of-control percentages only
  combined_table <- pct_table %>%
    mutate(param = param_name, .before = 1)

  list(charts = charts, table = combined_table)
}


# 5. EXECUTE ALL FOUR SHIFT SCENARIOS

cat("--> Simulating slope shift (beta)...\n")
beta_res <- simulate_scenario(
  beta_shifts,
  function(b) list(a0 = alpha_true, a1 = b,
                   eta_sd = sigma_eta_true, eps_sd = sigma_eps_true),
  "beta"
)

cat("--> Simulating intercept shift (alpha)...\n")
alpha_res <- simulate_scenario(
  alpha_shifts,
  function(a) list(a0 = a, a1 = beta_true,
                   eta_sd = sigma_eta_true, eps_sd = sigma_eps_true),
  "alpha"
)

cat("--> Simulating additive error shift (sigma_eps)...\n")
sigma_eps_res <- simulate_scenario(
  sigma_eps_shifts,
  function(s) list(a0 = alpha_true, a1 = beta_true,
                   eta_sd = sigma_eta_true, eps_sd = s),
  "sigma_eps"
)

cat("--> Simulating multiplicative error shift (sigma_eta)...\n")
sigma_eta_res <- simulate_scenario(
  sigma_eta_shifts,
  function(s) list(a0 = alpha_true, a1 = beta_true,
                   eta_sd = s, eps_sd = sigma_eps_true),
  "sigma_eta"
)


# 6. ASSEMBLE SUMMARY TABLES & EXPORT

results_table <- bind_rows(
  beta_res$table,
  alpha_res$table,
  sigma_eps_res$table,
  sigma_eta_res$table
) %>%
  mutate(reps = reps, .before = 1) %>%
  arrange(param, shift)

cat("\n========================= COMPLETE RESULTS SUMMARY =========================\n")
print(as.data.frame(results_table))

# Save summary table
write.csv(results_table, file.path(out_dir_tables, "results_summary.csv"), row.names = FALSE)


# 7. GENERATE PUBLICATION-QUALITY PLOT PANELS

make_title <- function(chart_expr, sym, value, avg_pct) {
  bquote(.(chart_expr) * ": " ~ .(sym) == .(value) * ",  " ~
           .(round(avg_pct, 1)) * "% out of control")
}

plot_chart <- function(charts, target_shift, chart_var, ylab,
                       center, limits_vec, title_expr) {
  d <- charts %>% filter(shift == target_shift)

  ggplot(d, aes(x = run, y = .data[[chart_var]], group = sim_id)) +
    geom_line(alpha = 0.10, color = "#0072B2") +
    (if (!is.null(center)) geom_hline(yintercept = center, linetype = 2) else NULL) +
    geom_hline(yintercept = limits_vec, linetype = 2, color = "firebrick") +
    geom_vline(xintercept = p1_runs + 0.5, linetype = 3) +
    labs(title = title_expr, x = "Run", y = ylab) +
    theme_classic()
}

# shift_vals: physical parameter values; mults: shift multipliers used in titles
plot_shift_set <- function(scenario, shift_vals, mults, sym) {
  tab <- scenario$table

  stacks <- map2(shift_vals, mults, function(v, m) {
    pct <- function(col) tab[[col]][tab$shift == v]
    p_ewma <- plot_chart(scenario$charts, v, "ewma", "EWMA", 0,
                         c(UCL_EWMA, LCL_EWMA),
                         make_title(quote(EWMA), sym, m, pct("avg_pct_EWMA")))
    p_r <- plot_chart(scenario$charts, v, "R", "Range",
                      CL_R, c(UCL_R, LCL_R),
                      make_title(quote(Range), sym, m, pct("avg_pct_R")))
    p_t2 <- plot_chart(scenario$charts, v, "T2", expression(T^2), NULL,
                       ucl_T2,
                       make_title(quote(T^2), sym, m, pct("avg_pct_T2")))
    p_ewma / p_r / p_t2
  })
  names(stacks) <- as.character(round(mults, 2))
  stacks
}

plots <- list(
  beta      = plot_shift_set(beta_res,      beta_shifts,      beta_mults,      quote(beta)),
  alpha     = plot_shift_set(alpha_res,     alpha_shifts,     alpha_mults,     quote(lambda)),
  sigma_eps = plot_shift_set(sigma_eps_res, sigma_eps_shifts, sigma_eps_mults, quote(gamma)),
  sigma_eta = plot_shift_set(sigma_eta_res, sigma_eta_shifts, sigma_eta_mults, quote(gamma[eta]))
)

cat("\n--> Saving plot figures to results/plots/ (PNG and Cairo PDF)...\n")
for (nm in names(plots)) {
  for (stack_name in names(plots[[nm]])) {
    p <- plots[[nm]][[stack_name]]

    # Clean standardized filename: e.g. beta_14.63.png, sigma_eps_86.6.png
    fname_clean <- file.path(out_dir_plots, paste0(nm, "_", stack_name))
    ggsave(paste0(fname_clean, ".png"), plot = p, width = 7.62, height = 7.83, units = "in", dpi = 300, bg = "white")
    ggsave(paste0(fname_clean, ".pdf"), plot = p, width = 7.62, height = 7.83, units = "in", device = cairo_pdf)
  }
}

cat("\nSimulation Study II execution complete. All tables and plots generated.\n")
