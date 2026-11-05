# =============================================================================
# tcem_helpers.R
# =============================================================================

# 1. Derived design constants
#    Computed once from the model parameters; reused by every replicate.
build_tcem_constants <- function() {
  X          <- cbind(1, x)
  XtX_inv_Xt <- solve(t(X) %*% X) %*% t(X)

  sigma_0_sq <- sigma_eps^2 * (1 / n + x_bar^2 / S_xx)
  sigma_1_sq <- sigma_eps^2 / S_xx
  sigma_01   <- -sigma_eps^2 * x_bar / S_xx
  Sigma_inv  <- solve(matrix(c(sigma_0_sq, sigma_01,
                               sigma_01,  sigma_1_sq), nrow = 2))

  list(
    X          = X,
    XtX_inv_Xt = XtX_inv_Xt,
    Sigma_inv  = Sigma_inv,
    U          = c(a0, a1)
  )
}



# 2. Data-generating process - One simulated linear profile with two-component (multiplicative + additive) measurement error.

simulate_profile <- function(a0_s, a1_s, sigma_eps_s, sigma_eta_s) {
  a0_s + a1_s * x * exp(rnorm(n, 0, sigma_eta_s)) + rnorm(n, 0, sigma_eps_s)
}



# 3. Per-replicate run lengths
#    Each returns the run length (capped at max_rl) for a single realisation
#    of the monitored process. These are the units of work handed to mclapply.


# EWMA on the mean residual + R chart on the residual range.
run_length_ewma_r <- function(a0_s, a1_s, sigma_eps_s, sigma_eta_s, L) {
  UCL_EWMA <-  L * sigma_eps * sqrt(theta / ((2 - theta) * n))
  LCL_EWMA <- -UCL_EWMA
  UCL_R    <- sigma_eps * (d2 + L * d3)
  LCL_R    <- max(0, sigma_eps * (d2 - L * d3))

  z  <- 0
  rl <- 0
  repeat {
    rl <- rl + 1
    y     <- simulate_profile(a0_s, a1_s, sigma_eps_s, sigma_eta_s)
    resid <- y - (a0 + a1 * x)

    z <- theta * mean(resid) + (1 - theta) * z
    R <- max(resid) - min(resid)

    if (z > UCL_EWMA || z < LCL_EWMA || R > UCL_R || R < LCL_R || rl >= max_rl) break
  }
  rl
}

# Hotelling T^2 on the fitted regression coefficients.
run_length_t2 <- function(a0_s, a1_s, sigma_eps_s, sigma_eta_s, UCL_T2) {
  rl <- 0
  repeat {
    rl <- rl + 1
    y         <- simulate_profile(a0_s, a1_s, sigma_eps_s, sigma_eta_s)
    Z         <- as.vector(XtX_inv_Xt %*% y)
    diff_beta <- Z - U
    t2        <- as.numeric(t(diff_beta) %*% Sigma_inv %*% diff_beta)

    if (t2 > UCL_T2 || rl >= max_rl) break
  }
  rl
}



# 4. Reporting

make_table <- function(res_matrix, shift_values, sigma_eta_levels) {
  tab <- data.frame(sigma_eta = sigma_eta_levels, round(res_matrix, 1))
  colnames(tab)[-1] <- shift_values
  tab
}



# 5. Shift study (Kang tables)
#    We move the true process away from in-control and see how fast each chart
#    reacts. The control limits and centre line stay at their in-control values.

# the four shift types and the shift values we scan for each
shift_grid <- function() {
  rbind(
    data.frame(shift_type = "intercept", shift_symbol = "lambda", shift = seq(0, 2.0, by = 0.2)),
    data.frame(shift_type = "slope",     shift_symbol = "beta",   shift = seq(0, 0.25, by = 0.025)),
    data.frame(shift_type = "sigma_eps", shift_symbol = "gamma",  shift = seq(1.0, 3.0, by = 0.2)),
    data.frame(shift_type = "combined",  shift_symbol = "delta",  shift = c(0, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0))
  )
}

# turn one shift into the true (shifted) process parameters
apply_shift <- function(shift_type, shift) {
  pars <- list(a0 = a0, a1 = a1, sigma_eps = sigma_eps)

  if (shift_type == "intercept") pars$a0 <- a0 + shift * sigma_eps
  if (shift_type == "slope")     pars$a1 <- a1 + shift * sigma_eps
  if (shift_type == "sigma_eps") pars$sigma_eps <- sigma_eps * shift
  if (shift_type == "combined") {
    pars$a1 <- a1 + shift * sigma_eps
    pars$a0 <- a0 - shift * sigma_eps * x_bar
  }
  pars
}

# ARL of both charts for one grid cell (one eta, one shift value).
# n_sim run lengths are averaged; this is the unit of work for mclapply.
shift_cell_arl <- function(eta, shift_type, shift, L, UCL_T2, n_sim) {
  p <- apply_shift(shift_type, shift)

  rl_ewma <- sapply(seq_len(n_sim), function(i)
    run_length_ewma_r(p$a0, p$a1, p$sigma_eps, eta, L))
  rl_t2 <- sapply(seq_len(n_sim), function(i)
    run_length_t2(p$a0, p$a1, p$sigma_eps, eta, UCL_T2))

  data.frame(
    sigma_eta  = eta,
    shift_type = shift_type,
    shift      = shift,
    chart      = c("EWMA-R", "T2"),
    ARL        = c(mean(rl_ewma), mean(rl_t2))
  )
}

# ARL-vs-shift plot, one panel per chart, one line per sigma_eta.
# Dashed line marks the in-control target (ARL0 = 200).
plot_shift <- function(results, shift_type, title, xlab, use_abs = FALSE, log_y = TRUE) {
  library(ggplot2)

  d <- results[results$shift_type == shift_type, ]
  d$sigma_eta <- factor(d$sigma_eta)
  d$x <- if (use_abs) abs(d$shift) else d$shift
  d$y <- if (log_y) log(d$ARL) else d$ARL

  ggplot(d, aes(x = x, y = y, linetype = sigma_eta, shape = sigma_eta, group = sigma_eta)) +
    geom_hline(yintercept = if (log_y) log(200) else 200, linetype = "dashed") +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ chart, nrow = 1) +
    labs(title = title, x = xlab,
         y = if (log_y) "log(ARL)" else "ARL",
         linetype = expression(sigma[eta]), shape = expression(sigma[eta])) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}



# 6. Joint intercept x slope shift
#    Both location parameters move at once, so the result is a 2D surface:
#    rows = intercept shift (lambda), columns = slope shift (beta), one per eta.

# ARL of both charts for one (lambda, beta) cell at a given eta
shift2d_cell_arl <- function(eta, lambda, beta, L, UCL_T2, n_sim) {
  a0_s <- a0 + lambda * sigma_eps
  a1_s <- a1 + beta   * sigma_eps

  rl_ewma <- sapply(seq_len(n_sim), function(i)
    run_length_ewma_r(a0_s, a1_s, sigma_eps, eta, L))
  rl_t2 <- sapply(seq_len(n_sim), function(i)
    run_length_t2(a0_s, a1_s, sigma_eps, eta, UCL_T2))

  data.frame(
    sigma_eta = eta,
    lambda    = lambda,
    beta      = beta,
    chart     = c("EWMA-R", "T2"),
    ARL       = c(mean(rl_ewma), mean(rl_t2))
  )
}

# heatmap of the intercept x slope ARL surface, one panel per chart
plot_shift2d <- function(results, eta) {
  library(ggplot2)

  d <- results[results$sigma_eta == eta, ]

  ggplot(d, aes(x = factor(lambda), y = factor(beta), fill = log(ARL))) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(ARL)), size = 2.5) +
    facet_wrap(~ chart, nrow = 1) +
    scale_fill_gradient(low = "steelblue", high = "white") +
    labs(title = bquote("Intercept x slope shift, " ~ sigma[eta] ~ "=" ~ .(eta)),
         x = "intercept shift (lambda)", y = "slope shift (beta)", fill = "log(ARL)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}