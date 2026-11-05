# shifts_sim.R
# Out-of-control study: how fast do the EWMA-R and T2 charts detect a shift,
# for each measurement-error level sigma_eta? We use the limits already found
# in calibration.R (read from the CSVs) so every chart starts at ARL0 = 200.
# Functions live in tcem_helpers.R. This file just drives them in parallel.

library(parallel)
library(ggplot2)

RNGkind("L'Ecuyer-CMRG")   # independent random stream per mclapply worker
set.seed(123)

n_sim  <- 1e5
max_rl <- 1e5
verbose <- TRUE


# Kang settings
x     <- c(2, 4, 6, 8)
n     <- length(x)
x_bar <- mean(x)
S_xx  <- sum((x - x_bar)^2)

a0        <- 3
a1        <- 2
sigma_eps <- 1
theta     <- 0.2
d2        <- 2.059
d3        <- 0.880

sigma_eta_levels <- c(0, 0.02, 0.05, 0.08, 0.10)

source("tcem_helpers.R")

const      <- build_tcem_constants() # this function
X          <- const$X
XtX_inv_Xt <- const$XtX_inv_Xt
Sigma_inv  <- const$Sigma_inv
U          <- const$U

n_cores <- max(1, detectCores() - 1)

log_msg <- function(...) {
  if (isTRUE(verbose)) {
    cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
    flush.console()
  }
}

run_parallel_batches <- function(indices, worker, label) {
  batches <- split(indices, ceiling(seq_along(indices) / n_cores))
  out <- vector("list", length(indices))

  for (b in seq_along(batches)) {
    batch <- batches[[b]]
    log_msg(
      sprintf(
        "%s batch %d/%d: cells %d-%d",
        label, b, length(batches), min(batch), max(batch)
      )
    )

    batch_out <- mclapply(batch, worker, mc.cores = min(n_cores, length(batch)))
    out[batch] <- batch_out

    log_msg(
      sprintf(
        "%s batch %d/%d done",
        label, b, length(batches)
      )
    )
  }

  out
}

# output folders (created if missing)
tables_dir <- "results/tables"
plots_dir  <- "results/plots"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

# calibrated limits per eta (from calibration.R or results)
find_calib_file <- function(fname) {
  candidates <- c(file.path(tables_dir, fname), fname)
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) stop("Missing calibration file: ", fname, ". Run calibration.R first.")
  found[1]
}
L_tab   <- read.csv(find_calib_file("calibrated_L.csv"))
UCL_tab <- read.csv(find_calib_file("calibrated_UCL_T2.csv"))
L_by_eta   <- setNames(L_tab$L_ewma_r, L_tab$sigma_eta)
UCL_by_eta <- setNames(UCL_tab$UCL_T2, UCL_tab$sigma_eta)

# full plan: every (eta, shift) cell we need to simulate
plan <- merge(data.frame(sigma_eta = sigma_eta_levels), shift_grid())

# run all cells in parallel, one worker per cell
log_msg("Running", nrow(plan), "one-dimensional shift cells on", n_cores, "cores")
cells <- run_parallel_batches(seq_len(nrow(plan)), function(k) {
  row <- plan[k, ]
  eta <- row$sigma_eta
  shift_cell_arl(
    eta        = eta,
    shift_type = row$shift_type,
    shift      = row$shift,
    L          = L_by_eta[[as.character(eta)]],
    UCL_T2     = UCL_by_eta[[as.character(eta)]],
    n_sim      = n_sim
  )
}, "1D shifts")

results <- do.call(rbind, cells)

# wide ARL table for one shift type and one chart (rows = eta, cols = shift)
shift_table <- function(shift_type, chart) {
  d <- results[results$shift_type == shift_type & results$chart == chart, ]
  m <- tapply(d$ARL, list(d$sigma_eta, d$shift), mean)
  tab <- data.frame(chart = chart, sigma_eta = as.numeric(rownames(m)), round(m, 3))
  colnames(tab)[-(1:2)] <- colnames(m)
  tab
}




shift_types <- c("intercept", "slope", "sigma_eps", "combined")
for (st in shift_types) {
  cat("\n== Shift:", st, "==\n")
  tab <- rbind(shift_table(st, "EWMA-R"), shift_table(st, "T2"))
  print(tab, row.names = FALSE)
  write.csv(tab, file.path(tables_dir, paste0("shift_arl_", st, ".csv")), row.names = FALSE)
}

# plots: log(ARL) vs shift, one panel per chart
plots <- list(
  plot_shift(results, "intercept", "Intercept shift", "lambda"),
  plot_shift(results, "slope",     "Slope shift", "beta"),
  plot_shift(results, "sigma_eps", "Additive SD multiplier", "gamma"),
  plot_shift(results, "combined",  "Combined slope-intercept shift", "|delta|", use_abs = TRUE)
)
names(plots) <- shift_types

for (st in shift_types) {
  ggsave(file.path(plots_dir, paste0("shift_plot_", st, ".png")), plots[[st]], width = 10, height = 4, dpi = 300)
}
cat("\nDone. Tables and plots written.\n")



# Joint intercept x slope shift: rows = intercept, cols = slope, one per eta.
# This is a full 2D grid so it is much bigger than the studies above.

lambda_vals <- seq(0, 2.0,  by = 0.2)    # intercept shifts (table rows)
beta_vals   <- seq(0, 0.25, by = 0.025)  # slope shifts (table columns)

# every (eta, lambda, beta) cell
plan2d <- merge(
  data.frame(sigma_eta = sigma_eta_levels),
  merge(data.frame(lambda = lambda_vals), data.frame(beta = beta_vals))
)

log_msg("Running", nrow(plan2d), "joint intercept x slope cells on", n_cores, "cores")
cells2d <- run_parallel_batches(seq_len(nrow(plan2d)), function(k) {
  row <- plan2d[k, ]
  eta <- row$sigma_eta
  shift2d_cell_arl(
    eta    = eta,
    lambda = row$lambda,
    beta   = row$beta,
    L      = L_by_eta[[as.character(eta)]],
    UCL_T2 = UCL_by_eta[[as.character(eta)]],
    n_sim  = n_sim
  )
}, "2D shifts")

results2d <- do.call(rbind, cells2d)

# wide ARL table for one eta and one chart (rows = intercept, cols = slope)
shift2d_table <- function(eta, chart) {
  d <- results2d[results2d$sigma_eta == eta & results2d$chart == chart, ]
  m <- tapply(d$ARL, list(d$lambda, d$beta), mean)
  tab <- data.frame(chart = chart, intercept = as.numeric(rownames(m)), round(m, 3))
  colnames(tab)[-(1:2)] <- colnames(m)   # slope values as column names
  tab
}

# one table per (eta, chart) and one heatmap per eta
for (eta in sigma_eta_levels) {
  for (chart in c("EWMA-R", "T2")) {
    cat("\n== Joint intercept x slope | eta =", eta, "| chart =", chart, "==\n")
    tab <- shift2d_table(eta, chart)
    print(tab, row.names = FALSE)
    write.csv(tab, file.path(tables_dir, paste0("shift2d_arl_eta", eta, "_", chart, ".csv")), row.names = FALSE)
  }
  ggsave(file.path(plots_dir, paste0("shift2d_plot_eta", eta, ".png")), plot_shift2d(results2d, eta),
         width = 11, height = 4.5, dpi = 300)
}
cat("\nDone. Joint intercept x slope tables and plots written.\n")
