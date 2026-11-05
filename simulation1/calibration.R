# =============================================================================
# calibration_L.R
# Check tcem_helpers.R for the functions
# Calibrate the EWMA-R control limit multiplier L so the in-control ARL hits
# =============================================================================

library(parallel) # to speed up the simulations with mclapply()

# when using parallel::mclapply(), each worker needs its own random stream. The
# L'Ecuyer-CMRG gives each forked mclapply worker an independent, non-
# overlapping random stream. Set the kind BEFORE seeding.

RNGkind("L'Ecuyer-CMRG")
set.seed(123)


n_sim  <- 1e3
max_rl <- 1e4 # to make sure the forever-running simulations 

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

const      <- build_tcem_constants() # precompute T^2 values
X          <- const$X
XtX_inv_Xt <- const$XtX_inv_Xt
Sigma_inv  <- const$Sigma_inv
U          <- const$U

n_cores <- max(1, detectCores() - 1) # leave one core free for other tasks (10 cores - 1)

arl_ewma_r <- function(sigma_eta_s, L, a0_s = a0, a1_s = a1, sigma_eps_s = sigma_eps) {
  rls <- unlist(mclapply(seq_len(n_sim), function(i) {
    run_length_ewma_r(a0_s, a1_s, sigma_eps_s, sigma_eta_s, L)
  }, mc.cores = n_cores))
  mean(rls)
}

arl_t2 <- function(sigma_eta_s, UCL_T2, a0_s = a0, a1_s = a1, sigma_eps_s = sigma_eps) {
  rls <- unlist(mclapply(seq_len(n_sim), function(i) {
    run_length_t2(a0_s, a1_s, sigma_eps_s, sigma_eta_s, UCL_T2)
  }, mc.cores = n_cores))
  mean(rls)
}

# coarse-to-fine search for L (target in-control ARL) 
# Two phases of bisection:
#   ROUGH  - bracket L quickly using an ARL tolerance (cheap, few iters).
#   FINE   - keep bisecting the surviving bracket until its width < 1e-4,
#            so the returned L is precise to 4 decimal places.
find_L <- function(eta, target_arl = 200, lower = 2, upper = 10, arl_tol = 3, L_tol = 1e-4, max_iter = 40) {

  mid     <- (lower + upper) / 2
  arl_mid <- arl_ewma_r(sigma_eta_s = eta, L = mid)
  phase   <- "rough"

  for (iter in seq_len(max_iter)) {
    mid     <- (lower + upper) / 2
    arl_mid <- arl_ewma_r(sigma_eta_s = eta, L = mid)

    #cat(sprintf("  [%-5s] iter %2d | L = %.4f | ARL = %8.1f ",phase, iter, mid, arl_mid))

    # Once the ARL is close enough, switch to refining L to 4 decimals.
    if (phase == "rough" && abs(arl_mid - target_arl) <= arl_tol) {phase <- "fine"}

    # Done only when L itself is pinned to 4 decimals.
    if (phase == "fine" && (upper - lower) < L_tol) break

    if (arl_mid < target_arl) lower <- mid else upper <- mid
  }

  list(L = round(mid, 4), ARL = arl_mid)
}

# run calibration across the sigma_eta grid 
L_calibrated   <- numeric(length(sigma_eta_levels))
L_achieved_ARL <- numeric(length(sigma_eta_levels))

for (i in seq_along(sigma_eta_levels)) {
  eta <- sigma_eta_levels[i]
  #cat("\nCalibrating L for eta =", eta, "\n")
  res <- find_L(eta)
  L_calibrated[i]   <- res$L
  L_achieved_ARL[i] <- res$ARL
}

L_table <- data.frame(
  sigma_eta    = sigma_eta_levels,
  L_ewma_r     = round(L_calibrated, 4),
  achieved_ARL = round(L_achieved_ARL, 1)
)

print(L_table)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(L_table, "calibrated_L.csv", row.names = FALSE)
write.csv(L_table, file.path("results/tables", "calibrated_L.csv"), row.names = FALSE)

# sigma_eta L_ewma_r achieved_ARL
# 1      0.00   3.1209        200.0
# 2      0.02   3.2429        199.9
# 3      0.05   3.8805        200.0
# 4      0.08   5.0470        200.0
# 5      0.10   6.0822        199.7

## dont change the above code


# same coarse-to-fine trick as find_L, but searching the T2 upper limit.
# bigger UCL_T2 -> longer ARL, so the bisection direction is the same.
find_UCL_T2 <- function(eta, target_arl = 200, lower = 2, upper = 40, arl_tol = 3, UCL_tol = 1e-4, max_iter = 40) {

  mid     <- (lower + upper) / 2
  arl_mid <- arl_t2(sigma_eta_s = eta, UCL_T2 = mid)
  phase   <- "rough"

  for (iter in seq_len(max_iter)) {
    mid     <- (lower + upper) / 2
    arl_mid <- arl_t2(sigma_eta_s = eta, UCL_T2 = mid)

    #cat(sprintf("  [%-5s] iter %2d | UCL_T2 = %.4f | ARL = %8.1f ", phase, iter, mid, arl_mid))

    # ARL close enough -> refine UCL_T2 to 4 decimals
    if (phase == "rough" && abs(arl_mid - target_arl) <= arl_tol) {phase <- "fine"}

    # stop only when UCL_T2 is pinned to 4 decimals
    if (phase == "fine" && (upper - lower) < UCL_tol) break

    if (arl_mid < target_arl) lower <- mid else upper <- mid
  }

  list(UCL_T2 = round(mid, 4), ARL = arl_mid)
}

UCL_T2_calibrated <- numeric(length(sigma_eta_levels))
UCL_T2_achieved_ARL <- numeric(length(sigma_eta_levels))

for (i in seq_along(sigma_eta_levels)) {
  eta <- sigma_eta_levels[i]
  #cat("\nCalibrating UCL_T2 for eta =", eta, "\n")
  res <- find_UCL_T2(eta)
  UCL_T2_calibrated[i] <- res$UCL_T2
  UCL_T2_achieved_ARL[i] <- res$ARL
}

# chi-square, so an implied alpha would no longer reflect a real false-alarm rate.
# control limits below were adjusted directly to hit ARL0 = 200.
UCL_table <- data.frame(
  sigma_eta = sigma_eta_levels,
  UCL_T2 = round(UCL_T2_calibrated, 4),
  achieved_ARL = round(UCL_T2_achieved_ARL, 1)
)

print(UCL_table)
write.csv(UCL_table, "calibrated_UCL_T2.csv", row.names = FALSE)
write.csv(UCL_table, file.path("results/tables", "calibrated_UCL_T2.csv"), row.names = FALSE)

# sigma_eta  UCL_T2 achieved_ARL
# 1      0.00 10.6195        199.9
# 2      0.02 11.0981        200.0
# 3      0.05 14.5088        199.6
# 4      0.08 21.2893        200.4
# 5      0.10 28.2535        200.1