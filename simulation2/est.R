# ==============================================================================
# est.R: complete calibration of the EWMA-R and T^2 charts (target ARL0 = 200)
# Conditional on the Phase I sample generated in main.R (seed 123).
# Run from the simulation2/ folder:  Rscript est.R
# ==============================================================================

# ---- Settings -----------------------------------------------------------------
truth  <- list(a0 = 670.0, a1 = 10.3, eta_sd = 0.226, eps_sd = 43.3)
X      <- c(0, 10, 20, 50, 100)   # one observation per level, n = 5
n      <- length(X)
m      <- 30                      # Phase I profiles
p      <- 2                       # T^2 dimension
theta  <- 0.2                     # EWMA smoothing
target <- 200                     # target in-control ARL0
d2     <- 2.326                   # R-chart constants for n = 5
d3     <- 0.864
max_rl <- 10000                   # run-length cap
n_sims <- 50000                    # replications for calibration
n_ver  <- 10000                   # replications for verification

# ---- Phase I sample and reference (same draws as main.R) ----------------------
set.seed(123)
Xrep <- rep(X, times = m)
Y    <- truth$a0 + truth$a1 * Xrep * exp(rnorm(m * n, 0, truth$eta_sd)) +
        rnorm(m * n, 0, truth$eps_sd)
Ymat <- matrix(Y, m, n, byrow = TRUE)

Xc     <- X - mean(X)
Sxx    <- sum(Xc^2)
b1     <- as.vector(Ymat %*% Xc) / Sxx
b0     <- rowMeans(Ymat) - b1 * mean(X)
resid  <- Ymat - outer(b0, rep(1, n)) - outer(b1, X)

a0_ref <- mean(b0)
a1_ref <- mean(b1)
sig    <- sqrt(mean(rowSums(resid^2) / (n - 2)))
zbar   <- c(a0_ref, a1_ref)
Si     <- solve(cov(cbind(b0, b1)))

# ---- Phase II simulation (vectorised over replications) -----------------------
draw <- function(k) {
  Xm  <- matrix(X, k, n, byrow = TRUE)
  eta <- matrix(rnorm(k * n, 0, truth$eta_sd), k, n)
  eps <- matrix(rnorm(k * n, 0, truth$eps_sd), k, n)
  truth$a0 + truth$a1 * Xm * exp(eta) + eps
}

ewma_rl <- function(L, seed, reps = n_sims) {
  set.seed(seed)
  ucl  <- L * sig * sqrt(theta / ((2 - theta) * n))
  uclR <- sig * (d2 + L * d3)
  lclR <- max(0, sig * (d2 - L * d3))
  mu   <- a0_ref + a1_ref * X
  z    <- numeric(reps)
  rl   <- rep(max_rl, reps)
  act  <- seq_len(reps)
  for (s in seq_len(max_rl)) {
    if (length(act) == 0) break
    E   <- sweep(draw(reps)[act, , drop = FALSE], 2, mu, "-")
    z[act] <- theta * rowMeans(E) + (1 - theta) * z[act]
    R   <- apply(E, 1, max) - apply(E, 1, min)
    hit <- abs(z[act]) > ucl | R > uclR | R < lclR
    rl[act[hit]] <- s
    act <- act[!hit]
  }
  rl
}

t2_rl <- function(ucl, seed, reps = n_sims) {
  set.seed(seed)
  rl  <- rep(max_rl, reps)
  act <- seq_len(reps)
  for (s in seq_len(max_rl)) {
    if (length(act) == 0) break
    Yk <- draw(reps)[act, , drop = FALSE]
    b1 <- as.vector(Yk %*% Xc) / Sxx
    b0 <- rowMeans(Yk) - b1 * mean(X)
    d0 <- b0 - zbar[1]
    d1 <- b1 - zbar[2]
    T2 <- Si[1, 1] * d0^2 + 2 * Si[1, 2] * d0 * d1 + Si[2, 2] * d1^2
    hit <- T2 > ucl
    rl[act[hit]] <- s
    act <- act[!hit]
  }
  rl
}

ucl_T2 <- function(alpha) (p * (m - 1) / (m - p)) * qf(1 - alpha, p, m - p)

# Bisection: finds the value where ARL0 = target. 'increasing' is TRUE when
# ARL0 grows with the argument (L), FALSE when it shrinks (alpha on log scale).
solve_bisect <- function(arl, lo, hi, tol, increasing) {
  while (hi - lo > tol) {
    mid <- (lo + hi) / 2
    if ((arl(mid) > target) == increasing) hi <- mid else lo <- mid
  }
  (lo + hi) / 2
}

# ==============================================================================
# RUN THE CALIBRATION
# ==============================================================================
seed_cal <- 2024

# 1. EWMA-R multiplier L (union chart)
L <- solve_bisect(function(L) mean(ewma_rl(L, seed_cal)), 1, 20, 1e-3, TRUE)

# 2. T^2 significance level alpha (searched on the log scale)
log_alpha <- solve_bisect(function(la) mean(t2_rl(ucl_T2(exp(la)), seed_cal)),
                          log(1e-5), log(0.2), 1e-4, FALSE)
alpha  <- exp(log_alpha)

# 3. Limits in the form used by main.R
UCL_T2   <- ucl_T2(alpha)
UCL_EWMA <- L * sig * sqrt(theta / ((2 - theta) * n))
UCL_R    <- sig * (d2 + L * d3)
LCL_R    <- max(0, sig * (d2 - L * d3))

# 4. Verification with an independent seed
rl_e  <- ewma_rl(L, 4096, n_ver)
rl_t  <- t2_rl(UCL_T2, 4096, n_ver)
check <- data.frame(
  chart = c("EWMA-R", "T2"),
  ARL0  = c(mean(rl_e), mean(rl_t)),
  SE    = c(sd(rl_e), sd(rl_t)) / sqrt(n_ver)
)

cat(sprintf("\nL = %.4f   alpha = %.6f   UCL_T2 = %.3f\n", L, alpha, UCL_T2))
print(check, row.names = FALSE)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  a0_ref = a0_ref, a1_ref = a1_ref, sigma_ref = sig,
  L = L, alpha = alpha, UCL_T2 = UCL_T2,
  UCL_EWMA = UCL_EWMA, UCL_R = UCL_R, LCL_R = LCL_R,
  ARL0_EWMA_R = check$ARL0[1], ARL0_T2 = check$ARL0[2]
), "results/tables/calibration_limits.csv", row.names = FALSE)
