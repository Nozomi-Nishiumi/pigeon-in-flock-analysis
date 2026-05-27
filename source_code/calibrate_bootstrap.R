# =============================================================================
# Block Bootstrap Calibration via Simulation
# =============================================================================
#
# Validates the calibration of the block bootstrap two-sample test under the
# null hypothesis (H0: no difference between groups), across a range of block
# sizes.
#
# For each block size, Ntest pairs of correlated time series are simulated from
# the same Ornstein-Uhlenbeck (OU) process (so H0 is true by construction).
# The block bootstrap test is applied to each pair and p-values are recorded.
#
# Under a perfectly calibrated test, p-values are uniform on [0,1], meaning
# the k-th quantile of p-values equals k. Deviations from the diagonal (y = x)
# reveal miscalibration: curves above the diagonal indicate a conservative test
# (true FPR < nominal), curves below indicate a liberal test (true FPR > nominal).
#
# Output: bootstrap_calibration.png
#   First panel: Each coloured line shows the empirical p-value quantiles vs. the
#   false-positive-rate threshold for one block size.
#
#   Second panel: Required p-value to control FP at 5% for each block size
#
# Dependencies: mvtnorm, viridis, Rcpp (all CRAN), bootstrap_fast.cpp (same directory)
# =============================================================================

library(mvtnorm)
library(viridis)
library(parallel)
library(Rcpp)

# Compile the fast C++ bootstrap (first run only; Rcpp caches the result)
sourceCpp("bootstrap_fast.cpp")


# =============================================================================
# Utility: scientific-notation y-axis
# =============================================================================

# Format a single number as a plotmath string in scientific notation.
# e.g. fmt_sci(0.00025) → "2.5 %*% 10^-4"  (renders as 2.5 × 10⁻⁴)
fmt_sci <- function(x) {
  if (x == 0) return("0")
  s <- formatC(x, format = "e", digits = 1)
  s <- sub("e\\+0*(\\d+)$", " %*% 10^\\1",  s)
  s <- sub("e-0*(\\d+)$",   " %*% 10^-\\1", s)
  s
}

# Draw a y-axis with plain coefficient tick labels and a shared "10^n" scale
# indicator positioned just above the top of the axis.
sci_yaxis <- function() {
  usr <- par("usr")[3:4]
  at  <- pretty(usr, n = 5)
  at  <- at[at >= usr[1] & at <= usr[2]]

  nonzero <- at[at != 0]
  exp_val <- floor(log10(max(abs(nonzero))))
  scale   <- 10^exp_val

  coef_labels <- formatC(at / scale, format = "g", digits = 2)
  axis(2, at = at, labels = coef_labels, las = 1)

  # Small scale indicator just above the top of the axis
  scale_label <- parse(text = paste0("x ~ 10^", exp_val))
  text(x = par("usr")[1], y = par("usr")[4],
       labels = scale_label,
       pos = 3, xpd = TRUE, cex = 0.75)
}


# =============================================================================
# Simulation parameters
# =============================================================================

rho     <- 0.99          # per-frame autocorrelation of the OU process
npoints <- 500               # number of time points per simulated series
Nbs     <- 10000             # bootstrap replicates per test
Ntest   <- 10000               # simulated series pairs evaluated per block size
bs      <- seq(10, 300, 10)  # block sizes to evaluate (in frames)

# Quantile grid: false-positive-rate thresholds from 0 % to 10 %
q_grid <- seq(0, 0.1, 0.001)


# =============================================================================
# Pre-compute the OU covariance matrix
# =============================================================================

# K[i,j] = rho^|i-j|, the covariance between time points i and j under an
# OU process with per-step autocorrelation rho.
# This depends only on rho and npoints, so it is computed once outside the loop.
idx <- 1:npoints
K   <- outer(idx, idx, function(i, j) rho^abs(i - j))


# =============================================================================
# Calibration sweep: one block size at a time
# =============================================================================

# Each block size is independent, so the outer loop is parallelised with
# mclapply. bstest_cpp() (the fast C++ implementation) runs inside each
# worker; no nested parallelism arises.
#
# Results per block size b:
#   $p  – Ntest p-values
#   $cb – quantile curve over q_grid
#   $pc – p-value threshold for a true 5 % FPR
run_one_blocksize <- function(b) {
  blocksize <- bs[b]

  # Simulate Ntest pairs of series from the same OU distribution (H0 true)
  X1 <- rmvnorm(Ntest, mean = rep(0, npoints), sigma = K)
  X2 <- rmvnorm(Ntest, mean = rep(0, npoints), sigma = K)

  p_b <- rep(NA, Ntest)
  for (i in 1:Ntest) {
    x1      <- as.vector(X1[i, ])
    x2      <- as.vector(X2[i, ])
    p_b[i]  <- bstest_cpp(x1, x2, blocksize, Nbs)
  }

  list(
    p  = p_b,
    cb = quantile(p_b, q_grid),
    pc = quantile(p_b, 0.05)
  )
}

results <- mclapply(seq_along(bs), run_one_blocksize,
                    mc.cores = 4)

# Unpack results into matrices (row = block size, column = replicate / quantile)
p  <- do.call(rbind, lapply(results, `[[`, "p"))
cb <- do.call(rbind, lapply(results, `[[`, "cb"))
pc <- sapply(results, `[[`, "pc")


# =============================================================================
# Figure: calibration curves
# =============================================================================

png("bootstrap_calibration.png", width = 16, height = 8, units = "cm", res = 300)

layout(matrix(1:2, 1, 2))
par(mar = c(4, 4, 1.5, 0.5))

cols <- viridis(length(bs))

# ---- Panel 1: full calibration curves (p-value quantile vs. false positive rate) ----

plot(q_grid,cb[1, ],
     type = "n", bty = "l",
     xlab = "", ylab = "",
     xlim = range(q_grid),
     ylim = c(0, max(cb, na.rm = TRUE)))
                            
title(xlab = "False positive rate", line = 2)
title(ylab = "p-value",             line = 2)

abline(v= 0.05, lty = 2, col = "grey50")    # conventional α = 0.05 threshold

for (b in seq_along(bs)) {
  lines(q_grid,cb[b, ], col = cols[b])
}

legend("topleft",
       col    = cols,
       legend = bs,
       lty    = 1,
       title  = "Block size (sec)",
       ncol   = 3,
       cex    = 0.4,
       bty    = "o")

# ---- Panel 2: α = 0.05 calibration threshold vs. block size -------------------
# pc[b] is the p-value threshold needed to achieve a true 5 % FPR at block size bs[b].
# A well-calibrated test should give pc ≈ 0.05 (dashed reference line).

plot(bs, pc,
     bty = "l", xlab = "", ylab = "",
     pch = 20, cex = 0.7,
     yaxt = "n")
#sfit = smooth.spline(bs, pc)
#lines(sfit$x, sfit$y, col='blue')
sci_yaxis()
title(xlab = "Block size (seconds)",          line = 2)
title(ylab = "p-value threshold (α = 0.05)",  line = 2.5)

dev.off()
