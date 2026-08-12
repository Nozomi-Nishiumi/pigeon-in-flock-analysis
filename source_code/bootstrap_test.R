# =============================================================================
# Bootstrap Analysis of Drone Flock Position Data
# =============================================================================
#
# Tests whether a drone fly-over caused a significant change in flock
# behaviour by comparing pre- and post-event position time series using
# a block bootstrap hypothesis test.
#
# The test is run across a range of block sizes (in frames, converted to
# seconds on the plot) to assess robustness: p-values that remain small
# and stable across block sizes provide stronger evidence of an effect.
#
# Four panels are produced:
#   Panel 1  – Combined multi-stream analysis (all 3 flocks pooled)
#   Panels 2–4 – Individual flock analyses (Flocks 1, 2, 3)
#
# Output: drone_bootstrap_analysis.png
#
# Dependencies: Rcpp (CRAN), bootstrap_fast.cpp (same directory)
# =============================================================================

library(parallel)
library(Rcpp)

# Compile the fast C++ bootstrap (first run only; Rcpp caches the result)
sourceCpp(paste0(basedir,"/source_code/bootstrap_fast.cpp"))

# Path to processed flock CSV files (rotated/transformed coordinates)
datapath <- paste0(basedir,"/csv/")


# =============================================================================
# Utility: scientific-notation axis labels for base-R plots
# =============================================================================

# Format a single numeric value as a plotmath string in scientific notation.
# Returns a string suitable for use with parse(text = ...).
# Example: fmt_sci(0.00025) → "2.5 %*% 10^-4"  (renders as 2.5 × 10⁻⁴)
fmt_sci <- function(x) {
  if (x == 0) return("0")
  s <- formatC(x, format = "e", digits = 1)
  s <- sub("e\\+0*(\\d+)$", " %*% 10^\\1",  s)   # e.g. e+04 → × 10^4
  s <- sub("e-0*(\\d+)$",   " %*% 10^-\\1", s)   # e.g. e-04 → × 10^-4
  s
}

# Draw a y-axis with a shared scale indicator at the top.
# Tick labels show plain coefficients (e.g. 0, 1, 2) and a single
# "× 10^n" label is placed at the top of the axis to give the common scale.
# This keeps labels short and avoids overlap.
sci_yaxis <- function() {
  usr <- par("usr")[3:4]                         # current y-axis limits
  at  <- pretty(usr, n = 3)
  at  <- at[at >= usr[1] & at <= usr[2]]         # drop ticks outside plot area

  # Derive the common exponent from the largest non-zero tick value
  nonzero  <- at[at != 0]
  exp_val  <- floor(log10(max(abs(nonzero))))
  scale    <- 10^exp_val

  # Plain coefficient labels (trailing zeros removed via "g" format)
  coef_labels <- formatC(at / scale, format = "g", digits = 2)
  axis(2, at = at, labels = coef_labels, las = 1)

  # Scale indicator: small superscript just above the top of the y-axis.
  # text() with xpd=TRUE lets us draw outside the plot region.
  # adj=c(1.1, -0.4) right-aligns it against the axis and floats it just above usr[2].
  scale_label <- parse(text = paste0("10^", exp_val))
  text(x = par("usr")[1], y = par("usr")[4],
       labels = scale_label,
       pos=3,
       xpd = TRUE, cex = 0.75)
}

# Build a plotmath expression for a "max p = X" annotation in scientific notation.
# Example: fmt_maxp(2e-4) → expression for "max p = 2.0 × 10⁻⁴"
fmt_maxp <- function(p) {
  s <- fmt_sci(p)
  parts <- strsplit(s, " %\\*% 10\\^")[[1]]
  coef <- parts[1]
  exp  <- parts[2]
  
  bquote("max p = " ~ .(coef) %*% 10^.(as.numeric(exp)))
}# ensure consistent scientific notation (NN)


# =============================================================================
# Data I/O
# =============================================================================

# Read position data for a given flock and experimental period.
#
# Args:
#   flockid      – flock ID (integer: 1, 2, or 3)
#   birdid       – bird  ID (string: #1-#4, or #5(centerpoint) )
#   prepost – "Day1_Pre_drone"  (before drone) or "Day3_Post_drone" (after drone)
#   type    – which coordinate to extract:
#               "z"   → z  (altitude)
#               "y"   → y  (front-back position)
#               "centre"   → Euclidean distance from flock origin, sqrt(x²+y²)
#
# Returns: numeric vector of the requested position measurements.
read_in_flock <- function(flockid, birdid, prepost = "Day1_Pre_drone", type = "height") {
  day <- strsplit(prepost, "_")[[1]][1]
    filename <- paste0(
      basedir,"/csv/",
      prepost,#Day
      "/FL0", flockid,#Flock
      "/flock_frame data/flockframe_", day,"_FL0", flockid, "_",birdid,".csv"
    )
 
  d <- read.csv(filename, header = TRUE)
  return(d[[type]])
}


# =============================================================================
# Bootstrap hypothesis tests
# =============================================================================

# Block bootstrap two-sample test for a single time series.
#
# Thin R wrapper around bstest_cpp() (bootstrap_fast.cpp).
# Uses the geometric block bootstrap with a prefix-sum implementation
# that avoids assembling resampled series, giving a large speed gain
# over tsboot() for the mean statistic.
#
# Args:
#   x1, x2    – numeric vectors (pre- and post-event time series)
#   blocksize – expected block length (geometric distribution parameter)
#   Nbs       – number of bootstrap replicates (default 10,000)
#
# Returns: two-sided p-value (scalar)
bstest <- function(x1, x2, blocksize, Nbs = 10000) {
  bstest_cpp(x1, x2, blocksize, Nbs)
}


# Block bootstrap two-sample test combining three simultaneous flock streams.
#
# Thin R wrapper around bstest_multistream_cpp() (bootstrap_fast.cpp).
# Each stream is bootstrapped independently; the test statistic is the
# sum of bootstrap means across all three flocks.
#
# Args:
#   X1, X2     – named lists with elements $a, $b, $c (one time series per flock)
#   blocksize, Nbs – as in bstest()
#
# Returns: two-sided p-value (scalar)
bstest_multistream <- function(X1, X2, blocksize, Nbs = 10000) {
  bstest_multistream_cpp(X1$a, X1$b, X1$c, X2$a, X2$b, X2$c, blocksize, Nbs)
}


# =============================================================================
# Block-size scan: sweep over block sizes to assess robustness of results
# =============================================================================

# Compute p-values across a range of block sizes (single flock stream).
#
# Block sizes are evaluated in parallel via mclapply. Each worker calls
# bstest() (the fast C++ implementation), so no nested parallelism arises.
#
# Returns a list:
#   $bs – block sizes tested (in frames)
#   $p  – corresponding two-sided p-values
scan_p_blocksize <- function(x1, x2, bs = seq(50, 1500, 50), Nbs = 10000) {
  p <- mclapply(bs, function(blocksize) bstest(x1, x2, blocksize, Nbs),
                mc.cores = min(length(bs), max(1, detectCores() - 1)))
  return(list(bs = bs, p = unlist(p)))
}


# As above, but for the combined multi-stream (all three flocks) analysis.
scan_p_blocksize_multistream <- function(X1, X2, bs = seq(50, 1500, 50),
                                         Nbs = 10000) {
  p <- mclapply(bs, function(blocksize) bstest_multistream(X1, X2, blocksize, Nbs),
                mc.cores = min(length(bs), max(1, detectCores() - 1)))
  return(list(bs = bs, p = unlist(p)))
}


# =============================================================================
# Main analysis and figure
# =============================================================================

# Run the full bootstrap block-size scan and produce a 4-panel PNG figure.
#
# Each panel plots p-value (y, scientific notation) vs block size in seconds (x)
# for a different analysis:
#   Panel 1: Overall – combined multi-stream (Flocks 1–3 pooled)
#   Panels 2–4: Individual Flocks 1, 2, 3
#
# A dashed horizontal line marks the conventional α = 0.05 threshold.
# The maximum p-value across block sizes is annotated below each panel.
#
# Args:
#   Nbs – number of bootstrap replicates (default 1,000,000)
full_analysis <- function(Nbs = 1e6, birdid="#1", type = "y") {
  set.seed(123)

  savedir=paste0(basedir,"/figures/statistics")
  dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
  setwd(savedir)
  graphics.off()
  png(paste0("drone_bootstrap_analysis_",type,".png"), width = 15, height = 5,
      units = "cm", res = 300)

  layout(matrix(1:4, 1, 4, byrow = TRUE))
  # mar: (bottom, left, top, right) inner margins in lines
  # oma: (bottom, left, top, right) outer margins in lines
  par(mar = c(3, 3, 1.5, 0.5), oma = c(2, 0, 1.1, 0))

  # ---- Panel 1: combined multi-stream analysis (all three flocks) ------------

  # Load pre-drone data for all three flocks
  flock1_pre  <- read_in_flock(1,birdid, "Day1_Pre_drone",  type = type)
  flock2_pre  <- read_in_flock(2,birdid, "Day1_Pre_drone",  type = type)
  flock3_pre  <- read_in_flock(3,birdid, "Day1_Pre_drone",  type = type)
  X1 <- list(a = flock1_pre, b = flock2_pre, c = flock3_pre)

  # Load post-drone data for all three flocks
  flock1_post <- read_in_flock(1,birdid, "Day3_Post_drone", type = type)
  flock2_post <- read_in_flock(2,birdid, "Day3_Post_drone", type = type)
  flock3_post <- read_in_flock(3,birdid, "Day3_Post_drone", type = type)
  X2 <- list(a = flock1_post, b = flock2_post, c = flock3_post)

  # Scan p-values over block sizes (block sizes divided by 5 to convert frames → seconds)
  tmp <- scan_p_blocksize_multistream(X1, X2, Nbs = Nbs)
  bs  <- tmp$bs
  p   <- tmp$p

  plot(bs / 5, p,
       bty = "l", xlab = "", ylab = "",
       pch = 20, cex = 0.7,
       ylim = c(0, 1.5 * max(p)), yaxs = "i",
       yaxt = "n")                          # suppress default axis; draw custom below
  sci_yaxis()
  title(xlab = "Block size (seconds)", line = 2)
  title(ylab = "p-value",              line = 2)
  abline(h = 0.05, lty = 2, col = "grey60")  # α = 0.05 reference line

  maxp_overall <- signif(max(p), 2)

  # ---- Panels 2–4: individual flock analyses ---------------------------------

  maxp_flock <- rep(NA, 3)

  for (flockid in 1:3) {
    x1 <- read_in_flock(flockid,birdid, "Day1_Pre_drone",  type = type)
    x2 <- read_in_flock(flockid,birdid, "Day3_Post_drone",  type = type)

    tmp <- scan_p_blocksize(x1, x2, Nbs = Nbs)
    bs  <- tmp$bs
    p   <- tmp$p

    plot(bs / 5, p,
         bty = "l", xlab = "", ylab = "",
         pch = 20, cex = 0.7,
         ylim = c(0, 1.5 * max(p)), yaxs = "i",
         yaxt = "n")
    sci_yaxis()
    title(xlab = "Block size (seconds)", line = 2)
    title(ylab = "p-value",              line = 2)

    maxp_flock[flockid] <- signif(max(p), 2)
  }

  # ---- Outer-margin annotations (panel titles and max-p values) --------------

  # Panel centres are evenly spaced: 1/8, 3/8, 5/8, 7/8 of figure width
  panel_mid <- seq(1/8, 7/8, 1/4) + 0.03

  # Panel titles (top outer margin)
  mtext("Overall", side = 3, outer = TRUE, at = panel_mid[1],
        line = 0.1, cex = 0.7, font = 2)
  mtext("Flock 1", side = 3, outer = TRUE, at = panel_mid[2],
        line = 0.1, cex = 0.7, font = 2)
  mtext("Flock 2", side = 3, outer = TRUE, at = panel_mid[3],
        line = 0.1, cex = 0.7, font = 2)
  mtext("Flock 3", side = 3, outer = TRUE, at = panel_mid[4],
        line = 0.1, cex = 0.7, font = 2)

  # Max p-value annotations (bottom outer margin), formatted in scientific notation
  mtext(fmt_maxp(maxp_overall),  side = 1, outer = TRUE, at = panel_mid[1],
        line = 0.5, cex = 0.7, font = 2)
  mtext(fmt_maxp(maxp_flock[1]), side = 1, outer = TRUE, at = panel_mid[2],
        line = 0.5, cex = 0.7, font = 2)
  mtext(fmt_maxp(maxp_flock[2]), side = 1, outer = TRUE, at = panel_mid[3],
        line = 0.5, cex = 0.7, font = 2)
  mtext(fmt_maxp(maxp_flock[3]), side = 1, outer = TRUE, at = panel_mid[4],
        line = 0.5, cex = 0.7, font = 2)

  dev.off()
}


