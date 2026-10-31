# Drone flock analysis: compare pre- and post-drone exposure GPS metrics
# across three flocks using ESS-adjusted Welch t-tests and Fisher's method.

library(coda)


# ── Data reading ──────────────────────────────────────────────────────────────


phase_dirs <- c(pre = "Day1_Pre_drone", post = "Day3_Post_drone")

# Read a processed GPS flock file and return the requested metric vector.
# type: focal_3d_dist = 3-D distance from flock centre; focal_z_dist = vertical
#        (z) distance from flock centre; focal_y_dist = front-back (y) distance from flock centre;
#       nt_ave_3d_dist = non-threatened birds' average 3-D distance from their centre.
read_in_flock <- function(flockid, birdid, prepost = "pre", type = "focal_z_dist") {
  if (!prepost %in% names(phase_dirs))
    stop("Unknown prepost: ", prepost, " (use 'pre' or 'post')")
  phase_dir <- phase_dirs[[prepost]]
  day <- strsplit(phase_dir, "_")[[1]][1]
  filename <- paste0(
    basedir,"/csv/",
    phase_dir,#Day
    "/FL0", flockid,#Flock
    "/flock_frame data/flockframe_", day,"_FL0", flockid, "_for_statistics.csv"
  )
  d <- read.csv(filename, header = TRUE)
  switch(type,
         focal_z_dist   = d$focal_z_dist,
         focal_y_dist = d$focal_y_dist,
         focal_3d_dist   = sqrt(d$focal_x^2 + d$focal_y_dist^2 + d$focal_z_dist^2),
         nt_ave_3d_dist      = d$nt_ave_3d_dist,
         stop("Unknown type: ", type)
  )
}

# ── Summary statistics ────────────────────────────────────────────────────────

# Compute ESS-adjusted mean and 95% CI for one time series.
summarise_series <- function(x) {
  n      <- effectiveSize(x)
  mu     <- mean(x)
  t_crit <- qt(0.975, df = n - 1)
  sem    <- sd(x) / sqrt(n)
  list(
    len    = length(x),
    acf5   = acf(x, plot = FALSE)$acf[6],  # lag-5 autocorrelation (acf[1] = lag 0)
    mean   = mu,
    ESS    = n,
    CIlow  = mu - t_crit * sem,
    CIhigh = mu + t_crit * sem
  )
}

get_data_properties <- function() {
  types  <- c("focal_3d_dist", "focal_z_dist", "focal_y_dist", "nt_ave_3d_dist")
  phases <- c("pre", "post")
  ids    <- 1:3
  
  rows <- vector("list", length(types) * length(ids) * length(phases))
  idx  <- 0L
  
  for (type in types) {
    for (id in ids) {
      for (prepost in phases) {
        idx <- idx + 1L
        x   <- read_in_flock(id, prepost = prepost, type = type)
        s   <- summarise_series(x)
        rows[[idx]] <- data.frame(
          type    = type,
          prepost = prepost,
          flockid = id,
          len     = s$len,
          acf     = s$acf5,
          mean    = s$mean,
          ESS     = s$ESS,
          CIlow   = s$CIlow,
          CIhigh  = s$CIhigh
        )
      }
    }
  }
  do.call(rbind, rows)
}







T_test_ESS_analysis<-function(){
  df <- get_data_properties()
  tabledir=paste0(basedir,"/tables/statistics")
  dir.create(tabledir, recursive = TRUE, showWarnings = FALSE)
  savedir=paste0(basedir,"/figures/statistics")
  dir.create(savedir, recursive = TRUE, showWarnings = FALSE)
  setwd(savedir)
  # ── Welch t-test on effective samples ─────────────────────────────────────────
  
  # Recover SE from the stored CI bounds (reverse of CI construction above).
  df$se <- ifelse(
    df$ESS >= 2,
    (df$CIhigh - df$CIlow) / (2 * qt(0.975, df = df$ESS - 1)),
    NA_real_
  )
  
  # Two-sided Welch t-test using ESS as the degrees-of-freedom proxy for each series.
  welch_ess_p <- function(row_pre, row_post) {
    if (is.na(row_pre$se) || is.na(row_post$se)) return(NA_real_)
    se2_pre  <- row_pre$se^2
    se2_post <- row_post$se^2
    t_stat   <- (row_pre$mean - row_post$mean) / sqrt(se2_pre + se2_post)
    # Welch-Satterthwaite degrees of freedom
    nu <- (se2_pre + se2_post)^2 /
      (se2_pre^2  / (row_pre$ESS  - 1) +
         se2_post^2 / (row_post$ESS - 1))
    2 * pt(-abs(t_stat), df = nu)
  }
  
  df$welch_p <- NA_real_
  for (t in unique(df$type)) {
    for (fid in unique(df$flockid)) {
      pre  <- df[df$type == t & df$flockid == fid & df$prepost == "pre",  ]
      post <- df[df$type == t & df$flockid == fid & df$prepost == "post", ]
      if (nrow(pre) == 0 || nrow(post) == 0) next
      # Write the same p-value to both pre and post rows for this type/flock
      df[df$type == t & df$flockid == fid, "welch_p"] <- welch_ess_p(pre, post)
    }
  }
  df$se <- NULL  # temporary column, no longer needed
  
  # Fisher's method: combine the three per-flock p-values into one omnibus test.
  # Use "pre" rows only to avoid counting each p-value twice (pre and post carry
  # the same value).
  fisher_p <- sapply(unique(df$type), function(t) {
    p_vals <- na.omit(df[df$type == t & df$prepost == "pre", "welch_p"])
    if (length(p_vals) == 0) return(NA_real_)
    pchisq(-2 * sum(log(p_vals)), df = 2 * length(p_vals), lower.tail = FALSE)
  })
  names(fisher_p) <- unique(df$type)
  df$omnibus_p <- unname(fisher_p[df$type])   # same value on every row of a metric
  
  # ── Table S1: per-series summary plus the test results shown in the figure ──
  # welch_p   - ESS-adjusted Welch t-test, pre vs post, per metric and flock
  #             (identical on the pre and post rows of a pair)
  # omnibus_p - Fisher combination of the three per-flock welch_p values
  #             (identical on all six rows of a metric)
  write.csv(df, file = file.path(tabledir, "TableS1_drone_data_properties.csv"), row.names = FALSE)
  
  # ── Plot ──────────────────────────────────────────────────────────────────────
  
  png("drone_analysis_plot.png", width = 15, height = 15, units = "cm", res = 500)
  
  types  <- unique(df$type)
  types=types[c(3,1,2,4)]
  flocks <- sort(unique(df$flockid))
  colors <- c(pre = "#2196F3", post = "#E53935")
  
  # layout() fills column-by-column, so panel order is top-left, bottom-left,
  # top-right, bottom-right — hence the A, C, B, D label sequence below.
  layout(matrix(1:4, 2, 2))
  par(mar = c(3.0, 3.5, 2.8, 0.8), mgp = c(2.2, 0.5, 0), tcl = -0.3,
      cex.axis = 0.7, cex.lab = 0.8)
  
  width  <- 0.15
  titles <- c(
    "A: The threatened bird's position 
    along front-back axis",
    "C: The threatened bird's distance 
  from the flock center",
    "B: The threatened bird's position 
    along high-low axis",
    "D: Non-threatened flockmates'
  distance from their own center"
  )
  
  fmt <- function(p) {
    if (is.na(p)) return("p=NA")
    if (p < 0.001) "p<0.001" else paste0("p=", formatC(p, digits = 2, format = "g"))
  }
  
  for (count in seq_along(types)) {
    t         <- types[count]
    sub       <- df[df$type == t, ]
    x_centres <- seq_along(flocks) - 1
    x_pre     <- x_centres - width
    x_post    <- x_centres + width
    
    # Base y-axis on finite CI values only, then add padding
    ci_vals <- c(sub$CIlow, sub$CIhigh)
    ci_vals <- ci_vals[is.finite(ci_vals) & !is.na(ci_vals)]
    ylim    <- range(c(sub$mean, ci_vals))
    ylim    <- ylim + diff(ylim) * c(-0.05, 0.18)
    
    plot(NA, xlim = c(-0.5, max(x_centres) + 0.5), ylim = ylim,
         xaxt = "n", xlab = "", ylab = "", bty = "l",
         main = titles[count], cex.main = 0.95)
    mtext(paste0("Fisher omnibus: ", fmt(fisher_p[t])),
          side = 3, line = -0.4, cex = 0.6, col = "grey30")
    
    abline(h = 0, lty = 2, col = "grey70", lwd = 0.8)
    axis(1, at = x_centres, labels = paste0("Flock ", flocks), cex.axis = 0.7)
    
    for (i in seq_along(flocks)) {
      fid      <- flocks[i]
      row_pre  <- sub[sub$flockid == fid & sub$prepost == "pre",  ]
      row_post <- sub[sub$flockid == fid & sub$prepost == "post", ]
      
      # Grey line connecting pre and post means for visual comparison
      segments(x_pre[i], row_pre$mean, x_post[i], row_post$mean,
               col = "grey60", lwd = 1, lty = 1)
      
      for (pp in c("pre", "post")) {
        row   <- if (pp == "pre") row_pre else row_post
        x     <- if (pp == "pre") x_pre[i] else x_post[i]
        ci_lo <- row$CIlow
        ci_hi <- row$CIhigh
        if (!is.na(ci_lo) && !is.na(ci_hi) && is.finite(ci_lo) && is.finite(ci_hi))
          arrows(x, max(ci_lo, ylim[1]), x, min(ci_hi, ylim[2]),
                 angle = 90, code = 3, length = 0.04, col = colors[pp], lwd = 1.5)
        points(x, row$mean, pch = 19, col = colors[pp], cex = 0.9)
      }
      
      # Per-flock p-value label near the top of the panel
      y_top <- ylim[2] - diff(ylim) * 0.02
      text(x_centres[i], y_top, fmt(row_pre$welch_p),
           cex = 0.6, adj = c(0.5, 1), col = "grey20")
    }
    
    if (t == types[1]) {
      legend("bottomright", legend = c("pre-drone", "post-drone"), col = colors,
             pch = 19, pt.cex = 0.9, bty = "n", cex = 0.8, y.intersp = 1.1)
    }
  }
  
  dev.off()
}
