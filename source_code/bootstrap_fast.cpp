// =============================================================================
// Fast geometric block bootstrap for the mean statistic
// =============================================================================
//
// Implements the stationary (geometric) block bootstrap of Politis & Romano
// (1994), matching the behaviour of tsboot(sim = "geom") from the R boot
// package, but restricted to the mean statistic.
//
// Speed gain comes from two sources:
//   1. Prefix-sum trick: cumulative sums are pre-computed so any circular
//      block sum is a single O(1) subtraction — no data copying needed.
//   2. C++ inner loop: avoids R interpreter overhead across all Nbs replicates.
//
// Exported functions (callable from R via Rcpp::sourceCpp):
//   bstest_cpp()             – two-sample test, single stream
//   bstest_multistream_cpp() – two-sample test, three streams combined
//
// =============================================================================

#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Internal helper: generate Nbs bootstrap means for one time series.
//
// Args:
//   x   – time series of length n
//   Nbs – number of bootstrap replicates
//   l   – mean block length (geometric parameter; block length = rgeom(1/l) + 1)
// -----------------------------------------------------------------------------
static NumericVector geom_block_means(const NumericVector& x, int Nbs, double l) {
  int    n    = x.size();
  double prob = 1.0 / l;

  // Prefix sums: cs[i] = sum(x[0], ..., x[i-1]), cs[0] = 0
  NumericVector cs(n + 1, 0.0);
  for (int i = 0; i < n; i++) cs[i + 1] = cs[i] + x[i];
  double x_total = cs[n];

  NumericVector means(Nbs);

  for (int b = 0; b < Nbs; b++) {
    double s     = 0.0;
    int    count = 0;

    while (count < n) {
      // Uniform random start (0-indexed) and geometric block length (mean = l)
      int start = (int)(R::unif_rand() * n);
      int len   = (int)R::rgeom(prob) + 1;
      int take  = (len < n - count) ? len : n - count;  // cap at remaining elements

      // O(1) circular block sum via prefix sums
      int end = start + take;
      if (end <= n) {
        s += cs[end] - cs[start];
      } else {
        // Wraps: tail of series + head of series
        s += (x_total - cs[start]) + cs[end - n];
      }
      count += take;
    }

    means[b] = s / n;
  }

  return means;
}

// -----------------------------------------------------------------------------
// Two-sided block bootstrap test comparing means of two time series.
// Statistical equivalent of bstest() in the R code.
//
// [[Rcpp::export]]
// -----------------------------------------------------------------------------
double bstest_cpp(NumericVector x1, NumericVector x2, int blocksize, int Nbs) {
  RNGScope rng;

  NumericVector z = geom_block_means(x1, Nbs, blocksize);
  NumericVector w = geom_block_means(x2, Nbs, blocksize);

  double p1 = 0.0, p2 = 0.0;
  for (int i = 0; i < Nbs; i++) {
    double d = w[i] - z[i];
    if (d < 0.0) p1 += 1.0;
    if (d > 0.0) p2 += 1.0;
  }
  p1 = (p1 + 1.0) / (Nbs + 1.0);
  p2 = (p2 + 1.0) / (Nbs + 1.0);

  return 2.0 * (p1 < p2 ? p1 : p2);
}

// -----------------------------------------------------------------------------
// Two-sided block bootstrap test combining three simultaneous flock streams.
// Statistical equivalent of bstest_multistream() in the R code.
// Streams are bootstrapped independently; the test statistic is the sum of
// means across all three flocks.
//
// [[Rcpp::export]]
// -----------------------------------------------------------------------------
double bstest_multistream_cpp(NumericVector x1a, NumericVector x1b, NumericVector x1c,
                               NumericVector x2a, NumericVector x2b, NumericVector x2c,
                               int blocksize, int Nbs) {
  RNGScope rng;

  NumericVector za = geom_block_means(x1a, Nbs, blocksize);
  NumericVector wa = geom_block_means(x2a, Nbs, blocksize);
  NumericVector zb = geom_block_means(x1b, Nbs, blocksize);
  NumericVector wb = geom_block_means(x2b, Nbs, blocksize);
  NumericVector zc = geom_block_means(x1c, Nbs, blocksize);
  NumericVector wc = geom_block_means(x2c, Nbs, blocksize);

  double p1 = 0.0, p2 = 0.0;
  for (int i = 0; i < Nbs; i++) {
    double d = (wa[i] + wb[i] + wc[i]) - (za[i] + zb[i] + zc[i]);
    if (d < 0.0) p1 += 1.0;
    if (d > 0.0) p2 += 1.0;
  }
  p1 = (p1 + 1.0) / (Nbs + 1.0);
  p2 = (p2 + 1.0) / (Nbs + 1.0);

  return 2.0 * (p1 < p2 ? p1 : p2);
}
