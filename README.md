# Data and analysis code

This repository contains the raw data and analysis code accompanying the manuscript **"Revisiting the selfish herd reveals threat-perceiving birds fly behind their flock, rather than at the center"**.

We attacked flocks of homing pigeons with a drone and examined where individuals that had experienced a previous attack positioned themselves within the flock.

---

## Directory layout

```
.
├── csv/                 # Raw and intermediate data from the experiment
│   ├── Day1_Pre_drone/        # Day before the drone attack (baseline)
│   │   ├── FL01/  FL02/  FL03/  # Three flocks (4 birds + a virtual centerpoint each)
│   │   │   ├── original GPS data/    # Raw GPS logger output (lon/lat/alt)
│   │   │   ├── earth_frame data/     # Time-synced, smoothed ECEF coordinates (.Rdata)
│   │   │   ├── flock_frame data/     # Positions rotated into the flock frame (.Rdata, .csv)
│   │   │   └── visualized data/      # Visualisation outputs (PNG / HTML)
│   ├── Day2_Drone_interception/  # Day of the drone attack (drone and pigeon trajectories)
│   │   └── 20230319_RT_drone.csv / *_pigeon.csv / *_nozomi.csv / *_takao.csv
│   └── Day3_Post_drone/         # Day after the attack (birds anticipating another attack)
│       └── FL01/  FL02/  FL03/  # Same structure as Day 1
│
├── derived_Rdata/      # Cached intermediate products of the pipeline
│   ├── flockframe_Day1_FL0*.Rdata
│   ├── flockframe_Day3_FL0*.Rdata
│   └── list_perDay.Rdata
│
├── figures/            # Figures used in the manuscript
│   ├── bootstrap/             # Block-bootstrap p-value plots
│   ├── closeness/             # Deviation from the flock centre (time series, boxplots)
│   ├── distribution/          # Position distributions in the flock frame (xy, yz)
│   └── trajectories/          # Flight trajectories on a map (fore-aft / altitude colour-coded)
│
├── materials/          # Map data used for rendering
│   └── map.RData
│
└── source_code/        # Analysis scripts
    ├── R_script.Rproj
    ├── main.R                       # Entry point (runs the full pipeline)
    ├── GPS_to_ECEF.R                # GPS → ECEF, with time syncing and smoothing
    ├── ECEF_to_Flock_frame.R        # ECEF → flock frame (fore-aft y / vertical z)
    ├── visualization.R              # Trajectory and position-distribution plots
    ├── aggregation_analysis.R       # Aggregation analysis of the non-experienced birds (#2–#4)
    ├── bootstrap_test.R             # Two-sample block-bootstrap test
    ├── bootstrap_fast.cpp           # Fast C++ implementation of the test (Rcpp)
    ├── calibrate_bootstrap.R        # Standalone calibration check for the block bootstrap (not part of main.R)
    ├── velocity_analysis.R          # Flight-speed analysis
    └── common_functions.R           # Shared utilities (coordinate transforms, rotations, plotting)
```

## Analysis pipeline

Open `source_code/R_script.Rproj` in RStudio first — this is required so that `here::here()` is anchored correctly (see [How to run](#how-to-run)). Then, using `source_code/main.R` as the entry point, the pipeline runs in the following order.

1. **Coordinate transformation**
   - `GPS_to_ECEF.R`: converts raw GPS (lon/lat/alt) into ECEF, time-syncing all streams to 5 fps and smoothing them.
   - `ECEF_to_Flock_frame.R`: rotates each flock's coordinates using its mean heading to produce the **flock frame**:
     - **x**: lateral (left–right with respect to the heading)
     - **y**: fore–aft (**positive = forward**)
     - **z**: vertical (**positive = up**)

     Each flock contains five series:
     - `#1` … the target bird that experienced the drone attack (drone-experienced)
     - `#2`–`#4` … the other birds in the same flock (flockmates)
     - `#5(centerpoint)` … the centroid of `#2`–`#4` (a virtual centre point)

     The output is written to `flock_frame data/` as CSV and Rdata.
2. **Visualisation**
   - `visualization.R`: plots flight trajectories on a map, xy/yz position distributions in the flock frame, and 3D scenes.
3. **Statistical testing** (`bootstrap_test.R::full_analysis()`)
   - Compares the target bird `#1` between Day 1 and Day 3 on **y** (fore–aft position) and **z** (altitude) using a block bootstrap.
   - Compares the aggregation of `#2`–`#4` (mean distance to their own centroid; `ave.dist_center_of_three`) between Day 1 and Day 3 in the same way.
   - Runs both the per-flock analysis (FL01/02/03) and a combined multi-stream analysis across the three flocks, sweeping block size to check the stability of the p-values.
4. **Additional aggregation analysis**
   - `aggregation_analysis.R`: plots the deviation from the flock centre over time, against distance to the drone interception area, and as boxplots.

### Supplementary: bootstrap calibration check

`calibrate_bootstrap.R` is **not** sourced by `main.R` and is **not required** to reproduce any result in the manuscript. It is provided for readers who want to verify the calibration of the block-bootstrap test themselves: it simulates pairs of correlated time series from an Ornstein–Uhlenbeck process under H₀ and checks that the resulting p-values are uniform on [0, 1] across a range of block sizes. It is heavy (Monte-Carlo simulation × many block sizes × `bstest_cpp`) and may take a long time. Run it standalone from `source_code/` only if you need this methodological cross-check.

---

## How to run

### Requirements

- R 4.x or later
- A C++ toolchain (needed by Rcpp; Xcode Command Line Tools on macOS, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, `build-essential` or equivalent on Linux). If you have ever built an Rcpp-using package on this machine, you already have it.

The required R packages are installed automatically (where missing) at the top of `source_code/common_functions.R`.

### Steps

1. Clone the repository.
   ```sh
   git clone <repo-url>
   ```
2. **Open `source_code/R_script.Rproj` in RStudio.**
   Each script resolves the repository root via `dirname(here::here())`. Opening the `.Rproj` anchors `here::here()` to `source_code/`, so relative paths to `csv/`, `derived_Rdata/`, `figures/` etc. resolve correctly. If you are not using RStudio, start R with `source_code/` as the working directory instead.
3. Open `source_code/main.R` and `Source` it (or run `Rscript source_code/main.R`).

Keep the directory layout as-is (`source_code/` at the repository root, with `csv/`, `derived_Rdata/`, `figures/` etc. as siblings).

The block bootstrap (`bootstrap_test.R::full_analysis()`) defaults to 10⁶ replicates and runs in parallel across (number of cores − 1) via `parallel::mclapply`. A single call to `full_analysis()` may take tens of minutes. For a quick sanity check, pass a smaller `Nbs` (e.g. `1e4`).

---

## Data provenance

- Location: around the University of Georgia, Georgia, USA
- Date: March 2023
- Flock size: 4 birds per flock (5 series including the virtual centerpoint `#5`)
- Drone attack: carried out against each flock on Day 2, targeting only the focal bird `#1`
- GPS: time-synced and smoothed to 5 Hz before analysis

The coordinate `drone_engage = (-83.3928, 33.8147)` marks the drone's interception point and is used as the reference for "distance to the drone area" in `aggregation_analysis.R` and elsewhere.

---

## Contact

Dr. Nozomi Nishiumi
Graduate School of Science and Technology, Niigata University
8050 Ikarashi 2 no-cho, Nishi-ku, Niigata City, Niigata 950-2181, Japan
nishiumi@create.niigata-u.ac.jp

Co-authors: Richard Mann, Bennett Caughran, Bingnan Li, Takao Sasaki
