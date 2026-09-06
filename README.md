# Nanobubble tracking with motion-predicted assignment

This repository contains a MATLAB tracker for intermittent ultrasound localization lists. It is designed for the tracking stage only: it receives per-frame `[z, x]` coordinates and connects detections into trajectories. It does not run SVD, beamforming, IQ processing, or localization.

The tracker combines a constant-velocity Kalman prediction with global Hungarian assignment. Before each frame, it predicts each active trajectory's position and uncertainty. Hungarian assignment then chooses exclusive detection-to-track links using the predicted innovation cost, with explicit unmatched alternatives. Tracks can survive short acquisition gaps, and a score-based evidence check suppresses trajectories formed by random detections. The original `munkres.m` implementation is included in `src/` so the core does not require MATLAB's Sensor Fusion and Tracking Toolbox.

## Result

On a frozen tracking-only holdout of 100 primary trials, using the same localization lists for every method, the new tracker reached **95.5/100 detection-conditioned IDF1**. The original tracker reached **71.4/100** when given a longer gap allowance and **17.3/100** with its current setting. The new tracker therefore improved the primary score by **24.0 percentage points** beyond that gap-only baseline (paired bootstrap 95% interval: 23.1–24.9 points).

These are simulated coordinates with known identities, not measured physical nanobubble accuracy. A 128-frame subset of real L22 2 kHz saved localizations was used only to verify the MATLAB workflow; it has no independent identity labels.

The detailed methods, condition-by-condition table, figures, ablations, failure cases, and validation record are in [`results/RESULTS.md`](results/RESULTS.md).

## Quick start in MATLAB

Add `src/` to the MATLAB path and pass the saved localization cells directly:

```matlab
addpath('/path/to/nanobubble-tracking/src');
[tracks, adjacency_tracks, details] = nbtracker(SR_Localizations, ...
    'FrameRateHz', 2000, ...
    'PixelSizeUm', [12.32 12.41446725317693]);
```

`SR_Localizations` must be a cell vector. Each cell is an `N-by-2` numeric array of `[z, x]` pixel coordinates for one frame. `adjacency_tracks` uses the same one-based global-index convention as the original SRU `simpletracker`: it indexes the vertical concatenation of all input detections. `tracks` contains one per-frame index vector per retained trajectory, with `NaN` during missed observations. Predictions through gaps are never counted as observations.

For saved MAT files, use the bounded wrapper. It loads only localization variables, preserves `meta.dsRange` or saved timestamps, and defaults to at most 192 frames:

```matlab
addpath('/path/to/nanobubble-tracking/src');
result = track_nb_file('/path/to/localizations.mat', ...
    'PixelSizeUm', [12.32 12.41446725317693], ...
    'FrameRateHz', 2000, ...
    'StartFrame', 1, ...
    'MaxFrames', 192, ...
    'OutputMAT', '/path/to/new_nb_tracks.mat');
```

The wrapper refuses to overwrite an existing output. For the assessed L22 grid, the axial spacing is 12.32 µm and the interpolated lateral spacing is 12.414467 µm. Use the calibration for the actual saved coordinate grid if it differs. The native lateral spacing before interpolation is 100 µm.

The example function accepts any saved localization file:

```matlab
addpath('/path/to/nanobubble-tracking/src');
addpath('/path/to/nanobubble-tracking/examples');
result = demo_nb_tracking('/path/to/localizations.mat');
```

The original SRU tuning UI is not modified; call `nbtracker` or `track_nb_file` explicitly when you want the new behavior.

## Tests

With MATLAB R2025a or a compatible release:

```matlab
addpath('/path/to/nanobubble-tracking/src');
addpath('/path/to/nanobubble-tracking/tests');
test_nbtracker;
test_file_workflow;
```

The tests cover empty and singleton frames, gap boundaries, irregular timestamps, anisotropic pixel spacing, physical velocity, crossing identities, immutable input coordinates, exclusive assignments, invalid inputs, and the bounded file wrapper.

The scorer has independent known-answer checks:

```bash
python -m pip install -r requirements.txt
python benchmark/test_scoring.py
```

## Reproducing the benchmark

The benchmark is deliberately separate from experimental data. It generates 140 fixed-localization scenarios (14 conditions × 10 seeds), evaluates the original and new trackers in MATLAB, and scores the output in Python.

```bash
python benchmark/generate_holdout.py
```

Run the native MATLAB evaluation from the repository's `benchmark/` folder:

```matlab
cd('/path/to/nanobubble-tracking/benchmark');
run_holdout;
```

Then audit and score the results:

```bash
python benchmark/analyze_holdout.py
```

The generator records a SHA-256 hash of `src/nbtracker.m` before generating the holdout. If the candidate changes, it refuses to reuse the frozen evaluation. The generated MAT files are intentionally ignored by Git because they are large; the frozen protocol, manifest, metrics, summary, and figures are included for review.

## Scope and limitations

The primary metric is identity preservation on supplied detections. It excludes missing localizations from the denominator so that it evaluates association separately from localization. It does not establish the fraction of physical bubbles detected or tracked in an experiment. Persistent stationary false peaks remain difficult when the tracker receives coordinates alone. The benchmark is not a FIELD II acoustic simulation and does not model nonlinear bubble response, shell dynamics, transmit pressure, or experimental clutter calibration.

The code assumes MATLAB's standard language features and does not require Parallel Computing, Image Processing, Signal Processing, Statistics, or Sensor Fusion and Tracking Toolbox licenses. The benchmark's Python dependencies are listed in [`requirements.txt`](requirements.txt).

## Attribution

The Hungarian assignment helper included in `src/munkres.m` is the Yi Cao vectorized MATLAB implementation (version 2.3, 2011), whose original header and reference are retained. The tracker and benchmark code in this repository were developed for the SRU nanobubble tracking assessment described in `results/RESULTS.md`.
