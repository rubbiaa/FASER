# Muon Spectrometer Reconstruction: GenFit vs Taubin — Implementation Report (Oct 26, 2025)

This report summarizes the implementation and validation of a dual-fit muon reconstruction pipeline in the FASER spectrometer codebase, including a Taubin-style circle fit alongside the existing GenFit-based Kalman fitter, expanded output, and plotting/diagnostics.

## Overview

- Added a dual reconstruction path that runs both:
  - GenFit Kalman fit (reference) and
  - Taubin-inspired circle fit in the bending plane.
- Implemented per-hit 100 μm smearing on (x, y) before station averaging to emulate measurement resolution.
- Introduced a geometry-aware magnetic field model with B on between stations and off near detector planes.
- Expanded the ROOT output (TTree "Hits") with explicit `genfit_*` and `taubin_*` branches.
- Added a plotting macro to compare GenFit vs Taubin and export a resolution CSV.

Primary entry points:
- `Display/FaserCalMuonSpectrometer.cxx` — event IO, reconstruction, and tree writing
- `CoreUtils/TMuTrack.hh/.cc` — fit implementations (GenFit + Taubin)
- `Display/plot_muon_genfit_taubin.C` — plots and resolution table

Output ROOT file: `scifi_hits_all.root` (written in the current working directory when running the spectrometer; typical usage is from `Display/build`).

## Build and Run

From the spectrometer build directory (important for runtime library resolution on macOS):

```bash
# Build (from Display/build)
cmake --build . -j8 --target FaserCalMuonSpectrometer

# Run over all runs and all events (recommended to run from Display/build)
./FaserCalMuonSpectrometer 0 0 -

# Run a quick test over a single event
./FaserCalMuonSpectrometer 0 1 -
```

Notes:
- Running the executable outside the build directory can trigger a dyld @rpath error for Boost (macOS). Run from `Display/build` or set up your `DYLD_LIBRARY_PATH` accordingly.
- CLI: `./FaserCalMuonSpectrometer <RunNumber> <MaxEvents> <Mask>`
  - `RunNumber`: 0 means all runs; a positive value filters that run only
  - `MaxEvents`: 0 means all events
  - `Mask`: nueCC, numuCC, nutauCC, nuNC, nuES, or `-` for all

## Magnetic Field Model

Implemented in `CoreUtils/GenMagneticField.hh/.cc` and configured in `FaserCalMuonSpectrometer.cxx`:
- Geometry-driven field gating: |B| ≈ 1.5 T between stations; B ≈ 0 near station planes.
- Direction: Bx with sign by region between stations.
- Per-event station Z positions are provided (from averaged hits) via `SetEventStationZsCm(...)` to robustly gate the field even if geometry is misaligned.
- Field values printed at the first averaged hit for sanity (Tesla printed by converting kGauss from GenFit).

Environment toggles:
- `MS_VERBOSE=1` — extra prints for fitting and field setup
- `MS_DIAG_BDL=1` — enable ΣB⊥·dl diagnostic along the averaged polyline with a quick p estimate

## Measurement Smearing and Station Averaging

- Each raw hit is smeared with a Gaussian σ = 0.1 mm in (x, y).
- Hits belonging to the same station are averaged (position-wise) to form one measurement per station.
- Tracks with fewer than 3 distinct stations are skipped (no fit attempted).

## Fitting Methods

Both fits are executed when `nstations >= 3`:

1) GenFit Kalman
- Uses GenFit with `RKTrackRep`, `KalmanFitterRefTrack`, TGeo material, and the geometry-aware field.
- Seeds at p ≈ 10 GeV/c; initial direction from the first two averaged hits.
- Success criteria: p-value > 0.01.

2) Taubin-style circle fit (algebraic)
- Fits a circle in the bending plane (y–z) to averaged station positions.
- Momentum estimate: p = 0.3 |B| R (GeV/c), reading |B| from the field (kG) with fallbacks.
- Uncertainty p_err from radial residuals (heuristic); charge from sagitta sign × sign(Bx).
- Provides px/py/pz by projecting the local tangent (px close to zero by construction), plus chi2/ndf and a heuristic p-value.

## Output Data (TTree "Hits")

Per event, the following branches are written. Units are noted per branch family.

Raw per-hit (units from simulation):
- `eventID` (int), `trackID` (vector<int>), `pdg` (vector<int>)
- `stationID` (vector<int>), `layerID` (vector<int>)
- `x`, `y`, `z` (vector<double>) — hit positions (mm)
- `px`, `py`, `pz` (vector<double>) — hit momenta (MeV/c)

Muon reconstruction (per reconstructed muon track; individual vectors are aligned by index):
- `muon_trackID` (vector<int>)
- `muon_pdg` (vector<int>)
- `muon_truth_px`, `muon_truth_py`, `muon_truth_pz` (vector<double>) — truth momentum components (GeV/c)
- `muon_truth_p` (vector<double>) — truth momentum magnitude (GeV/c)
- `muon_reco_px`, `muon_reco_py`, `muon_reco_pz` (vector<double>) — legacy reco (mirrors GenFit) (GeV/c)
- `muon_reco_p` (vector<double>) — legacy reco p (GeV/c)
- `muon_nhits` (vector<int>) — number of averaged positions (stations) used by this track
- `muon_nstations` (vector<int>) — distinct stations for this track
- `muon_chi2`, `muon_ndf` (vector<double>) — legacy chi2/ndf (from GenFit)
- `muon_fit_success` (vector<bool>) — legacy fit success flag (from GenFit)

GenFit results (per track):
- `genfit_px`, `genfit_py`, `genfit_pz`, `genfit_p` (vector<double>) — GeV/c
- `genfit_chi2`, `genfit_ndf`, `genfit_pval` (vector<double>)
- `genfit_fit_success` (vector<bool>)

Taubin results (per track):
- `taubin_px`, `taubin_py`, `taubin_pz`, `taubin_p` (vector<double>) — GeV/c
- `taubin_p_err` (vector<double>) — GeV/c
- `taubin_chi2`, `taubin_ndf`, `taubin_pval` (vector<double>)
- `taubin_charge` (vector<int>) — inferred charge sign
- `taubin_fit_success` (vector<bool>)

Important: We fixed a branch-filling bug so `muon_truth_p` and `muon_reco_p` are cleared and filled every event. Earlier, the plotting loop that iterated over `muon_truth_p` would see zero-length vectors for early events despite fits being successful, resulting in “Tracks used: 0”. That is now resolved.

## Plotting and Resolution Table

The macro `Display/plot_muon_genfit_taubin.C` reads `scifi_hits_all.root` and produces:
- `build/muon_p_truth_vs_reco_genfit_taubin.(png|pdf)` — truth vs reco p for both methods
- `build/muon_p_resolution_genfit_taubin.(png|pdf)` — relative p resolution distributions
- `build/muon_components_resolution_genfit_taubin.(png|pdf)` — px/py/pz residuals
- `build/muon_resolution_table.csv` — CSV with entries, mean bias, and sigma for each metric and method

Run from the `Display/` directory:

```bash
root -l -b -q plot_muon_genfit_taubin.C
```

The macro auto-falls back to `scifi_hits_all.root` in the current dir if `build/scifi_hits_all.root` is absent.

Example after the branch-filling fix (quick run):

```
Tracks used: GenFit=4, Taubin=4
```

Counts will increase when running over all events and successful fits.

## Diagnostics and Extras

- ΣB⊥·dl diagnostic (optional): integrates the perpendicular field along the averaged hit polyline, estimates a deflection angle and a quick momentum estimate p_est ≈ 0.3·(ΣB⊥·dl)/θ. Enabled with `MS_DIAG_BDL=1`.
- Verbosity: set `MS_VERBOSE=1` to see additional prints for field and fitting.
- Geometry: `TGeoManager::Import("../../GeomGDML/geometry.gdml")` is used; ensure the relative path is valid from your run location.

## Known Issues and Edge Cases

- macOS runtime: running the executable outside `Display/build` may fail due to `@rpath` dylib lookup for Boost. Prefer running from the build dir.
- Tracks with < 3 stations are skipped (no fit attempted); these appear with default values and `fit_success=false`.
- Early files or special events may have sparse stations; expect fewer successful fits there.

## What Changed (Highlights)

- Dual-fit pipeline: always running both Taubin and GenFit (when feasible), recording both to the output tree.
- New branches added for explicit per-method results (`genfit_*`, `taubin_*`). Legacy `muon_*` mirrors GenFit for backward compatibility.
- Per-hit smearing and station averaging were implemented and logged.
- Geometry-aware magnetic field with per-event gating and sanity prints.
- New plotting macro and a CSV resolution export.
- Fixed vector alignment/branch filling for `muon_truth_p` and `muon_reco_p`.

## Suggested Next Steps

- Add unit tests or quick ROOT-based checks that validate branch sizes align across truth and reco per event.
- Improve Taubin momentum uncertainty estimation and charge inference robustness.
- Consider matching truth ↔ reco by `muon_trackID` in plotting to be resilient to any future vector alignment hiccups.
- Package a small script to run full reconstruction and plots end-to-end.

---

If anything needs refinement (thresholds, selection definitions, CSV schema), we can iterate quickly now that the full pipeline is wired and producing outputs.
