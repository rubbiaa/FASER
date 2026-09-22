# Simulate → reconstruct → compare regression tests

A second kind of test, alongside the `Tests/` gtests (`ctest`): those are
fast, pure-C++ unit tests of individual formulas/geometry (the muon field
model, the magnet-slit probe) and don't touch Geant4 at runtime. This suite
instead runs the *real* pipeline end to end — `faserps` (Geant4 simulation)
then `batchreco.exe` (reconstruction) — for a handful of event types, and
checks the output against a stored "golden" baseline. It's meant to catch
the class of bug this whole `run_faserps.py`/`run_batchreco.py`/`fb` effort
has been about all along: a code change that silently alters simulation or
reconstruction behavior, with no error, no crash, just a quietly different
answer.

**Status: designed, not yet run.** Everything below is grounded in reading
the actual source (tree/branch names, output file naming, the exact code
path each `run_faserps.py` flag takes) — not guessed — but I can't compile
or run Geant4/ROOT/PyROOT from where I work, so the two scripts below need
a first real run on your Mac before they're trustworthy. Treat the first
`--record` run as the actual validation step, not this write-up.

## Why simulate+reconstruct regression tests, on top of the gtests

The gtests check that two *formulas* agree (Geant4's field vs GenFit's
field) or that a geometry probe returns the right numbers for a shape built
by hand. They can't catch a bug in, say, `PrimaryGeneratorAction`'s event
generation, or a change to `BatchReco.cc`'s reconstruction logic, or a
regression introduced by touching `TPOEvent::kinematics_event()` — those
only show up once you actually run the pipeline on real events and look at
what comes out the other end.

## Determinism this relies on

`faserps.cc` hardcodes its random seed (`G4Random::setTheSeed(123456789)`,
`CLHEP::MTwistEngine`) — it's not configurable via macro or CLI. Given the
same code, the same input, the same `--n-events`, and single-threaded mode
(`/run/numberOfThreads 1`, `run_faserps.py`'s default), a run should be
**exactly** reproducible. That's what makes a golden-reference comparison
meaningful here rather than just noisy: a real code-behavior change should
show up as an exact mismatch, not something you have to squint at through a
statistical tolerance. If FASER ever exposes the seed as a CLI option,
that's a reason to revisit this (the regression suite must always pin it).
Multi-threaded mode (`--n-threads > 1`) is deliberately not covered here —
Geant4 MT's per-event RNG-stream determinism under multiple threads hasn't
been verified for this codebase and shouldn't be assumed.

## The pipeline, and what it actually writes

- `run_faserps.py` → `faserps`: one input mode (the default, reading
  neutrino-interaction events from a `TPOEvent`-format ROOT file) or two
  generator-only modes (`--muons`, `--muondis`, both `wantMuonBackground`).
  Writes **one ROOT file per event** — `FASERG4-Tcalevent_<run>_<event>.root`
  under `$FASERDATA/faserG4/` — each with a `calEvent` tree of exactly one
  entry, branch `event` holding the truth-level `TPOEvent` (see
  `CoreUtils/TcalEvent.cc`).
- `run_batchreco.py` → `batchreco.exe`: reads that range of per-event
  files back in (optionally filtered by `--mask`) and writes **one file per
  invocation** — `Batch-TPORecevent_<run>_<min>_<max>[_<mask>].root` under
  `$FASERDATA/batch/` — with a `RecoEvent` tree, branch `TPORecoEvent`, one
  entry per reconstructed event (see `Batch/BatchReco.cc`).
- Both `TPOEvent` and `TPORecoEvent` (and everything they reference —
  `TPORec`, `DigitizedTrack`, `PO`, ...) are covered by one ROOT
  dictionary, `CoreUtils/libCoreUtilsDict.so`, built automatically by the
  normal CMake build (see `CoreUtils/CMakeLists.txt`). `CoreUtils/dumpReco.C`
  is the existing (ROOT-macro) example of loading it and reading a
  `RecoEvent` tree — the summarizer script below is a PyROOT translation of
  exactly that pattern, not a new one.

### `run_number` per case — read from source, not guessed

- Default (file-input) mode: `run_number` comes straight from the input
  file's own stored `TPOEvent.run_number` (`PrimaryGeneratorAction.cc`
  reads each input entry directly into `fTPOEvent`). For the default
  sample (`$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root`), this is expected to be
  `10000` — matching the filename and the run number `run_batchreco.py`'s
  own docstring already uses as its basic example (`--run 10000`). **Not
  independently verified against the file's actual contents** (I can't
  open a ROOT file from here) — the first `--record` run will either
  confirm this or fail loudly with "0 events matched", which is the
  signal to fix this constant.

  (Historical note: before commit `32d4358`, `BatchReco.cc`'s per-event retry
  loop couldn't tell "this event isn't of the requested mask" apart from
  "truth file not written yet", and aborted the whole run the first time a
  masked truth file was missing -- typically at event 0, with `Total elapsed
  time: 0 ms`. That's now fixed: a missing masked truth file is treated as
  "not this event's type, move on", so a masked reco run correctly scans every
  event and only reports a real "0 events matched" when the sample genuinely
  has none of that type.)
- `--muons` and `--muondis`: both take the `wantMuonBackground` branch in
  `PrimaryGeneratorAction.cc`, which hardcodes `run_number = 999` — **the
  same number for both.** They must not share `$FASERDATA`, or one
  overwrites the other's per-event files. The orchestrator below gives
  every case its own `$FASERDATA` under `Tests/regression/work/<case>/`
  for exactly this reason, not just for tidiness.
- `nueCC`/`numuCC`/`nutauCC`/`nuNC` are not separate simulation runs at
  all — `TPOEvent::EncodeEventMask`/`BatchReco.cc` apply the mask at
  *reconstruction* time, filtering which of the same simulated run's
  events get reconstructed and written out. So the suite runs `faserps`
  **once** (the default neutrino case, run 10000) and reconstructs it
  **four times**, once per mask. A real risk worth flagging: with only a
  handful of simulated events, a given mask could easily match zero of
  them, depending on the input sample's actual composition (which I
  can't inspect). If a `--record` run comes back with 0 reconstructed
  events for some mask, that's the sample being too small/homogeneous,
  not a bug — bump `--n-events` for the neutrino case, or point at an
  input sample you know mixes interaction types.

## Layout

```
run_regression_tests.py          # orchestrator (repo root, alongside run_faserps.py/run_batchreco.py)
Tests/regression/
  summarize_output.py            # PyROOT: TPOEvent/TPORecoEvent -> a small JSON summary dict
  golden/                        # committed baselines, one JSON file per case
    neutrino_nueCC.json
    neutrino_numuCC.json
    neutrino_nutauCC.json
    neutrino_nuNC.json
    muons.json
    muondis.json
  work/                          # gitignored scratch: isolated $FASERDATA per case, logs
```

## Golden JSON schema

One file per case. `meta` records how it was produced (for humans reading
a diff, and so `--record` can warn if you're recording against a different
`--n-events` than the file expects); `truth`/`reco` are the actual
comparison targets — small dicts of aggregate numbers, not full per-event
dumps, so a diff is readable in a PR.

```json
{
  "meta": {
    "case": "neutrino_numuCC",
    "faserps_args": ["--n-events", "100"],
    "batchreco_args": ["--mask", "numuCC"],
    "run_number": 10000,
    "n_events_simulated": 100
  },
  "truth": {
    "n_events": 100,
    "mean_Evis": 12.34,
    "mean_n_particles": 7.2,
    "mean_nuE": 45.6,
    "mean_Q2": 3.1,
    "mean_xBj": 0.21
  },
  "reco": {
    "n_events_reconstructed": 41,
    "mean_n_PORecs": 3.4,
    "mean_n_TKTracks": 2.1,
    "mean_n_TKVertices": 1.0,
    "mean_n_MuTracks": 0.3,
    "mean_total_Evis_reco": 11.8,
    "mean_total_Ecompensated": 12.9
  }
}
```

`muons`/`muondis` omit the truth DIS-kinematics fields where they're not
meaningful (`nuE`/`Q2`/`xBj` are zero/unset for a muon primary that never
went through `TPOEvent::kinematics_event()`'s DIS branch) — the summarizer
only includes a field if it's actually applicable to that case, rather than
padding with meaningless zeros that would falsely "pass" a comparison.

## Comparison

Since the pipeline should be exactly reproducible (see "Determinism"
above), the default tolerance is deliberately tight — a relative
difference of `1e-9` for floating-point aggregates (room for harmless
compiler/platform floating-point differences, nothing else), and an exact
match for counts. Any larger difference is reported as a regression, with
old vs. new values printed side by side. `--record` overwrites the golden
file for the case(s) given instead of comparing — use it deliberately, the
same way you'd review any other change to committed test expectations, not
as a way to silence a failing comparison.

## Running it

```bash
# Build first, as always:
fb   # or: cmake --build build -j

# First time (or after a deliberate behavior change): record the baseline
python3 run_regression_tests.py --record

# Every other time: compare against the committed baseline
python3 run_regression_tests.py

# Just one case, e.g. while iterating on MuonDIS:
python3 run_regression_tests.py --case muondis
python3 run_regression_tests.py --case muondis --record
```

## Open items (need your input, not something to silently decide)

1. **CI needs an input sample it doesn't have -- mechanism now exists,
   CI itself not wired up yet.** The default neutrino case's input file,
   `$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root` (4.6 MB,
   moved out of `FASERG4/` into its own `GENIE/` subdirectory of the data
   dir since it's an input sample, not source code), is still
   `*.root`-gitignored (via the existing blanket `/data/` rule) and was
   never committed. It's now fetched from a public CERNBox link
   automatically instead: `fetch_data.py`'s `REMOTE_FILES` manifest (a
   public link, not a personal EOS token -- see that file's own docstring
   for why, and the commit that introduced it) has an entry for it, with a
   known sha256 checked after every download so a bad/expired link fails
   loudly instead of writing garbage into `data/GENIE/`. `run_faserps.py`
   calls this automatically when its default `--input-file` is missing, so
   a fresh `git clone` + `python3 run_faserps.py` now just works without a
   manual fetch step. **Still unverified**: cernbox.cern.ch wasn't
   reachable from either sandbox this was written in (proxy
   allowlist, not a CERNBox problem -- confirmed github.com worked fine
   from the same shells), so the actual download has only been logic-tested
   against a local HTTP server standing in for CERNBox, not the real URL --
   see `fetch_data.py`'s own docstring. Still open: wiring an explicit
   `python3 fetch_data.py` (or just letting `run_faserps.py`'s
   auto-fetch handle it) into `.github/workflows/build.yml` -- not done
   yet, and blocked on item 2 below anyway (CI doesn't run anything past
   Configure/Build yet). `--muons`/`--muondis` never needed this (they
   generate primaries directly).
2. **CI doesn't run `ctest` at all yet**, let alone this suite --
   `.github/workflows/build.yml` only configures and builds. Wiring in
   even the cheap gtests is a separate, smaller first step worth doing
   before adding a much heavier simulate+reconstruct job.
3. **Event count vs. mask coverage** (see above) -- may need tuning once
   you see real numbers from a `--record` run.

## Status of this file's own recommendations

- **Committing:** not yet -- you asked to validate first. Run `fb` then
  `python3 run_regression_tests.py --record` on this Mac; once that
  succeeds (or once you've fixed whatever it finds broken), the three new
  files (`run_regression_tests.py`, `Tests/regression/summarize_output.py`,
  this doc, plus `Tests/regression/golden/*.json` once recorded) are ready
  to `git add`/commit.
- **Branch:** directly on `main`, same as the docs-move commit -- no
  separate branch needed for this.
