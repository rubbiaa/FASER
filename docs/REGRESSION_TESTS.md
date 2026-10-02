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

This is a different knob from `batchreco.exe`'s own `-mt` flag
(`TPORecoEvent::multiThread`, parallelizing `Reconstruct3DPS_2`'s
per-module voxel reconstruction across one `std::thread` per detector
module) — `run_regression_tests.py` runs with `-mt` **on by default**
(`--multi-thread`/`--no-multi-thread`). The existing `golden/*.json` were
recorded single-threaded, and the one comparison run with `-mt` on so far
passed against them — evidence, not proof, that the per-module work is
independent and race-free. A future run that fails with the default but
passes again under `--no-multi-thread` would point at a real race in
`Reconstruct3DPS_2`/`reconstruct3DPS_module`, not noise.

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
  reads each input entry directly into `fTPOEvent`). The neutrino case's
  sample is resolved by `run_regression_tests.py`'s own
  `resolve_neutrino_input_file()`, which reuses `run_faserps.py`'s CVGENIE
  discovery helper (`resolve_cvgenie_po_file()`) to find
  `$FASERDATA/CVGENIE/Run10000/FASERMC-PO-Run10000-..._3DCAL.root` —
  before `run_group()` isolates `$FASERDATA` to that case's own scratch
  directory — and passes it to `run_faserps.py` explicitly via
  `--input-file`, since there's no implicit "default sample" any more
  (see `run_faserps.py`'s CVGENIE auto-discovery/`--input-file` in
  `docs/HOWTO.md`). `run_number: 10000`/`cvgenie_detector: "3DCAL"` are
  the `SIMULATION_GROUPS` entry's own fields, matching the run number
  `run_batchreco.py`'s own docstring already uses as its basic example
  (`--run 10000`).
- `--muons` and `--muondis`: both take the `wantMuonBackground` branch in
  `PrimaryGeneratorAction.cc`, which hardcodes `run_number = 999` — **the
  same number for both.** They must not share `$FASERDATA`, or one
  overwrites the other's per-event files. The orchestrator below gives
  every case its own `$FASERDATA` under `Tests/regression/work/<case>/`
  for exactly this reason, not just for tidiness.
- `nueCC`/`numuCC`/`nutauCC`/`nuNC` are not separate simulation **or**
  **reconstruction** runs: the suite runs `faserps` **once** (the default
  neutrino case, run 10000) and `batchreco.exe` **once**, fully unmasked,
  over the whole unbiased sample. The four golden cases are a *post-hoc
  split in Python*, not four separate `--mask` invocations.

  This replaces an earlier, incorrect design that reconstructed the same
  run four times with `--mask nueCC`/`numuCC`/`nutauCC`/`nuNC` and relied
  on `Batch/BatchReco.cc`'s per-event retry loop to skip indices that
  didn't match. That design assumed truth files get tagged with a mask at
  *simulation* time (`TcalEvent`'s constructor appends `_<mask>` to the
  filename `if(event_mask>0)`), but in this workflow `event_mask` is never
  anything but 0: `TPOEvent::SetEventMask()` is only ever called from
  `ConvertFASERMC.cc`'s conversion path, not from
  `FASERG4`/`FASERCalProtoG4`'s `ParticleManager.cc`, which is what
  actually constructs the `TcalEvent` that writes
  `FASERG4-Tcalevent_*.root` for this (CVGENIE-based) neutrino case. So
  every truth file in this sample is unmasked, every masked
  `Load_event()` lookup misses for *every* event index, not just "most"
  of them, and a `--mask`-filtered `batchreco.exe` run here always
  produces an empty reco file — regardless of how the input sample is
  actually composed. (The historical note this replaced, about commit
  `32d4358` fixing `BatchReco.cc`'s retry loop to tell "wrong mask" apart
  from "not written yet", is still accurate as far as it goes — it's a
  real fix, just for a scenario this suite's own neutrino case never
  actually exercises.)

  Instead, `run_regression_tests.py`'s `SIMULATION_GROUPS` entry for
  "neutrino" has a `split_by_reaction` dict (golden name — reaction
  string) instead of a `reco_cases` list, and `run_group()` calls
  `summarize_output.py --split-by-reaction` once on the single unmasked
  truth+reco output. That script buckets every truth event (by its own
  `TPOEvent::reaction_desc()`) and every reconstructed event (by its own
  truth `TPOEvent`, read back via `TPORecoEvent::GetPOEvent()` —
  `Batch/BatchReco.cc` constructs every `TPORecoEvent` with the exact
  `TPOEvent*` for that event, and that pointer is a persisted member of
  `TPORecoEvent`, so a reco entry carries its own classification with it,
  no separate truth-file lookup needed to classify it) into
  `{reaction: {"truth": {...}, "reco": {...}}}`, and `run_group()`
  maps each of the four golden names to its reaction's bucket (an empty
  `{"n_events": 0}`/`{"n_events_reconstructed": 0}` if the sample
  happens to have none of that flavor — see the risk noted below). A
  real risk worth flagging: with only a handful of simulated events, a
  given reaction could easily account for zero of them, depending on the
  input sample's actual composition. If a `--record` run comes back with
  0 events for some case, that's the sample being too small/homogeneous,
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
  results/                       # gitignored scratch: last run's numbers + comparison.csv
```

**`golden/*.json` is committed directly in git, not on CERNBox via
`fetch_data.py`.** This has come up more than once, so it's worth settling
here: `fetch_data.py` exists for *inputs* too big/binary to commit (and
`.gitignore`d for exactly that reason) -- the golden files are the
opposite on both counts, six files at ~4KB each, plain diffable text,
and `.gitignore` deliberately carves `golden/` out as *not* ignored.
More importantly, a golden file isn't an input, it's a baseline *output*
of a specific commit of the reconstruction code -- keeping it in git
means a PR that changes `BatchReco.cc`/`TPORecoEvent.cc` shows the
resulting baseline diff right alongside the code change, for a reviewer
to question. Hosting it externally and fetching it would break that
link: code could drift out of sync with the baseline with nothing in
the repo's history to show it. (The one real CERNBox-fetch candidate
left is the neutrino case's raw CVGENIE input sample -- see "Open
items" below -- which is genuinely large, binary, and has no
code-review value in its bytes; that's a different problem.)

Keeping it in git only helps if a PR touching `golden/` actually gets a
human look before merging -- CI can't verify a baseline change is correct
(that needs the real ROOT/Geant4 pipeline, see "Open items" below), so
there are two GitHub-native guards instead of one:

- **`.github/CODEOWNERS`** names a required reviewer for
  `Tests/regression/golden/` (and the code that produces it,
  `run_regression_tests.py`/`summarize_output.py`) -- but it only takes
  effect once **"Require review from Code Owners"** (and "Require a pull
  request before merging") is turned on under this repo's Settings ->
  Branches -> branch protection rule for `main`. That toggle has to be
  set in GitHub's own settings; nothing in the repo can turn it on for
  you.
- **`.github/workflows/golden-baseline-guard.yml`** runs on every PR that
  touches `golden/**` and posts a `::warning::` annotation per changed
  file -- purely a visibility aid (it never fails the build; CI has no
  way to tell a deliberate `--record` from an accidental one), so a
  baseline change can't quietly slip past in a large diff unnoticed.

Both of those only cover GitHub -- `run_regression_tests.py --record`
itself also refuses to overwrite a golden file that already exists unless
you pass `--force` too (a brand-new case's first `--record` needs no
`--force`, since there's nothing to overwrite yet). So even running
`--record` locally without thinking can't silently clobber a committed
baseline; `--force` is there for exactly the deliberate, reviewed update
this whole section is about -- never to make a failing comparison go
away.

## Golden JSON schema

One file per case. `meta` records how it was produced (for humans reading
a diff, and so `--record` can warn if you're recording against a different
`--n-events` than the file expects); `truth`/`reco` are the actual
comparison targets — small dicts of aggregate numbers (`mean_*` and,
alongside every one of them, `rms_*` -- the spread around that mean, ROOT's
`TH1::GetRMS()` convention), not full per-event dumps, so a diff is
readable in a PR. This is the real, literal content of the committed
`golden/neutrino_numuCC.json` (reordered meta/truth/reco for readability;
the file itself is alphabetical within each section):

```json
{
  "meta": {
    "case": "neutrino_numuCC",
    "reaction": "numuCC",
    "faserps_args": [],
    "run_number": 10000,
    "n_events_simulated": 100
  },
  "truth": {
    "n_events": 25,
    "mean_Evis": 595.731565396213,
    "rms_Evis": 401.1094242759732,
    "mean_n_particles": 17.88,
    "rms_n_particles": 6.519631891449087,
    "mean_nuE": 257.89267689473013,
    "rms_nuE": 206.50216825257777,
    "mean_Q2": 77.26041397106434,
    "rms_Q2": 58.725522257556854,
    "mean_xBj": 0.2018162065692521,
    "rms_xBj": 0.14032735412008787,
    "mean_yInel": 0.5048221131155411,
    "rms_yInel": 0.2758859573043792
  },
  "reco": {
    "n_events_reconstructed": 25,
    "mean_n_PORecs": 14.76,
    "rms_n_PORecs": 4.9094195176212025,
    "mean_n_TKTracks": 0.0,
    "rms_n_TKTracks": 0.0,
    "mean_n_TKVertices": 0.0,
    "rms_n_TKVertices": 0.0,
    "mean_n_MuTracks": 0.72,
    "rms_n_MuTracks": 0.6013318551349164,
    "mean_total_Evis_reco": 72.62853511678694,
    "rms_total_Evis_reco": 66.27062278442045,
    "mean_total_Ecompensated": 72.62853511678694,
    "rms_total_Ecompensated": 66.27062278442045
  }
}
```

For `neutrino_*`, `truth.n_events`/`reco.n_events_reconstructed` are how
many of the sample's 100 events actually turned out to be that reaction
— not 100, since the sample is unbiased and mixes flavors (see
`split_by_reaction` above; `25` here means 25 of the 100 simulated events
were `numuCC`). For `muons`/`muondis` (which don't split by reaction,
every event already being the same type) those counts do equal
`n_events_simulated` once batchreco successfully reconstructs every one.

`muons`/`muondis` omit the truth DIS-kinematics fields where they're not
meaningful (`nuE`/`Q2`/`xBj`/`yInel` are zero/unset for a muon primary that
never went through `TPOEvent::kinematics_event()`'s DIS branch) — the
summarizer only includes a field (and its `rms_` partner) if it's actually
applicable to that case, rather than padding with meaningless zeros that
would falsely "pass" a comparison. See `golden/muondis.json` in this repo
for a real example of that narrower shape.

## Comparison

Since the pipeline should be exactly reproducible (see "Determinism"
above), the default tolerance is deliberately tight — a relative
difference of `1e-9` for floating-point aggregates (room for harmless
compiler/platform floating-point differences, nothing else), and an exact
match for counts. Any larger difference is reported as a regression, with
old vs. new values printed side by side. `--record` writes the golden
file for the case(s) given instead of comparing — use it deliberately, the
same way you'd review any other change to committed test expectations, not
as a way to silence a failing comparison. It refuses to overwrite a golden
file that already exists unless `--force` is also given; see the
baseline-protection note above.

## Running it

The orchestrator script itself (`run_regression_tests.py`) doesn't need
PyROOT -- only the `summarize_output.py` subprocess it spawns does, and
that has to be run with a Python whose major.minor version matches the
one your ROOT build was linked against (`thisroot.sh` doesn't repoint
`python3` itself, so an active conda/venv on top easily shadows it with a
mismatched version -- PyROOT's `import ROOT` fails loudly, naming both
versions, if so). `setup.sh` exports `$FASER_PYTHON` for the sites where
this matters (currently André's Mac, pointing at the Homebrew Python 3.14
ROOT was built against there), and `--python` defaults to it when set, so
sourcing `setup.sh` is normally enough -- no need to pass `--python` by
hand. Adding a new site whose ROOT needs a non-default interpreter: export
`FASER_PYTHON` in that site's `setup.sh` branch, the same way
`GEANT4_INSTALL`/`PYTHIA8` are already done there.

```bash
# Build first, as always:
fb   # or: cmake --build build -j

# First time for a case with no golden/<case>.json yet: record the baseline
python3 run_regression_tests.py --record

# Every other time: compare against the committed baseline
python3 run_regression_tests.py

# After a deliberate, reviewed behavior change: re-record (needs --force,
# since a golden file already exists -- see the baseline-protection note
# above; --record without --force refuses rather than silently overwriting)
python3 run_regression_tests.py --record --force

# Just one case, e.g. while iterating on MuonDIS:
python3 run_regression_tests.py --case muondis
python3 run_regression_tests.py --case muondis --record --force

# $FASER_PYTHON not set (not sourcing setup.sh) or you need a different
# interpreter one-off: point --python at it explicitly --
python3 run_regression_tests.py --python /opt/homebrew/bin/python3.14

# Iterating on reconstruction-only code (BatchReco.cc, TPORecoEvent, ...)?
# Re-simulating every time is pure overhead once the truth sample exists --
# run once normally, then skip straight to batchreco.exe on later runs:
python3 run_regression_tests.py --case muondis           # first run: simulates + reconstructs
python3 run_regression_tests.py --case muondis --skip-faserps   # later runs: reconstructs only

# Iterating on summarize_output.py or the golden-comparison logic itself?
# Re-running batchreco.exe every time is pure overhead once the reco file
# exists -- run once normally, then skip straight to summarize on later runs
# (combine with --skip-faserps to skip both and go straight to summarize):
python3 run_regression_tests.py --case muondis --skip-reco      # later runs: summarizes only
python3 run_regression_tests.py --case muondis --skip-faserps --skip-reco  # summarize-only, no sim or reco
```

`--skip-faserps` reuses whatever truth sample is already sitting in that
case's `work/<key>/data/faserG4` -- it's on you to know the sample is
still valid for what you're testing (it doesn't hash `faserps_args` or
otherwise detect staleness); it fails fast with a clear error if no truth
files are there yet. Note this is a *group*-level reuse: the "neutrino"
group's four golden cases (`neutrino_nueCC`/`numuCC`/`nutauCC`/`nuNC`)
share one simulated sample, so `--case neutrino_nueCC --skip-faserps`
reuses the same sample `--case neutrino_numuCC` would have produced.

`--skip-reco` is the same idea, one step further down the pipeline: it
reuses whatever reco file is already sitting in that case's
`work/<key>/data/batch` -- same caveat (it's on you to know it's still
valid), same fail-fast behavior if the file isn't there yet, and same
group-level reuse for the "neutrino" case's four golden names. It's
independent of `--skip-faserps`; pass both together to jump straight to
the summarize step.

Every run -- `--record` or a plain comparison -- also writes each case's
freshly-computed summary under `Tests/regression/results/<case>.json`:
the actual numbers (`mean_*`/`rms_*`/etc.), not just the PASS/FAIL line
printed to stdout. Unlike `golden/`, `results/` isn't a committed
baseline (it's gitignored) -- it's just "what did the last run actually
compute", there so you don't have to re-run with `--record` (clobbering
the real baseline) just to look at the numbers.

A plain comparison run (not `--record`) additionally writes
`Tests/regression/results/comparison.csv`, one row per field per case:
`case,section,field,golden,current,abs_diff,rel_diff,status`, for every
field -- not only the FAILs the terminal prints. `--rel-tol` (default
`1e-9`) is what `status` is based on, and it's worth being explicit about
what that number means: it's near machine precision, not a tolerance for
expected statistical fluctuation. This suite's whole premise (see
"Determinism this relies on" below) is that a run is **exactly**
reproducible, so `PASS` means "matches to ~15 significant figures", and
any real `FAIL` is a signal that something actually changed -- not noise
to average away or a tolerance to loosen.

`faserps`/`batchreco.exe`'s own stdout+stderr (Geant4 init, MDT geometry
scans, per-event progress, ...) no longer streams to the terminal --
it's redirected to `Tests/regression/work/<key>/logs/{faserps,batchreco*
}.log`, with just a one-line pointer printed in its place. A failure
still prints the last ~40 lines of the relevant log immediately, so
there's no need to go open the file just to see what broke.

## Open items (need your input, not something to silently decide)

1. **CI needs an input sample it doesn't have, and there's no auto-fetch
   any more.** The neutrino case needs a real, already-converted CVGENIE
   sample on disk at `$FASERDATA/CVGENIE/Run10000/FASERMC-PO-Run10000-..._3DCAL.root`
   (see `docs/HOWTO.md`'s `run_convertgenie.py` section for how that's
   produced from raw GENIE output). The old CERNBox auto-fetch path
   described here previously (a `fetch_data.py` `REMOTE_FILES` entry for
   the final converted PO file, called automatically by `run_faserps.py`
   when its old hardcoded default `--input-file` was missing) no longer
   exists — it was removed when `GENIE/` was reorganized into `CVGENIE/`
   and `run_faserps.py` switched to auto-discovering
   `$FASERDATA/CVGENIE/Run<run>/` directories instead of pointing at one
   fixed default path (see `run_faserps.py`'s own docstring/CLI help).
   `fetch_data.py`'s manifest now only covers the *raw* GENIE output
   (`GENIE/fasercal.Aki2024.v10.{charm,light}.0.gfaser.root`) and the
   script that produced it — not the converted PO file itself, since
   converting is a build-and-run step (`ConvertGENIE.exe` via
   `run_convertgenie.py`), not something to auto-fetch. So: on a machine
   that doesn't already have `$FASERDATA/CVGENIE/Run10000/.../..._3DCAL.root`
   on disk, `run_regression_tests.py`'s neutrino case fails fast with a
   clear "run `python3 run_convertgenie.py --run 10000` first" error
   (`resolve_neutrino_input_file()`) rather than silently fetching
   anything. Wiring this into CI needs either committing a small
   converted sample (still `*.root`-gitignored, same blanket `/data/`
   rule as before) or running `fetch_data.py` + `run_convertgenie.py` as a
   CI step — not done yet, and blocked on item 2 below anyway (CI
   doesn't run anything past Configure/Build yet). `--muons`/`--muondis`
   never needed this (they generate primaries directly).
2. **CI doesn't run `ctest` at all yet**, let alone this suite --
   `.github/workflows/build.yml` only configures and builds. Wiring in
   even the cheap gtests is a separate, smaller first step worth doing
   before adding a much heavier simulate+reconstruct job.
3. **Event count vs. reaction-type coverage** (see "run_number per case"
   above) -- may need tuning once you see real numbers from a `--record`
   run.

## Status of this file's own recommendations

- **Committing:** not yet -- you asked to validate first. Run `fb` then
  `python3 run_regression_tests.py --record` on this Mac; once that
  succeeds (or once you've fixed whatever it finds broken), the three new
  files (`run_regression_tests.py`, `Tests/regression/summarize_output.py`,
  this doc, plus `Tests/regression/golden/*.json` once recorded) are ready
  to `git add`/commit.
- **Branch:** directly on `main`, same as the docs-move commit -- no
  separate branch needed for this.


## Example: what a typical run looks like

A plain `python3 run_regression_tests.py` (no flags -- compares every case
against the committed `golden/*.json`, `-mt` on by default) looks roughly
like this. The paths below are illustrative (yours will show your own
`$FASERDATA`/build dir and the real neutrino input file resolved by
`resolve_neutrino_input_file()`); the terminal output *format* and the
numbers in the file excerpts further down are real, from the baselines
actually committed in this repo.

```
[run_regression_tests] neutrino: faserps: /opt/homebrew/bin/python3.14 run_faserps.py --input-file /Users/rubbiaa/data/CVGENIE/Run10000/FASERMC-PO-Run10000-..._3DCAL.root --n-events 100 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build
[run_regression_tests] neutrino: faserps: output -> Tests/regression/work/neutrino/logs/faserps.log
[run_regression_tests] neutrino: batchreco: /opt/homebrew/bin/python3.14 run_batchreco.py --run 10000 --max-event 100 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build --multi-thread
[run_regression_tests] neutrino: batchreco: output -> Tests/regression/work/neutrino/logs/batchreco.log
[run_regression_tests] neutrino: summarize: /opt/homebrew/bin/python3.14 Tests/regression/summarize_output.py --truth-dir .../work/neutrino/data/faserG4 --run 10000 --n-events 100 --reco-file .../work/neutrino/data/batch/Batch-TPORecevent_10000_0_100.root --split-by-reaction
[run_regression_tests] neutrino_nueCC: PASS
[run_regression_tests] neutrino_numuCC: PASS
[run_regression_tests] neutrino_nutauCC: PASS
[run_regression_tests] neutrino_nuNC: PASS
[run_regression_tests] muons: faserps: /opt/homebrew/bin/python3.14 run_faserps.py --muons --n-events 20 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build
[run_regression_tests] muons: faserps: output -> Tests/regression/work/muons/logs/faserps.log
[run_regression_tests] muons: batchreco: /opt/homebrew/bin/python3.14 run_batchreco.py --run 999 --max-event 20 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build --multi-thread
[run_regression_tests] muons: batchreco: output -> Tests/regression/work/muons/logs/batchreco_muons.log
[run_regression_tests] muons: summarize: /opt/homebrew/bin/python3.14 Tests/regression/summarize_output.py --truth-dir .../work/muons/data/faserG4 --run 999 --n-events 20 --reco-file .../work/muons/data/batch/Batch-TPORecevent_999_0_20.root
[run_regression_tests] muons: PASS
[run_regression_tests] muondis: faserps: /opt/homebrew/bin/python3.14 run_faserps.py --muondis --n-events 20 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build
[run_regression_tests] muondis: faserps: output -> Tests/regression/work/muondis/logs/faserps.log
[run_regression_tests] muondis: batchreco: /opt/homebrew/bin/python3.14 run_batchreco.py --run 999 --max-event 20 --build-dir /Users/rubbiaa/MACDEV/FASERV9/FASER/build --multi-thread
[run_regression_tests] muondis: batchreco: output -> Tests/regression/work/muondis/logs/batchreco_muondis.log
[run_regression_tests] muondis: summarize: /opt/homebrew/bin/python3.14 Tests/regression/summarize_output.py --truth-dir .../work/muondis/data/faserG4 --run 999 --n-events 20 --reco-file .../work/muondis/data/batch/Batch-TPORecevent_999_0_20.root
[run_regression_tests] muondis: PASS
[run_regression_tests] full per-case results (the actual computed numbers, not just PASS/FAIL) written under Tests/regression/results/
[run_regression_tests] field-by-field comparison against golden/ for every case (148 fields, 0 FAIL) written to Tests/regression/results/comparison.csv
```

Exit code `0` (every case `PASS`). The `summarize` steps print no "output ->"
line -- they're still `capture=True` (their stdout is the JSON `run_group()`
parses), not redirected to a log file like `faserps`/`batchreco` are.

### What a FAIL looks like

If a code change actually shifted a value -- say `muondis`'s
`mean_total_Ecompensated` came out `47.1` instead of the golden
`46.52888077839226` -- the terminal would show:

```
[run_regression_tests] muondis: FAIL
    reco.mean_total_Ecompensated: golden=46.52888077839226 current=47.1
```

and the exit code would be `1`. Only the FAILing field(s) print to the
terminal; `Tests/regression/results/comparison.csv` would still have a row
for every field, FAILs included:

```
muondis,reco,mean_total_Ecompensated,46.52888077839226,47.1,0.5711192216077379,0.012274510197824535,FAIL
```

(`rel_diff` here is ~1.2%, far above `--rel-tol`'s default `1e-9` -- see
"Determinism this relies on" above for why that's treated as a real
regression, not noise.)

### `results/<case>.json` and `comparison.csv`, for real

`Tests/regression/results/muondis.json` (written on every run, `--record`
or not) would currently read exactly like `golden/muondis.json` does,
since the last recorded run and the current golden agree on every field --
see `golden/muondis.json` in this repo for the full file. Its
`comparison.csv` rows (real output, from running `diff_rows()` against the
actual committed `golden/muondis.json`) look like this (truncated to a
few representative fields -- the real file has one row per field, 26 for
`muondis` alone, 148 across all six cases):

```
case,section,field,golden,current,abs_diff,rel_diff,status
muondis,truth,mean_Evis,540.231094535469,540.231094535469,0.0,0.0,PASS
muondis,truth,n_events,20,20,0,0.0,PASS
muondis,truth,rms_Evis,601.0602382579357,601.0602382579357,0.0,0.0,PASS
muondis,reco,mean_n_PORecs,1.0,1.0,0.0,0.0,PASS
muondis,reco,mean_n_TKTracks,0.0,0.0,0.0,,PASS
muondis,reco,mean_total_Ecompensated,46.52888077839226,46.52888077839226,0.0,0.0,PASS
muondis,reco,rms_total_Ecompensated,153.55849336513893,153.55849336513893,0.0,0.0,PASS
```

(`rel_diff` is blank when `golden` is `0` -- see `_compare_value()`'s
docstring: the check falls back to an absolute one against `--rel-tol`
instead, since a relative difference against zero is undefined.)
