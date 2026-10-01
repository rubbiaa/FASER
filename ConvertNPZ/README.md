# ConvertNPZ

ConvertNPZ turns FASER events that are already simulated and reconstructed
into NPZ files, one per event. It does not generate, simulate or reconstruct
anything: it reads two kinds of existing files.

| Input | Produced by | Provides |
| --- | --- | --- |
| RECO, `Batch-TPORecevent_<run>_<first>_<last>.root` | BatchReco | 3DCal voxels, event truth record, muon spectrometer fit |
| TCAL, `FASERG4-Tcalevent_<run>_<event>.root` | FASERG4 | 3DCal truth hits, rear ECAL and AHCAL deposits, MDT truth |

Events are matched by run and event number, so both inputs must come from the
same production. A RECO event without its TCAL file stops the conversion with
an error.

## Before you start: build FASER

ConvertNPZ is not standalone. It reads FASER's C++ classes through PyROOT with
the library that the FASER build creates in this checkout,
`CoreUtils/libCoreUtilsDict.so` (`.dylib` on macOS). You need:

- all of FASER's requirements, including ROOT with PyROOT; see
  [INSTALL.md](../docs/INSTALL.md), which also covers machines that
  `setup.sh` does not recognise;
- FASER built in this checkout;
- the Python that ROOT was built with, with `numpy` and `tqdm`.

From the top of the FASER checkout:

1. Set up FASER's environment (ROOT, Geant4, GenFit and `FASERDATA`):

   ```bash
   source setup.sh
   ```

2. Build FASER. The build also generates the library ConvertNPZ loads:

   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
   cmake --build build -j
   ```

   If the library is missing or out of date, for example after pulling
   changes to `CoreUtils`, regenerate only the library:

   ```bash
   ./ConvertNPZ/generateLib.sh
   ```

3. Install the Python packages:

   ```bash
   python3 -m pip install -r ConvertNPZ/requirements.txt
   ```

4. Check the setup:

   ```bash
   ls CoreUtils/libCoreUtilsDict.*
   python3 -c "import ROOT, numpy, tqdm"
   ```

In every new shell, run `source setup.sh` before using ConvertNPZ.

## Convert

Once the RECO and TCAL files exist, choose the sample type:

```bash
python3 ConvertNPZ/convert.py --mode neutrino
# Or, for single-particle data:
python3 ConvertNPZ/convert.py --mode single-particle
```

By default it uses the same directories as FASER's simulation and
reconstruction:

| Files | Default directory |
| --- | --- |
| RECO files | `$FASERDATA/batch` |
| TCAL files | `$FASERDATA/faserG4` |
| NPZ output | `$FASERDATA/npz` |
| Particle geometry | `$FASERDATA/GDML/FASERCAL_V10.gdml` |

`FASERDATA` defaults to `data/` inside the checkout, as in FASER. Each saved
file is named `run_<run>_event_<event>.npz`, with run and event numbers taken
from ROOT. Existing files are overwritten unless `--skip-existing` is supplied.

## Choosing the input files

Convert one RECO file and save at most five events:

```bash
python3 ConvertNPZ/convert.py --mode neutrino \
  --reco-file /path/to/Batch-TPORecevent_10000_0_100.root \
  --tcal-dir /path/to/tcal --output-dir /path/to/npz --max-events 5
```

Use `--reco-dir` instead of `--reco-file` for a directory. RECO and TCAL
directories can be flat or split into `chunk_*` subdirectories. Empty RECO
files produced by masked runs are skipped; missing or unreadable inputs report
an error.

The older production layout, with `FASERCALDATA_*` and `FASERCALRECODATA_*`
directories, is also supported:

```bash
python3 ConvertNPZ/convert.py --mode neutrino \
  --base-path /path/to/FASERDATA_v10 --version v10.0_10000 \
  --output-dir /path/to/npz
```

Single-particle mode also needs the GDML geometry exported by the simulation of
the input production: pass `--geometry /path/to/geometry.gdml` or set
`FASER_NPZ_GEOMETRY`. Optional `--energy 50` writes under a `50GeV/`
subdirectory.

## Interactions in the rear ECAL or AHCAL

Every NPZ contains the rear calorimeter deposits from TCAL:

| Field | Content |
| --- | --- |
| `ecal_hits` | (N, 4) float32: `ix`, `iy`, `layer`, energy in MeV, one row per ECAL cell |
| `ahcal_hits` | (N, 4) float32: the same for the AHCAL |

The energy is the simulated charged-particle deposit in the scintillator of each
cell; BatchReco copies it unchanged. FASERG4 does not record which particle made
an ECAL or AHCAL deposit, so `true_hits`, `reco_hits` and `seg_labels` cover the
3DCal only.

When the vertex is in the ECAL or AHCAL (ConvertGENIE run with `ECAL`, `AHCAL`
or `ALL`), disable the 3DCal cut:

```bash
python3 ConvertNPZ/convert.py --mode neutrino --min-non-ghost 0 \
  --reco-dir /path/to/reco --tcal-dir /path/to/tcal --output-dir /path/to/npz
```

These events leave few or no 3DCal voxels. In the v9 test samples, 63% of ECAL
and 99% of AHCAL events had none, and the default cut would have removed 98%
and 100% of them. Events without 3DCal voxels are saved with empty `reco_hits`
and `seg_labels`; `true_hits` keeps any simulated 3DCal deposits.

## Output

Neutrinos and single particles share the same conversion code and keep the
fields, array shapes and dtypes of the original production NPZ files:

| Mode | Default fields | Additional fields beyond the shared 36 |
| --- | --- | --- |
| `neutrino` | 41 | `muspec_true`, `muspec_reco_info`, `muspec_reco_extra`, `muspec_reco_diagnostics`, `muspec_reco_trackid` |
| `single-particle` | 37 | `primary_momentum` |

The calculations come from the production scripts by Dr. Saul Alonso-Monsalve
and Fabio Cufino. The original CAL loading and truth-hit conventions are
preserved: the neutrino mode excludes channel ID 0, while particle mode retains
its original handling. Neutrino CAL lookup accepts interaction suffixes such
as `_nueCC.root` and matches both run and event IDs.

## Options

- `--max-events`: limit saved events, after cuts.
- `--min-non-ghost`: defaults to 20 non-ghost PS voxels; use 0 to disable the
  cut, as needed for ECAL and AHCAL vertices.
- `--skip-existing`: skip existing output filenames.
- `--number` / `--chunks`: split the RECO file list; defaults are 0 / 1.
- `--include-views`: add the PS 2D views.

Use `--help` for all options.

To check that RECO files written by another FASER version can be read by this
build, run the preflight check on one file before a large conversion:

```bash
python3 ConvertNPZ/preflight_check_reco.py /path/to/Batch-TPORecevent_10000_0_100.root
```

`FASER_NPZ_LIBRARY` points ConvertNPZ to a library built from another FASER
version, for productions made with it.
