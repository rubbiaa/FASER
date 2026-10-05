# condor/ - run faserps and batchreco on the CERN batch farm

Submits chunked `run_faserps.py` / `run_batchreco.py` jobs to HTCondor on
lxplus, either directly on lxplus or from your Mac. It uses the repo as it
is: the lxplus checkout's own `setup.sh`, `build/` and `$FASERDATA`
layout. Nothing is copied to the worker nodes except one small wrapper
script - batch nodes see the checkout (`/afs` or `/eos`), `/cvmfs` and the
shared data area directly.

| File | Runs on | Purpose |
|---|---|---|
| `submit_from_mac.sh` | Mac | rsync `condor/` to lxplus, run `submit_jobs.sh` there over ssh; `status`, `sync`, `pull`, `ssh` helpers |
| `submit_jobs.sh` | lxplus | resolve the input, build the chunk list (`make_jobs_list.py`), `condor_submit` |
| `submit_faserps.sub`, `submit_batchreco.sub` | lxplus | HTCondor descriptions (per-run values come in via `-append`) |
| `faserps_chunk.sh`, `batchreco_chunk.sh` | batch node | one chunk: set up, run, verify, publish |
| `make_jobs_list.py` | lxplus | `(chunk, start, nevt)` list |

## One-time setup (on lxplus)

```bash
git clone https://github.com/rubbiaa/FASER.git      # or use your existing checkout
cd FASER && source setup.sh && fb                   # build faserps + batchreco.exe
```

On the Mac, tell the script where that checkout is (defaults: host `lxplus`,
checkout `/afs/cern.ch/work/r/rubbiaa/FASER` - adjust to yours):

```bash
cat > ~/.faser_condor.conf <<'EOF'
LXPLUS_HOST=lxplus
LXPLUS_FASER=/afs/cern.ch/work/r/rubbiaa/FASER
EOF
```

On lxplus, point the jobs at a shared data area on EOS (not in the AFS
checkout - event files are ~200 kB each):

```bash
mkdir -p /eos/user/r/rubbiaa/FASERDATA
echo 'FASER_SHARED_DATA=/eos/user/r/rubbiaa/FASERDATA' > ~/.faser_condor.conf
```

(or pass `--faserdata DIR` every time). That area is the "shared
`$FASERDATA`": inputs are read from it and results are published into it -
`CVGENIE/Run<N>/` (converted GENIE sample, produced by
`run_convertgenie.py` or copied up from your Mac's `data/CVGENIE`),
`faserG4/`, `GDML/`, `batch/`.

Tip: put `ControlMaster auto` / `ControlPersist 4h` for the lxplus host in
`~/.ssh/config` so the script's several ssh/rsync calls need one 2FA only.

## Usage

```bash
cd ~/MACDEV/FASERV9/FASER/condor

# 1. simulate: 350000 events in chunks of 5000 from the 3DCAL sample of CVGENIE run 10000
./submit_from_mac.sh faserps --run 10000 --detector 3DCAL --total-events 350000 --chunk-size 5000 --dry-run
./submit_from_mac.sh faserps --run 10000 --detector 3DCAL --total-events 350000 --chunk-size 5000

# other faserps flavours: any run_faserps.py flag goes after a literal --
./submit_from_mac.sh faserps --total-events 100000 -- --muons --muon-momentum-gev 250
./submit_from_mac.sh faserps --input-file /eos/.../FASERMC-PO-....root --total-events 50000 -- --tilt-deg -4.5

# 2. watch
./submit_from_mac.sh status

# 3. when faserps is done, reconstruct everything it produced for that run
./submit_from_mac.sh batchreco --run 10000 --chunk-size 2000
```

On lxplus itself, skip the Mac script: `cd FASER/condor && ./submit_jobs.sh ...`
takes the same arguments. Job lists end up in `condor/jobs/`, logs in
`condor/logs/` (on lxplus).

## What the jobs do

* **Job-local `$FASERDATA`.** Each job exports `FASERDATA=<scratch>/data`
  *before* sourcing the repo's `setup.sh` (which only defaults it if unset),
  runs the repo's own `run_faserps.py` / `run_batchreco.py`, checks the
  result, and only then copies it to the shared area. Concurrent jobs never
  write the same file (faserps rewrites `GDML/FASERCAL_V10.gdml` on every
  run), and a failed job publishes nothing.
* **faserps** writes `<shared>/faserG4/FASERG4-Tcalevent_<run>_<evt>.root`
  (named by the event's own run/event ids, so chunks don't collide). It fails
  if no event file was produced, warns if the count differs from `--n-events`,
  and verifies the copy by size. The first finished chunk also publishes
  `<shared>/GDML/FASERCAL_V10.gdml` (only if absent - delete it if you change
  geometry flags such as `--tilt-deg`).
* **batchreco** stages the chunk's event files plus the geometry into
  scratch, runs events `[start, start+nevt)` (its `--max-event` is
  exclusive), and publishes `<shared>/batch/Batch-TPORecevent_<run>_<start>_<end>.root`.
  It fails if any input event file is missing. Without `--total-events` it
  reconstructs everything faserps produced for `--run` (highest event id + 1).
* **Seeds.** Chunk *k* uses seed `--seed + k` (default base `123456789`).
  Chunk 0 therefore reproduces a plain `run_faserps.py` run; the old
  scripts used one fixed seed for every chunk.
* **One single-threaded process per job.** For more parallelism use more,
  smaller chunks rather than `run_batchreco.py --split` (that launches
  background processes and returns immediately, which doesn't suit a batch
  job). Resources: faserps 1 cpu / 8 GB / 10 GB disk, batchreco 1 cpu / 4 GB /
  6 GB disk - the batchreco memory figure is a guess from a short local run;
  raise it in `submit_batchreco.sub` if jobs get held.
* `KEEP_WORKDIR=1` in a job's environment keeps its scratch directory for
  debugging.

## Not covered / not verified

* **Not tested on lxplus/HTCondor.** The wrappers and drivers were exercised
  against a stand-in checkout (fake `setup.sh`, `run_*.py`, `condor_submit`,
  `ssh`, `rsync`) covering the normal path and the failure paths, not against
  real CVMFS, your build, or `condor_submit`. Run one small chunk first
  (`--total-events 10 --chunk-size 10`) and read `logs/*.out`.
* `setup.sh`'s lxplus branch has a Pythia8 shim in `~/.cache/faser/`. The
  wrappers reuse yours if the node can see your home, else build one in
  scratch (an extra CVMFS lookup per job).
* The `ssh` host/path defaults and the EOS location are guesses - set them as
  above.
* `condor_genie` (GENIE event generation) isn't part of this; the repo's
  `run_convertgenie.py` / `fetch_data.py` handle the GENIE inputs.
