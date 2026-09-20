# FASERCAL - electronic calorimeter for FASER Run4

FASERCAL code to simulate and analyse events in the FASERCAL detector

## Lxplus Installation

```bash
source setup.sh # auto-detects lxplus and sets up root, geant4, and the LD_LIBRARY_PATH
make pythia8
make clhep
make rave
make googletest
make genfit
cp GenFit-build/bin/*.pcm GenFit-install/lib64/ # this is necessary to avoid errors of the type ".pcm file not found" in ROOT
```

### Example Usage (Batch reconstruction)
```bash
cd Batch
ln -sfn /eos/project-f/faser-upgrade/fasercal/FASERCALDATA_v5.1 input
make
./batchreco.exe 130 0 100 -mt # -mt option activates multi-threading
```
