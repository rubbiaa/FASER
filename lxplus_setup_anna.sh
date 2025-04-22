export HOMEFASER=/afs/cern.ch/work/a/amascell/FASERCal/debug/FASER

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
source /cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb/bin/geant4.sh

export PYTHIA8=$HOMEFASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH