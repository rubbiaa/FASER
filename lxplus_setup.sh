export HOMEFASER=/eos/home-r/rubbiaa/FASER

source $HOMEFASER/root-install/bin/thisroot.sh
pushd .
cd /cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb/bin/; source geant4.sh
popd

export PYTHIA8=$HOMEFASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
