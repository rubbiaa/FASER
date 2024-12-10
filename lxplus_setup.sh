export HOMEFASER=/eos/home-r/rubbiaa/FASER

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.00/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh
pushd .
cd /cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb/bin/; source geant4.sh
popd

# export PYTHIA8=/home/rubbiaa/ROOT/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
