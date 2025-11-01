export HOMEFASER=$PWD

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source $HOMEFASER/root-install/bin/thisroot.sh
echo "Root installed in $HOMEFASER/ROOT/root_install"

pushd .
cd /cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb/bin/; source geant4.sh
popd
echo "GEANT4 installed in $GEANT4_INSTALL"

export PYTHIA8=$HOMEFASER/pythia8312
echo "Pythia8 installed in $PYTHIA8"

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
