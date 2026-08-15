#HOMEFASER=/Users/rubbiaa/Documents/GitHub/FASER
HOMEFASER=/Users/rubbiaa/MACDEV/fasermuondis

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

GEANT4_INSTALL=/Users/rubbiaa/Documents/GitHub/FASER/geant4-install 

source $GEANT4_INSTALL/bin/geant4.sh
echo "GEANT4 installed in $GEANT4_INSTALL"
echo "Setting up environment for FASER simulation"

export PYTHIA8=$HOMEFASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH

