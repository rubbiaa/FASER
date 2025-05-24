export HOMEFASER=$PWD

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source /home/rubbiaa/ROOT/root_install_v6.32.02/bin/thisroot.sh

echo "Root installed in $HOMEFASER/ROOT/root_install"

GEANT4_INSTALL=/home/rubbiaa/geant4-install/

source $GEANT4_INSTALL/bin/geant4.sh
echo "GEANT4 installed in $GEANT4_INSTALL"

# export PYTHIA8=/home/rubbiaa/ROOT/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
