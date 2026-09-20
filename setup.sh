# Site setup: Ubuntu box (rubbiaa, Ryzen). Sets the handful of things that
# are genuinely specific to this machine, then hands off to
# common_setup.sh for everything shared across all sites.
export HOMEFASER=$PWD

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source /home/rubbiaa/ROOT/root_install_v6.32.02/bin/thisroot.sh
echo "Root installed in $HOMEFASER/ROOT/root_install"

export GEANT4_INSTALL=/home/rubbiaa/geant4-install/
source $GEANT4_INSTALL/bin/geant4.sh
echo "GEANT4 installed in $GEANT4_INSTALL"

export PYTHIA8=/home/rubbiaa/ROOT/pythia8312
echo "Pythia8 installed in $PYTHIA8"

source $HOMEFASER/common_setup.sh
