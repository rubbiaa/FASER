# Site setup: André's Mac (Apple Silicon). Sets the handful of things that
# are genuinely specific to this machine, then hands off to
# common_setup.sh for everything shared across all sites.
export HOMEFASER=$PWD

source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

export PYTHIA8=$HOMEFASER/pythia8312

export GEANT4_INSTALL=/Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.2-install/
source $GEANT4_INSTALL/bin/geant4.sh

source $HOMEFASER/common_setup.sh
