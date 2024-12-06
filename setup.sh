HOMEFASER=/home/rubbiaa/FASER

source /home/rubbiaa/ROOT/root_install_v6.32.02/bin/thisroot.sh
source /home/rubbiaa/geant4-install/bin/geant4.sh
export PYTHIA8=/home/rubbiaa/ROOT/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
