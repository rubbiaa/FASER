HOMEFASER=/Users/rubbiaa/Documents/GitHub/FASER

source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

export PYTHIA8=$HOMEFASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH

