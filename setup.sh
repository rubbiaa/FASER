HOMEFASER=/data/sw/FASERCAL/FASER


##scl enable devtoolset-11 bash

##source /data/sw/FASERCAL/FASER/root_install/bin/thisroot.sh

export PYTHIA8=/data/sw/FASERCAL/FASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=/data/sw/kits/GenFit-install/lib64:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH
