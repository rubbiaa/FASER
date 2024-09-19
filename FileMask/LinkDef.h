#include <Rtypes.h>

#include "TPORecoEvent.hh"
#include "TTauSearch.hh"
#include "TParticleGun.hh"
#include "TPSCluster.hh"
#include "TTKTrack.hh"
#include "TPSTrack.hh"

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;

#pragma link C++ class TPORecoEvent+;
#pragma link C++ class TPORec+;
#pragma link C++ struct TPORec::CALENERGIES+; 
#pragma link C++ function TPORec::TotalEvis+;
#pragma link C++ function TPORec::TotalET+;

#pragma link C++ class TPOEvent+;
#pragma link C++ class TcalEvent+;
#pragma link C++ class DigitizedTrack+;
#pragma link C++ struct PO+;

#pragma link C++ class TTauSearch;
#pragma link C++ struct TTauSearch::KINEMATICS;

#pragma link C++ class TParticleGun;
#pragma link C++ struct TParticleGun::FEATURES;

#pragma link C++ class TPSCluster;
#pragma link C++ struct TPSCluster::PSCLUSTERLONGPROFILE;

#pragma link C++ class TTKTrack;
#pragma link C++ class TPSTrack;

#pragma link C++ struct TPORecoEvent::REARCALS;

#endif