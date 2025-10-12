#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::map<std::string, std::vector<float>>+;
#pragma link C++ class std::vector<std::map<std::string, std::vector<float>>>+;

#pragma link C++ namespace display;
#pragma link C++ class display::FaserCalDisplay+;
//#pragma link C++ class display::FaserCalData+;
#pragma link C++ function display::FaserCalDisplay::SetEventNumber; // Explicit linkage

#pragma link C++ class TPOEvent+;
#pragma link C++ class TcalEvent+;
#pragma link C++ class TPORecoEvent+;
#pragma link C++ class MuTagTrack+;


#pragma link C++ struct PO+;
#pragma link C++ class TPORec+;
#pragma link C++ class DigitizedTrack+;
#pragma link C++ struct TPORecoEvent::PSHIT2D;


#pragma link C++ class MagnetTrack+;
#pragma link C++ class TTKTrack;
#pragma link C++ class TPSTrack;
#pragma link C++ class TPSCluster+;
#pragma link C++ class DBScan+;
#pragma link C++ class TMuTrack+;

#pragma link C++ struct TcalEvent::REARCALDEPOSIT+;

#endif


