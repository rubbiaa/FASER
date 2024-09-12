//This file is used to generate the required dictionary for the classes that use the ROOT I/O system
//This file gets read by rootcint and the dictionary is generated
//If new custom classes are added to the code, with the intention to write them to a ROOT file, they must be added here
//Their headers then must be included in the CMakeLists.txt file at the ROOT_GENERATE_DICTIONARY command
#include <vector>
#include <string>
#include <map>
#ifdef __CLING__
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
// #pragma link C++ class Track+;
// #pragma link C++ class vector<Track>+;
#pragma link C++ enum Geant4Process+;
#pragma link C++ enum vector<Geant4Process>+;
#pragma link C++ function operator<<(std::ostream&, const Geant4Process&);
#pragma link C++ class map<std::string, Geant4Process>+;
#pragma link C++ class map<Geant4Process, std::string>+;
#pragma link C++ class map<int, int>+;
#pragma link C++ class DigitizedTrack+;
#pragma link C++ class vector<DigitizedTrack*>+;
#pragma link C++ class TcalEvent+;
#pragma link C++ class TPOEvent+;
#pragma link C++ struct PO+;
#pragma link C++ struct TcalEvent::REARCALDEPOSIT+;
#endif
