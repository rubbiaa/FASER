#include "Track.hh"

#include "G4ios.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"

// ClassImp(Track)

Track::Track() {
	fPosition.clear();
	fMomentum.clear();
	fTime.clear();
	fEnergyDeposit.clear();
//	fProcess.clear();
	fVolume.clear();
}


Track::Track(int TrackID, int ParentID, int PDG, XYZVector Position, XYZVector Momentum, 
						double Time, double EnergyDeposit,
						std::string Volume, int CopyVolume) : Track()
{
    fTotalEDep = EnergyDeposit;
	fTrackID = TrackID;
	fParentID = ParentID;
	fPDG = PDG;

	fPosition.push_back(Position);
	fMomentum.push_back(Momentum);
	fTime.push_back(Time);
	fEnergyDeposit.push_back(EnergyDeposit);
//	fProcess.push_back(Process);
	fVolume.push_back(Volume);

	if(EnergyDeposit > 0) {
		const DetectorConstruction *detector = static_cast<const DetectorConstruction *>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		const XYZVector pos = Position;
		G4long chanID = detector->getChannelIDfromXYZ(Volume, CopyVolume, pos);

		auto it = m_hitIDMap.find(chanID);
		if (it != m_hitIDMap.end())
		{
			it->second += EnergyDeposit;
		}
		else
		{
			m_hitIDMap[chanID] = EnergyDeposit;
		}
	}
}

Track::~Track()
{
	
}


std::vector<XYZVector> Track::getPositionVector() const
{
	std::vector<XYZVector> positionVector;
	for (auto position : fPosition) {
		positionVector.push_back(position);
	}
	return positionVector;
}
std::vector<XYZVector> Track::getMomentumVector() const
{
	std::vector<XYZVector> momentumVector;
	for (auto momentum : fMomentum) {
		momentumVector.push_back(momentum);
	}
	return momentumVector;
}
std::vector<double> Track::getTimeVector() const
{
	std::vector<double> timeVector;
	for (auto time : fTime) {
		timeVector.push_back(time);
	}
	return timeVector;
}
std::vector<double> Track::getEnergyDepositVector() const
{
	std::vector<double> energyDepositVector;
	for (auto energyDeposit : fEnergyDeposit) {
		energyDepositVector.push_back(energyDeposit);
	}
	return energyDepositVector;
}

#if 0
std::vector<Geant4Process> Track::getProcessVector() const
{
	std::vector<Geant4Process> processVector;
	for (auto process : fProcess) {
		processVector.push_back(process);
	}
	return processVector;
}
#endif

std::vector<std::string> Track::getVolumeVector() const
{
	std::vector<std::string> volumeVector;
	for (auto volume : fVolume) {
		volumeVector.push_back(volume);
	}
	return volumeVector;
}
int Track::getTrackID() const { return fTrackID; }
int Track::getParentID() const { return fParentID; }
int Track::getPDG() const { return fPDG; }
double Track::getTotalEnergyDeposit() const { return fTotalEDep; }

XYZVector Track::getPosition(int i) const
{
	if (i < fPosition.size()) {
		return fPosition.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return XYZVector(0, 0, 0);
	}
}

XYZVector Track::getMomentum(int i) const
{
	if (i < fMomentum.size()) {
		return fMomentum.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return XYZVector(0, 0, 0);
	}
}

double Track::getTime(int i) const
{
	if (i < fTime.size()) {
		return fTime.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return 0;
	}
}

double Track::getEnergyDeposit(int i) const
{
	if (i < fEnergyDeposit.size()) {
		return fEnergyDeposit.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return 0;
	}
}

#if 0
Geant4Process Track::getProcess(int i) const
{
	if (i < fProcess.size()) {
		return fProcess.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return NoProcess;
	}
}

Geant4Process Track::getProcess() const
{
	if (fProcess.size() > 0) {
		return fProcess.at(fProcess.size() - 1);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return NoProcess;
	}
}
#endif

std::string Track::getVolume(int i) const
{
	if (i < fVolume.size()) {
		return fVolume.at(i);
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
		return "";
	}
}

void Track::setPositionVector(std::vector<XYZVector> position)
{
	fPosition.clear();
	// copy the vector and then copy the pointer
	//std::vector<XYZVector> *temp = new std::vector<XYZVector>(position);
	//fPosition = temp;
	for (auto pos : position) {
		fPosition.push_back(pos);
	}
}

void Track::setMomentumVector(std::vector<XYZVector> momentum)
{
	fMomentum.clear();
	// copy the vector and then copy the pointer
	//std::vector<XYZVector> *temp = new std::vector<XYZVector>(momentum);
	//fMomentum = temp;
	for (auto mom : momentum) {
		fMomentum.push_back(mom);
	}
}

void Track::setTimeVector(std::vector<double> time)
{
	fTime.clear();
	// copy the vector and then copy the pointer
	//std::vector<double> *temp = new std::vector<double>(time);
	//fTime = temp;
	for (auto t : time) {
		fTime.push_back(t);
	}
}

void Track::setEnergyDepositVector(std::vector<double> energyDeposit)
{
	fEnergyDeposit.clear();

	// copy the vector and then copy the pointer
	//std::vector<double> *temp = new std::vector<double>(energyDeposit);
	//fEnergyDeposit = temp;
	for (auto edep : energyDeposit) {
		fEnergyDeposit.push_back(edep);
	}
}

#if 0
void Track::setProcessVector(std::vector<Geant4Process> process)
{
	fProcess.clear();
	// copy the vector and then copy the pointer
	//std::vector<Geant4Process> *temp = new std::vector<Geant4Process>(process);
	//fProcess = temp;
	for (auto proc : process) {
		fProcess.push_back(proc);
	}
}
#endif

void Track::setVolumeVector(std::vector<std::string> volume)
{
	fVolume.clear();
	// copy the vector and then copy the pointer
	//std::vector<std::string> *temp = new std::vector<std::string>(volume);
	//fVolume = temp;
	for (auto vol : volume) {
		fVolume.push_back(vol);
	}
}

void Track::setTrackID(int trackID) { fTrackID = trackID; }

void Track::setParentID(int parentID) { fParentID = parentID; }

void Track::setPDG(int pdg) { fPDG = pdg; }

void Track::addPosition(XYZVector position) { fPosition.push_back(position); }

void Track::addMomentum(XYZVector momentum) { fMomentum.push_back(momentum); }

void Track::addTime(double time) { fTime.push_back(time); }

void Track::addEnergyDeposit(double energyDeposit) { fEnergyDeposit.push_back(energyDeposit); }

#if 0
void Track::addProcess(Geant4Process process) { fProcess.push_back(process); }
#endif

void Track::addVolume(std::string volume) { fVolume.push_back(volume); }

void Track::update(XYZVector position, XYZVector momentum, double time, double energyDeposit, Geant4Process process, std::string volume, int CopyVolume)
{
	//addPosition(position);
	//addMomentum(momentum);
	//addTime(time);
	//addEnergyDeposit(energyDeposit);
	//addProcess(process);
	//addVolume(volume);

	const DetectorConstruction* detector = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
	const XYZVector pos = position;
	G4long chanID = detector->getChannelIDfromXYZ(volume, CopyVolume, pos);

	auto it = m_hitIDMap.find(chanID);
	if (it != m_hitIDMap.end()) {
			it->second += energyDeposit;
	}
	else {
		m_hitIDMap[chanID] = energyDeposit;
	}
}

void Track::addTotalEnergyDeposit(double EDep) { fTotalEDep += EDep; }

#if 0
void Track::printProcessVector() const
{
	std::cout << "Process Vector: ";
	for (auto process : fProcess) {
		std::cout << process << " ";
	}
	std::cout << std::endl;
}

void Track::printProcess(int i) const
{
	if (i < fProcess.size()) {
		std::cout << "Process: " << fProcess.at(i) << std::endl;
	}
	else {
		std::cout << "Error: Index out of range" << std::endl;
	}
}
#endif

int Track::getNumberOfSteps() const { return fPosition.size(); }

void Track::Dump() {
	G4cout << "Dumping track TrackID = " << fTrackID << " PDG: " << fPDG << " ParentID = " << fParentID << G4endl;

	size_t npositions = fPosition.size();
	size_t nenergydepos = fEnergyDeposit.size();

	if (npositions != nenergydepos) {
        G4cerr << "Warning: Number of positions (" << npositions << ") does not match number of energy deposits (" << nenergydepos << ")." << G4endl;
    }

#if 0
    // Loop over positions and energy deposits
    for (size_t i = 0; i < std::min(npositions, nenergydepos); ++i) {
        G4cout << "Position[" << i << "]: (" << fPosition[i].X() << ", " << fPosition[i].Y() << ", " << fPosition[i].Z() << ")"
               << ", EnergyDeposit[" << i << "]: " << fEnergyDeposit[i] << G4endl;
    }
#endif

	size_t nhits = m_hitIDMap.size();
	G4cout << "Number of hits " << nhits << G4endl;

	for (const auto& hit : m_hitIDMap) {
    	G4cout << "Hit ID: " << hit.first << ", Deposited Energy: " << hit.second << G4endl;
	}

}
