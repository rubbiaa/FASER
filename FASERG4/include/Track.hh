#ifndef PARTICLETRAJECTORY_H
#define PARTICLETRAJECTORY_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Geant4Process.hh"
#include "TObject.h"
#include "TVector3.h"
#include "typedef.h"

/**
 * @class Track
 * @brief Object to store the information of a single particle track
 * @details This class is used to store the information of a particle track.
 * The information stored is two-fold. Once the static information of the particle like the TrackID, its ParentId and the PDG code.
 * The second part is the dynamic information of the particle like the position, momentum, energy, process, Volumen name and time.
 * It contains the corresponding getter and setter functions.
 */

class Track : public TObject {
    private:
	std::vector<XYZVector> fPosition;	 ///< Position of the particle at each step
	std::vector<XYZVector> fMomentum;	 ///< Momentum of the particle at each step
	std::vector<double> fTime;		 ///< Time of the particle at each step
	std::vector<double> fEnergyDeposit;	 ///< Energy deposit of the particle at each step
//	std::vector<Geant4Process> fProcess;	 ///< Process that caused the particle to stop
	std::vector<std::string> fVolume;	 ///< Volume in which the particle is at each step
	int fTrackID;					 ///< Track ID of the particle
	int fParentID;					 ///< Track ID of the parent particle
	int fPDG;					 ///< PDG code of the particle
	double fTotalEDep;					 ///< Total energy deposit of the particle

	std::map<long, double> m_hitIDMap;  ///< Map to store hit ID and corresponding deposited energy

    public:
	/**
	 * @brief Default constructor, all values will be set to be empty
	 * @note Potential undefined behaviour when using this constructor, and calling the getter methods without setting the values first.
	 */
	Track();

	/**
	 * @brief Constructor with all values
	 * @param TrackID Track ID of the particle
	 * @param ParentID Track ID of the parent particle
	 * @param PDG PDG code of the particle
	 * @param Position Position of the particle at each step
	 * @param Momentum Momentum of the particle at each step
	 * @param Time Time of the particle at each step
	 * @param EnergyDeposit Energy deposit of the particle at each step
	 * @param Volume Volume in which the particle is at each step
	 * @param CopyVolume volume replica number (if any)
	 */
	Track(int TrackID, int ParentID, int PDG, XYZVector Position, XYZVector Momentum, double Time, double EnergyDeposit,
	      std::string Volume, int CopyVolume);

	~Track();  ///< Default destructor

	/**
	 * @brief Get the full vector of all positions
	 * @return Vector of all positions std::vector<XYZVector>
	 */
	std::vector<XYZVector> getPositionVector() const;

	/**
	 * @brief Get the full vector of all momenta
	 * @return Vector of all momenta std::vector<XYZVector>
	 */
	std::vector<XYZVector> getMomentumVector() const;

	/**
	 * @brief Get the full vector of all times
	 * @return Vector of all times std::vector<double>
	 */
	std::vector<double> getTimeVector() const;

	/**
	 * @brief Get the full vector of all energy deposits
	 * @return Vector of all energy deposits std::vector<double>
	 */
	std::vector<double> getEnergyDepositVector() const;

	/**
	 * @brief Get the full vector of all processes
	 * @return Vector of all processes std::vector<Geant4Process>
	 */
	std::vector<Geant4Process> getProcessVector() const;

	/**
	 * @brief Get the full vector of all volumes
	 * @return Vector of all volumes std::vector<std::string>
	 */
	std::vector<std::string> getVolumeVector() const;

	/**
	 * @brief get the TrackID of the particle
	 * @return TrackID int
	 */
	int getTrackID() const;
	/**
	 * @brief get the ParentID of the particle
	 * @return ParentID int
	 */
	int getParentID() const;
	/**
	 * @brief get the PDG code of the particle
	 * @return PDG int
	 */
	int getPDG() const;

	/**
	 * @brief get the total energy deposit of the particle
	 * @return TotalEnergyDeposit double
	 */
	double getTotalEnergyDeposit() const;

	/**
	 * @brief get the position of the particle at a specific step
	 * @param i step of the particle
	 * @return Position XYZVector
	 */
	XYZVector getPosition(int i) const;

	/**
	 * @brief get the last step position of the particle
	 * @return Position XYZVector
	 */
	XYZVector getPosition() const;

	/**
	 * @brief get the momentum of the particle at a specific step
	 * @param i step of the particle
	 * @return Momentum XYZVector
	 */
	XYZVector getMomentum(int i) const;

	/**
	 * @brief get the last step momentum of the particle
	 * @return Momentum XYZVector
	 */
	XYZVector getMomentum() const;

	/**
	 * @brief get the time of the particle at a specific step
	 * @param i step of the particle
	 * @return Time double
	 */
	double getTime(int i) const;

	/**
	 * @brief get the last step time of the particle
	 * @return Time double
	 */
	double getTime() const;

	/**
	 * @brief get the energy deposit of the particle at a specific step
	 * @param i step of the particle
	 * @return EnergyDeposit double
	 */
	double getEnergyDeposit(int i) const;

	/**
	 * @brief get the last step energy deposit of the particle
	 * @return EnergyDeposit double
	 */
	double getEnergyDeposit() const;

	/**
	  * @brief get the process of the particle at a specific step
	  * @param i step of the particle
	  * @return Process Geant4Process
`	  * @note Returns the process of the track at step i, returns Geant4Process::NoProcess if i is out of range, which however could also be a
legit process. Prints warning to cout
	  */
	Geant4Process getProcess(int i) const;

	/**
	 * @brief get the last step process of the particle
	 * @return Process Geant4Process
	 */
	Geant4Process getProcess() const;

	/**
	 * @brief get the volume of the particle at a specific step
	 * @param i step of the particle
	 * @return Volume std::string
	 */
	std::string getVolume(int i) const;

	/**
	 * @brief get the last step volume of the particle
	 * @return Volume std::string
	 */
	std::string getVolume() const;

	/**
	 * @brief get the number of steps of the particle
	 * @return Number of steps int
	 */
	int getNumberOfSteps() const;

	std::map<long, double> getmhitIDMap() const { return m_hitIDMap;};

	/**
	 * @brief Set the Position vector
	 * @param position std::vector<XYZVector>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setPositionVector(std::vector<XYZVector> position);

	/**
	 * @brief Set the Momentum vector
	 * @param momentum std::vector<XYZVector>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setMomentumVector(std::vector<XYZVector> momentum);

	/**
	 * @brief Set the Time vector
	 * @param time std::vector<double>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setTimeVector(std::vector<double> time);

	/**
	 * @brief Set the Energy Deposit vector
	 * @param energyDeposit std::vector<double>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setEnergyDepositVector(std::vector<double> energyDeposit);

	/**
	 * @brief Set the Process vector
	 * @param process std::vector<Geant4Process>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setProcessVector(std::vector<Geant4Process> process);

	/**
	 * @brief Set the TrackID
	 * @param trackID int
	 * @note Public function but is not really intended to be used by the user
	 */
	void setTrackID(int trackID);

	/**
	 * @brief Set the Parent ID of the track
	 * @param parentID int
	 * @note Public function but is not really intended to be used by the user
	 */
	void setParentID(int parentID);

	/**
	 * @brief Set the PDG code of the track
	 * @param pdg int
	 * @note Public function but is not really intended to be used by the user
	 */
	void setPDG(int pdg);

	/**
	 * @brief Set the Volume vector
	 * @param volume std::vector<std::string>
	 * @note Public function but is not really intended to be used by the user
	 */
	void setVolumeVector(std::vector<std::string> volume);

    private:
	/**
	 * @brief Appends the position vector
	 * @param position XYZVector
	 */
	void addPosition(XYZVector position);

	/**
	 * @brief Appends the momentum vector
	 * @param momentum XYZVector
	 */
	void addMomentum(XYZVector momentum);

	/**
	 * @brief Appends the time vector
	 * @param time double
	 */
	void addTime(double time);

	/**
	 * @brief Appends the energy deposit vector
	 * @param energyDeposit double
	 */
	void addEnergyDeposit(double energyDeposit);

	/**
	 * @brief Appends the process vector
	 * @param process Geant4Process
	 */
	void addProcess(Geant4Process process);

	/**
	 * @brief Appends the volume vector
	 * @param volume std::string
	 */
	void addVolume(std::string volume);

    public:
	/**
	 * @brief Updates the track with the information of a new step
	 * @details This method is used to update the track with the information of a new step. It is called by the ParticleManager class.
	 * @param position Position of the particle at the new step (XYZVector)
	 * @param momentum Momentum of the particle at the new step (XYZVector)
	 * @param time Time of the particle at the new step (double)
	 * @param energyDeposit Energy deposit of the particle at the new step (double)
	 * @param volume Volume in which the particle is at the new step (std::string)
	 * @param CopyVolume volume replica number (if any)
	 */
	void update(XYZVector position, XYZVector momentum, double time, double energyDeposit, 
			std::string volume, int CopyVolume);

	/**
	 * @brief Add to the total energy deposit of the particle
	 * @param totalEDep double
	 */
	void addTotalEnergyDeposit(double totalEDep);

	/**
	 * @brief Prints the Process vector of the track
	 * @note Required as the Process enum is not known to ROOT
	 */
	void printProcessVector() const;

	/**
	 * @brief Prints the process at step i of the track
	 * @param i Step of the track
	 * @note Required as the Process enum is not known to ROOT
	 */
	void printProcess(int i) const;

	void Dump();

//	ClassDef(Track, 1);
};

#endif
