#ifndef G4PROCESS_H
#define G4PROCESS_H 1
#pragma once

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

/**
 * @enum Geant4Process
 * @brief Enumeration of Geant4 processes to store in the output tree
 * @note Should be updated when new proccesses are added or simulated
 */
enum Geant4Process {
	NoProcess,	       ///< No process
	primary,	       ///< Primary particle creation
	eIoni,		       ///< Ionisation
	msc,		       ///< Multiple scattering
	compt,		       ///< Compton scattering
	phot,		       ///< Photoelectric effect
	eBrem,		       ///< Bremsstrahlung
	ionIoni,	       ///< Ionisation by ions
	hIoni,		       ///< Hadronic ionisation
	RadioactiveDecayBase,  ///< Radioactive decay
	CoulombScat,	       ///< Coulomb scattering
	Rayl,		       ///< Rayleigh scattering
	Transportation,	       ///< Transportation between Volumes
	annihil,	       ///< Annihilation
	conv,		       ///< Conversion
	hadElastic,	       ///< Hadronic elastic scattering
	nCapture,	       ///< Neutron capture
	neutronInelastic,      ///< Neutron inelastic scattering
	photonNuclear,	       ///< Photonuclear interaction
	protonInelastic,       ///< Proton inelastic scattering
	dInelastic,	       ///< Deuteron inelastic scattering
	Decay,		       ///< Decay
	RadioactiveDecay,      ///< Radioactive decay
	Scinitllation,	       ///< Scintillation light generation
	Cerenkov	       ///< Cerenkov light generation
};

/**
 * @brief A map between the Geant4Process enum and a string representation
 */
static std::map<std::string, Geant4Process> Geant4ProcessMap = {{"NoProcess", NoProcess},
								{"primary", primary},
								{"eIoni", eIoni},
								{"msc", msc},
								{"compt", compt},
								{"phot", phot},
								{"eBrem", eBrem},
								{"ionIoni", ionIoni},
								{"hIoni", hIoni},
								{"RadioactiveDecayBase", RadioactiveDecayBase},
								{"CoulombScat", CoulombScat},
								{"Rayl", Rayl},
								{"Transportation", Transportation},
								{"annihil", annihil},
								{"conv", conv},
								{"hadElastic", hadElastic},
								{"nCapture", nCapture},
								{"neutronInelastic", neutronInelastic},
								{"photonNuclear", photonNuclear},
								{"protonInelastic", protonInelastic},
								{"dInelastic", dInelastic},
								{"Decay", Decay},
								{"RadioactiveDecay", RadioactiveDecay},
								{"Scinitllation", Scinitllation},
								{"Cerenkov", Cerenkov}};

/**
 * @brief A map between the Geant4Process enum and a string representation
 */
static std::map<Geant4Process, std::string> Geant4ProcessMapInverse = {{NoProcess, "NoProcess"},
								       {primary, "primary"},
								       {eIoni, "eIoni"},
								       {msc, "msc"},
								       {compt, "compt"},
								       {phot, "phot"},
								       {eBrem, "eBrem"},
								       {ionIoni, "ionIoni"},
								       {hIoni, "hIoni"},
								       {RadioactiveDecayBase, "RadioactiveDecayBase"},
								       {CoulombScat, "CoulombScat"},
								       {Rayl, "Rayl"},
								       {Transportation, "Transportation"},
								       {annihil, "annihil"},
								       {conv, "conv"},
								       {hadElastic, "hadElastic"},
								       {nCapture, "nCapture"},
								       {neutronInelastic, "neutronInelastic"},
								       {photonNuclear, "photonNuclear"},
								       {protonInelastic, "protonInelastic"},
								       {dInelastic, "dInelastic"},
								       {Decay, "Decay"},
								       {RadioactiveDecay, "RadioactiveDecay"},
								       {Scinitllation, "Scinitllation"},
								       {Cerenkov, "Cerenkov"}};

/**
 * @brief A print function for the Geant4Process enum
 */
inline std::ostream& operator<<(std::ostream& os, const Geant4Process& p)
{
	if (Geant4ProcessMapInverse.find(p) != Geant4ProcessMapInverse.end()) {
		os << Geant4ProcessMapInverse[p];
	}
	else {
		os << "Unknown Geant4Process";
	}
	return os;
}

#endif
