//////////(o^o)///////////
#include "GenMagneticField.hh"
#include <iostream>
#include <cstdlib>
#include <string>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>

namespace {
void TraverseAndCollect(TGeoNode* node,
												std::vector<std::pair<double,double>>& magnet_z_ranges_cm,
												std::vector<double>& scifi_layer_z_cm,
												int verbose,
												const std::string& indent = "") {
	if (!node) return;
	TGeoVolume* vol = node->GetVolume();
	if (!vol) return;
	const char* cname = node->GetName();
	std::string name = cname ? std::string(cname) : std::string();
	TGeoBBox* shape = dynamic_cast<TGeoBBox*>(vol->GetShape());
	if (shape) {
		const double zc_mm = node->GetMatrix()->GetTranslation()[2];
		const double dz_mm = shape->GetDZ();
		const double zc_cm = zc_mm / 10.0;
		const double dz_cm = dz_mm / 10.0;
		if (name.find("Magnet") != std::string::npos) {
			magnet_z_ranges_cm.push_back({zc_cm - dz_cm, zc_cm + dz_cm});
			if (verbose >= 2) std::cout << "[Geom] Magnet range (cm): " << zc_cm - dz_cm << " to " << zc_cm + dz_cm << " name=" << name << std::endl;
		}
		if (name.find("SciFiLayer") != std::string::npos || name.find("InitialSciFiLayer") != std::string::npos || name.find("SciFi") != std::string::npos) {
			scifi_layer_z_cm.push_back(zc_cm);
			if (verbose >= 3) std::cout << "[Geom] SciFi layer z (cm): " << zc_cm << " name=" << name << std::endl;
		}
	}
	// Recurse into daughters
	TObjArray* daughters = vol->GetNodes();
	if (!daughters) return;
	for (int i = 0; i < daughters->GetEntries(); ++i) {
		TGeoNode* child = (TGeoNode*) daughters->At(i);
		TraverseAndCollect(child, magnet_z_ranges_cm, scifi_layer_z_cm, verbose, indent + "  ");
	}
}
}
//////////(o^o)///////////
// Build geometry-driven magnetic field maps
void GenMagneticField::BuildGeometryMaps(int verbose) {
	magnet_z_ranges_cm_.clear();
	scifi_layer_z_cm_.clear();
	int venv = 0;
	if (const char* s = std::getenv("MS_VERBOSE")) { try { venv = std::stoi(s); } catch(...) { venv = 1; } }
	if (venv > verbose) verbose = venv;

	if (!gGeoManager) {
		if (verbose >= 1) std::cout << "[Geom] No TGeoManager; geometry-driven field maps disabled" << std::endl;
		return;
	}

	TGeoNode* top = gGeoManager->GetTopNode();
	if (!top) {
		if (verbose >= 1) std::cout << "[Geom] No top node in geometry" << std::endl;
		return;
	}

	TraverseAndCollect(top, magnet_z_ranges_cm_, scifi_layer_z_cm_, verbose);
	std::sort(magnet_z_ranges_cm_.begin(), magnet_z_ranges_cm_.end());
	std::sort(scifi_layer_z_cm_.begin(), scifi_layer_z_cm_.end());

	if (verbose >= 1) {
		std::cout << "[Geom] Collected " << magnet_z_ranges_cm_.size() << " magnet ranges and "
							<< scifi_layer_z_cm_.size() << " SciFi layers" << std::endl;
	}
	if (verbose >= 2) {
		std::cout << "[Geom] Magnet Z-ranges (cm):" << std::endl;
		for (auto& r : magnet_z_ranges_cm_) std::cout << "  " << r.first << " to " << r.second << std::endl;
		std::cout << "[Geom] SciFiLayer Z positions (cm):";
		for (auto z : scifi_layer_z_cm_) std::cout << " " << z;
		std::cout << std::endl;
	}
}

bool GenMagneticField::IsInMagnet(double z_cm) const {
	for (const auto& rng : magnet_z_ranges_cm_) if (z_cm > rng.first && z_cm < rng.second) return true;
	return false;
}

bool GenMagneticField::IsNearStation(double z_cm, double eps_cm) const {
	for (double zl : scifi_layer_z_cm_) if (std::abs(zl - z_cm) <= eps_cm) return true;
	return false;
}
//////////(o^o)///////////

