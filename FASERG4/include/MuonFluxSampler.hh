#ifndef MuonFluxSampler_h
#define MuonFluxSampler_h 1

#include <string>
#include <vector>

/// @brief Samples muon charge (mu-/mu+) and energy from a tabulated cosmic/beam-induced
/// muon flux at FASER.
///
/// The flux is expected in the FASERnu Run 3 FLUKA muon flux grid format bundled with the
/// muDIS charm-asymmetry study (see the FASER muDIS Dropbox folder,
/// faser-muon-dis-charm-asymmetry/data/muon_flux_FASERv_Run3_var2_0000.dat): a repurposed
/// LHAPDF "lhagrid1" grid, where the "x" variable is E_mu / kReferenceBeamEnergyGeV (the
/// nominal 7 TeV LHC beam energy) and the stored quantity is x*flux(x) for PDG codes -13/13
/// (mu+/mu-), following the same x*f(x) convention as an ordinary PDF grid. The grid's Q
/// dimension is a format placeholder only (duplicated identical rows) and carries no physics;
/// it is read but otherwise ignored here.
///
/// This is a singleton, mirroring the MuonDISPythiaGenerator pattern used elsewhere in
/// FASERG4: load the grid once (loadFromFile is a no-op if the same path is already loaded),
/// then call sample() once per generated event.
class MuonFluxSampler {
public:
  static MuonFluxSampler& instance();

  /// Load (or reload) the flux grid from the given lhagrid1-format .dat file. No-op if the
  /// grid at this exact path is already loaded. Aborts the run (G4Exception, FatalException)
  /// if the file cannot be opened or does not parse as expected -- silently falling back would
  /// let a run continue with a mismatched/nonsensical flux.
  void loadFromFile(const std::string& path);

  /// Sample a muon PDG code (+-13) and energy (GeV) from the loaded flux, weighting mu-/mu+
  /// by their respective integrated flux. Returns false (leaving outputs untouched) if no grid
  /// has been loaded yet.
  bool sample(int& pdgId, double& energyGeV) const;

  bool isLoaded() const { return m_loaded; }

private:
  MuonFluxSampler() = default;

  struct Species {
    std::vector<double> x;         // grid points, x = E_mu / kReferenceBeamEnergyGeV
    std::vector<double> cumWeight; // cumulative integral of f(x) = (stored x*f(x))/x over x
    double totalWeight = 0.0;      // cumWeight.back(), or 0 if this species has no flux
  };

  double sampleX(const Species& sp) const;

  Species m_muMinus; // PDG 13
  Species m_muPlus;  // PDG -13
  bool m_loaded = false;
  std::string m_loadedPath;

  // Nominal LHC beam energy (GeV) the flux grid's "x" variable is normalized to -- matches
  // the convention used when this grid was built for the muDIS charm-asymmetry study.
  static constexpr double kReferenceBeamEnergyGeV = 7000.0;
};

#endif
