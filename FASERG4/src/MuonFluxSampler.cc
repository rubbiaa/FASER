#include "MuonFluxSampler.hh"

#include "G4Exception.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace {

std::vector<double> parseDoubles(const std::string& line) {
  std::vector<double> values;
  std::istringstream iss(line);
  double v;
  while (iss >> v) {
    values.push_back(v);
  }
  return values;
}

}  // namespace

MuonFluxSampler& MuonFluxSampler::instance() {
  static MuonFluxSampler sampler;
  return sampler;
}

void MuonFluxSampler::loadFromFile(const std::string& path) {
  if (m_loaded && m_loadedPath == path) {
    return;
  }

  std::ifstream in(path);
  if (!in.is_open()) {
    G4Exception("MuonFluxSampler::loadFromFile", "MuonFluxSamplerNoFile", FatalException,
                ("Cannot open muon flux grid file: " + path).c_str());
    return;
  }

  // Skip the lhagrid1 header down to the first "---" separator.
  std::string line;
  bool foundSeparator = false;
  while (std::getline(in, line)) {
    if (line.find("---") != std::string::npos) {
      foundSeparator = true;
      break;
    }
  }
  if (!foundSeparator) {
    G4Exception("MuonFluxSampler::loadFromFile", "MuonFluxSamplerBadFormat", FatalException,
                ("No lhagrid1 '---' separator found in: " + path).c_str());
    return;
  }

  std::getline(in, line);
  const std::vector<double> xs = parseDoubles(line);  // x = E_mu / kReferenceBeamEnergyGeV

  std::getline(in, line);
  const std::vector<double> qs = parseDoubles(line);  // format placeholder, not used

  std::getline(in, line);
  const std::vector<double> flavorsRaw = parseDoubles(line);
  const std::vector<int> flavors(flavorsRaw.begin(), flavorsRaw.end());

  int idxMuMinus = -1;
  int idxMuPlus = -1;
  for (size_t i = 0; i < flavors.size(); ++i) {
    if (flavors[i] == 13) idxMuMinus = static_cast<int>(i);
    if (flavors[i] == -13) idxMuPlus = static_cast<int>(i);
  }
  if (idxMuMinus < 0 || idxMuPlus < 0 || xs.empty() || qs.empty()) {
    G4Exception("MuonFluxSampler::loadFromFile", "MuonFluxSamplerNoMuonFlavor", FatalException,
                ("Muon flux grid does not contain both +-13 flavors: " + path).c_str());
    return;
  }

  // Data rows: len(xs)*len(qs) rows, one per (x,Q) pair, len(flavors) columns each. The Q
  // dimension carries no physics here (identical rows duplicated across Q, see header
  // comment in MuonFluxSampler.hh) -- only the first Q slice (qIndex == 0) is kept.
  const size_t nRows = xs.size() * qs.size();
  std::vector<double> xfMuMinus(xs.size(), 0.0);
  std::vector<double> xfMuPlus(xs.size(), 0.0);
  for (size_t r = 0; r < nRows; ++r) {
    if (!std::getline(in, line)) {
      G4Exception("MuonFluxSampler::loadFromFile", "MuonFluxSamplerTruncated", FatalException,
                  ("Unexpected end of file while reading flux grid: " + path).c_str());
      return;
    }
    const size_t qIndex = r % qs.size();
    if (qIndex != 0) {
      continue;
    }
    const size_t xIndex = r / qs.size();
    const std::vector<double> row = parseDoubles(line);
    const size_t maxIdx = static_cast<size_t>(std::max(idxMuMinus, idxMuPlus));
    if (row.size() <= maxIdx) {
      continue;
    }
    xfMuMinus[xIndex] = row[static_cast<size_t>(idxMuMinus)];
    xfMuPlus[xIndex] = row[static_cast<size_t>(idxMuPlus)];
  }

  auto buildSpecies = [&xs](const std::vector<double>& xf) {
    Species sp;
    sp.x = xs;
    sp.cumWeight.assign(xs.size(), 0.0);
    double cum = 0.0;
    double fPrev = 0.0;
    for (size_t i = 0; i < xs.size(); ++i) {
      // The grid stores x*f(x) (PDF-style convention); undo it to recover the physical
      // density f(x) = dN/dx before integrating.
      const double f = (xs[i] > 0.0) ? xf[i] / xs[i] : 0.0;
      if (i > 0) {
        cum += 0.5 * (f + fPrev) * (xs[i] - xs[i - 1]);  // trapezoidal integral
      }
      sp.cumWeight[i] = cum;
      fPrev = f;
    }
    sp.totalWeight = cum;
    return sp;
  };

  m_muMinus = buildSpecies(xfMuMinus);
  m_muPlus = buildSpecies(xfMuPlus);
  m_loaded = true;
  m_loadedPath = path;

  G4cout << "MuonFluxSampler: loaded flux grid from " << path << " (mu- weight="
         << m_muMinus.totalWeight << ", mu+ weight=" << m_muPlus.totalWeight << ")" << G4endl;
}

double MuonFluxSampler::cumAtX(const Species& sp, double x) const {
  if (sp.x.empty()) {
    return 0.0;
  }
  if (x <= sp.x.front()) {
    return 0.0;
  }
  if (x >= sp.x.back()) {
    return sp.totalWeight;
  }
  size_t hi = 1;
  while (hi < sp.x.size() && sp.x[hi] < x) {
    ++hi;
  }
  const size_t lo = hi - 1;
  const double xlo = sp.x[lo];
  const double xhi = sp.x[hi];
  const double wlo = sp.cumWeight[lo];
  const double whi = sp.cumWeight[hi];
  const double frac = (xhi > xlo) ? (x - xlo) / (xhi - xlo) : 0.0;
  return wlo + frac * (whi - wlo);
}

double MuonFluxSampler::sampleXAbove(const Species& sp, double cumLo) const {
  const double target = cumLo + G4UniformRand() * (sp.totalWeight - cumLo);
  size_t hi = 1;
  while (hi < sp.cumWeight.size() && sp.cumWeight[hi] < target) {
    ++hi;
  }
  if (hi >= sp.cumWeight.size()) {
    hi = sp.cumWeight.size() - 1;
  }
  const size_t lo = hi - 1;
  const double wlo = sp.cumWeight[lo];
  const double whi = sp.cumWeight[hi];
  const double frac = (whi > wlo) ? (target - wlo) / (whi - wlo) : 0.0;
  return sp.x[lo] + frac * (sp.x[hi] - sp.x[lo]);
}

bool MuonFluxSampler::sample(int& pdgId, double& energyGeV, double minEnergyGeV) const {
  if (!m_loaded) {
    return false;
  }

  // Restrict each species' sampling range to the portion of its flux at/above the cutoff,
  // so both the mu-/mu+ mix and the energy spectrum are correctly renormalized to the
  // surviving (>= minEnergyGeV) population, rather than rejection-sampling the full flux
  // (which would be arbitrarily wasteful for a high cutoff that discards most of the flux).
  const double minX = (minEnergyGeV > 0.0) ? minEnergyGeV / kReferenceBeamEnergyGeV : 0.0;

  const double cumMinusLo = cumAtX(m_muMinus, minX);
  const double cumPlusLo = cumAtX(m_muPlus, minX);
  const double weightMinus = m_muMinus.totalWeight - cumMinusLo;
  const double weightPlus = m_muPlus.totalWeight - cumPlusLo;
  const double totalWeight = weightMinus + weightPlus;
  if (totalWeight <= 0.0) {
    return false;  // no flux at/above this cutoff
  }

  const bool isMuMinus = (G4UniformRand() * totalWeight) < weightMinus;
  const Species& sp = isMuMinus ? m_muMinus : m_muPlus;
  const double cumLo = isMuMinus ? cumMinusLo : cumPlusLo;

  const double x = sampleXAbove(sp, cumLo);
  pdgId = isMuMinus ? 13 : -13;
  // Guard against the cutoff landing exactly on a grid point plus floating-point roundoff
  // nudging the sampled x a hair below minX.
  energyGeV = std::max(x, minX) * kReferenceBeamEnergyGeV;
  return true;
}
