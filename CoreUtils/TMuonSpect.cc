#include "TMuonSpect.hh"

#include <TTree.h>

ClassImp(TMuonSpectrometer);

void TMuonSpectrometer::Create_Sel_Tree(TTree *t) {
  t->Branch("ntracks", &features.ntracks);
  t->Branch("charge", &features.charge, "charge[ntracks]/F");
  t->Branch("npoints", &features.npoints, "npoints[ntracks]/I");
  t->Branch("px", &features.px, "px[ntracks]/F");
  t->Branch("py", &features.py, "py[ntracks]/F");
  t->Branch("pz", &features.pz, "pz[ntracks]/F");
  t->Branch("p", &features.p, "p[ntracks]/F");
  t->Branch("chi2", &features.chi2, "chi2[ntracks]/F");
  t->Branch("nDoF", &features.nDoF, "nDoF[ntracks]/I");
  t->Branch("pval", &features.pval, "pval[ntracks]/F");
  t->Branch("fpErr", &features.fpErr, "fpErr[ntracks]/F");
  t->Branch("fipErr", &features.fipErr, "fipErr[ntracks]/F");
};
