// create signal and background summaries
// A. Rubbia/June 2024
//

#include <TRandom3.h>

TRandom3 randomGen(0); // 0 means seed with current time

// summary of event
struct SUMEVENT {
  Bool_t         isCC;
  Bool_t         isES;
  Bool_t         istau;
  Float_t        prim_vx;
  Float_t        prim_vy;
  Float_t        prim_vz;
  Bool_t         infid;
  Int_t          tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
  Int_t          n_particles;
  Int_t          n_charged;
  Int_t          in_lepton_pdgid;
  Float_t        vis_spx;   // total visible energy
  Float_t        vis_spy;
  Float_t        vis_spz;
  Float_t        jetpx;     // jet energy
  Float_t        jetpy;
  Float_t        jetpz;
  Float_t        tauvis_px;  // tau decay products visible energy
  Float_t        tauvis_py;
  Float_t        tauvis_pz;
  Float_t        tautracklength;
  Float_t        Evis;
  Float_t        ptmiss;
  Float_t        cost;
  Float_t        cosf;
  Float_t        plep;
  Float_t        ptlep;
  Float_t        qtlep;

    // selection
  Float_t        taupi_cand[3];
  Float_t        taurho_cand[3];

} sumevent;

struct stats {
  size_t infid;
  size_t nueCC;
  size_t numuCC;
  size_t nutauCC;
  size_t nutauCClost;
  size_t NC;
  size_t ES;
} stats;

TChain *event_tree;

// tau->e
TFile *nueCC_signal_file;
TTree *nueCC_signal_tree;
TFile *nueCC_bkg_file;
TTree *nueCC_bkg_tree;

// tau->mu
TFile *numuCC_signal_file;
TTree *numuCC_signal_tree;
TFile *numuCC_bkg_file;
TTree *numuCC_bkg_tree;

// tau->pi
TFile *nutaupi_signal_file;
TTree *nutaupi_signal_tree;
TFile *nutaupi_bkg_file;
TTree *nutaupi_bkg_tree;

// tau->rho
TFile *nutaurho_signal_file;
TTree *nutaurho_signal_tree;
TFile *nutaurho_bkg_file;
TTree *nutaurho_bkg_tree;

double compute_momentum(float p[3]) {
  return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}
double compute_pt(float p[3]) {
  return sqrt(p[0]*p[0]+p[1]*p[1]);
}
double compute_qtlep(float p[3], float jet[3]) {
  // projection of momentum on hadronic jet
  double pproj = (p[0]*jet[0] + p[1]*jet[1] + p[2]*jet[2])/sqrt(jet[0]*jet[0]+jet[1]*jet[1]+jet[2]*jet[2]);
  double p2 = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  return sqrt(p2-pproj*pproj);
}

double compute_cost(float plep[3], float jet[3]) {
  return (plep[0]*jet[0] + plep[1]*jet[1])/(sqrt(plep[0]*plep[0]+plep[1]*plep[1])*sqrt(jet[0]*jet[0]+jet[1]*jet[1]));
}

void create_sel_tree(TTree *t) {
  t->Branch("infid",&sumevent.infid);
  t->Branch("evis",&sumevent.Evis);
  t->Branch("plep",&sumevent.plep);
  t->Branch("ptlep",&sumevent.ptlep);
  t->Branch("qtlep",&sumevent.qtlep);
  t->Branch("ptmiss",&sumevent.ptmiss);
  t->Branch("cost",&sumevent.cost);
  t->Branch("cosf",&sumevent.cosf);
};

void open_tuples() {
  event_tree = new TChain("event_tree");  
  event_tree->Add("event_data_test.root");
  
  std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

  // Set branch addresses.
  
   event_tree->SetBranchAddress("isCC",&sumevent.isCC);
   event_tree->SetBranchAddress("isES",&sumevent.isES);
   event_tree->SetBranchAddress("istau",&sumevent.istau);
   event_tree->SetBranchAddress("tau_decaymode",&sumevent.tau_decaymode);
   event_tree->SetBranchAddress("n_particles",&sumevent.n_particles);
   event_tree->SetBranchAddress("n_charged",&sumevent.n_charged);
   event_tree->SetBranchAddress("prim_vx", &sumevent.prim_vx);
   event_tree->SetBranchAddress("prim_vy", &sumevent.prim_vy);
   event_tree->SetBranchAddress("prim_vz", &sumevent.prim_vz);
   event_tree->SetBranchAddress("in_lepton_pdgid",&sumevent.in_lepton_pdgid);
   event_tree->SetBranchAddress("vis_spx",&sumevent.vis_spx);
   event_tree->SetBranchAddress("vis_spy",&sumevent.vis_spy);
   event_tree->SetBranchAddress("vis_spz",&sumevent.vis_spz);
   event_tree->SetBranchAddress("jetpx",&sumevent.jetpx);
   event_tree->SetBranchAddress("jetpy",&sumevent.jetpy);
   event_tree->SetBranchAddress("jetpz",&sumevent.jetpz);
   event_tree->SetBranchAddress("tauvis_px",&sumevent.tauvis_px);
   event_tree->SetBranchAddress("tauvis_py",&sumevent.tauvis_py);
   event_tree->SetBranchAddress("tauvis_pz",&sumevent.tauvis_pz);
   event_tree->SetBranchAddress("tautracklength",&sumevent.tautracklength);
   event_tree->SetBranchAddress("Evis",&sumevent.Evis);
   event_tree->SetBranchAddress("ptmiss",&sumevent.ptmiss);
   event_tree->SetBranchAddress("taupi_cand[3]",sumevent.taupi_cand);
   event_tree->SetBranchAddress("taurho_cand[3]",sumevent.taurho_cand);

   // tau -> e
   nueCC_signal_file = new TFile("nuecc_signal.root","RECREATE");
   nueCC_signal_tree = new TTree("nuecc_signal","Event data");
   create_sel_tree(nueCC_signal_tree);
   
   nueCC_bkg_file = new TFile("nuecc_bkg.root","RECREATE");
   nueCC_bkg_tree = new TTree("nuecc_bkg","Event data");
   create_sel_tree(nueCC_bkg_tree);

   // tau -> mu
   numuCC_signal_file = new TFile("numucc_signal.root","RECREATE");
   numuCC_signal_tree = new TTree("numucc_signal","Event data");
   create_sel_tree(numuCC_signal_tree);
   
   numuCC_bkg_file = new TFile("numucc_bkg.root","RECREATE");
   numuCC_bkg_tree = new TTree("numucc_bkg","Event data");
   create_sel_tree(numuCC_bkg_tree);

   // tau -> pi
   nutaupi_signal_file = new TFile("nutaupi_signal.root","RECREATE");
   nutaupi_signal_tree = new TTree("nutaupi_signal","Event data");
   create_sel_tree(nutaupi_signal_tree);
   nutaupi_bkg_file = new TFile("nutaupi_bkg.root","RECREATE");
   nutaupi_bkg_tree = new TTree("nutaupi_bkg","Event data");
   create_sel_tree(nutaupi_bkg_tree);

   // tau -> rho
   nutaurho_signal_file = new TFile("nutaurho_signal.root","RECREATE");
   nutaurho_signal_tree = new TTree("nutaurho_signal","Event data");
   create_sel_tree(nutaurho_signal_tree);
   nutaurho_bkg_file = new TFile("nutaurho_bkg.root","RECREATE");
   nutaurho_bkg_tree = new TTree("nutaurho_bkg","Event data");
   create_sel_tree(nutaurho_bkg_tree);
}

void fill_histos() {

  if (std::isnan(sumevent.Evis) || std::isnan(sumevent.ptmiss) || std::isnan(sumevent.cost) || std::isnan(sumevent.cosf)) {
    std::cout << " Nan found evis:" << sumevent.Evis << " ptmiss:"<<sumevent.ptmiss<<" cost:"<<sumevent.cost << " cosf:"<<sumevent.cosf<<std::endl;																		  
    std::cout << "jet " << sumevent.jetpx << " " << sumevent.jetpy << std::endl;
    std::cout << "tauvis " << sumevent.tauvis_px << " " << sumevent.tauvis_py << std::endl;
    std::cout << "istau " << sumevent.in_lepton_pdgid << std::endl;
    return;
  }
  
  //
  // tau->e channel
  //
  // BACKGROUND
  if(abs(sumevent.in_lepton_pdgid) == 12 && sumevent.isCC) {
    nueCC_bkg_tree->Fill();
  }
  // signal
  if(sumevent.istau && sumevent.isCC && sumevent.tau_decaymode == 1) {
    nueCC_signal_tree->Fill();
  }
  
  //
  // tau->mu channel
  //
  // BACKGROUND
  if(abs(sumevent.in_lepton_pdgid) == 14 && sumevent.isCC) {
    numuCC_bkg_tree->Fill();
  }
  // signal
  if(sumevent.istau && sumevent.isCC && sumevent.tau_decaymode == 2) {
    numuCC_signal_tree->Fill();
  }
  
  //
  // tau->pi channel
  //
  // BACKGROUND
  if(!sumevent.isCC && !sumevent.isES) {
    double taupi_cand_p = compute_momentum(sumevent.taupi_cand);
    if(taupi_cand_p>0) {
      float jet[3];
      jet[0] = sumevent.jetpx - sumevent.taupi_cand[0];
      jet[1] = sumevent.jetpy - sumevent.taupi_cand[1];
      jet[2] = sumevent.jetpz - sumevent.taupi_cand[2];
      sumevent.plep = taupi_cand_p;
      sumevent.ptlep = compute_pt(sumevent.taupi_cand);
      sumevent.qtlep = compute_qtlep(sumevent.taupi_cand, jet);
      sumevent.cost = compute_cost(sumevent.taupi_cand, jet);
      float ptmiss[3];
      ptmiss[0] = -sumevent.jetpx;
      ptmiss[1] = -sumevent.jetpy;
      sumevent.cosf = compute_cost(sumevent.taupi_cand, ptmiss);
      nutaupi_bkg_tree->Fill();
    }
  }
  // signal
  if(sumevent.istau && sumevent.isCC && sumevent.tau_decaymode == 3) {
    nutaupi_signal_tree->Fill();
  }
  
  //
  // tau->rho channel
  //
  // BACKGROUND
  if(!sumevent.isCC && !sumevent.isES) {
    double taurho_cand_p = compute_momentum(sumevent.taurho_cand);
    if(taurho_cand_p>0) {
      float jet[3];
      jet[0] = sumevent.jetpx - sumevent.taurho_cand[0];
      jet[1] = sumevent.jetpy - sumevent.taurho_cand[1];
      jet[2] = sumevent.jetpz - sumevent.taurho_cand[2];
      sumevent.plep = taurho_cand_p;
      sumevent.ptlep = compute_pt(sumevent.taurho_cand);
      sumevent.qtlep = compute_qtlep(sumevent.taurho_cand, jet);
      sumevent.cost = compute_cost(sumevent.taurho_cand, jet);
      float ptmiss[3];
      ptmiss[0] = -sumevent.jetpx;
      ptmiss[1] = -sumevent.jetpy;
      sumevent.cosf = compute_cost(sumevent.taurho_cand, ptmiss);
      nutaurho_bkg_tree->Fill();
    }
  }
  // signal
  if(sumevent.istau && sumevent.isCC && sumevent.tau_decaymode == 4) {
    nutaurho_signal_tree->Fill();
  }

}

void close_tuples(){
  std::cout << "Closing tuples..." << std::endl;
  nueCC_signal_file->cd();
  nueCC_signal_tree->Write();
  nueCC_signal_file->Close();

  nueCC_bkg_file->cd();
  nueCC_bkg_tree->Write();
  nueCC_bkg_file->Close();

  numuCC_signal_file->cd();
  numuCC_signal_tree->Write();
  numuCC_signal_file->Close();

  numuCC_bkg_file->cd();
  numuCC_bkg_tree->Write();
  numuCC_bkg_file->Close();

  nutaupi_signal_file->cd();
  nutaupi_signal_tree->Write();
  nutaupi_signal_file->Close();

  nutaupi_bkg_file->cd();
  nutaupi_bkg_tree->Write();
  nutaupi_bkg_file->Close();

  nutaurho_signal_file->cd();
  nutaurho_signal_tree->Write();
  nutaurho_signal_file->Close();

  nutaurho_bkg_file->cd();
  nutaurho_bkg_tree->Write();
  nutaurho_bkg_file->Close();
}

void smear_p(int pdgid, double px, double py, double pz, double *spx, double *spy, double *spz) {

  double m = 0;
  double p2 = px * px + py * py + pz * pz;
  double p = sqrt(p2);
  double e = sqrt(p2+m*m);

  // default smearing factors
  double eres_stoch = 0.5;
  double eres_const = 0.2;
  double angle_res = 1.0; // in degrees

  // Define smearing factors based on particle ID
  switch (abs(pdgid)) {
  case 11: // Electron
  case 22: // Gamma
  case 111: // pi0
    eres_stoch = 0.10;
    eres_const = 0.05;
    break;
  case 13: // Muon
    eres_stoch = 0;
    eres_const = 0.2;
    break;
  default:
    // Keep default values
    break;
  }
  
  // Apply Gaussian smearing to energy
  double eneres = sqrt(pow(eres_stoch/sqrt(e),2)+pow(eres_const,2));
  double randomValue = randomGen.Gaus(0.0, eneres);
  double smeared_e = e*(1.0+randomValue);

  // Angular smearing
  double degree_to_radian = M_PI / 180.0;
  double angle_stddev = angle_res * degree_to_radian; // 1 degree in radians
  double delta_theta = randomGen.Gaus(0.0, angle_stddev);
  double delta_phi = randomGen.Gaus(0.0, angle_stddev);
  
  // Current direction in spherical coordinates
  
  double theta = acos(pz / p);
  double phi = atan2(py, px);
  
  // Smear angles
  theta += delta_theta;
  phi += delta_phi;

  double smeared_p = smeared_e*smeared_e - m*m;
  if(smeared_p>0) {
    smeared_p = sqrt(smeared_p);
  } else {
    smeared_p = 0.1;
  }
  
  // Convert back to Cartesian coordinates
  *spx = smeared_p * sin(theta) * cos(phi);
  *spy = smeared_p * sin(theta) * sin(phi);
  *spz = smeared_p * cos(theta);
}

void smear_event() {
  double lep_px = sumevent.vis_spx-sumevent.jetpx;
  double lep_py = sumevent.vis_spy-sumevent.jetpy;
  double lep_pz = sumevent.vis_spy-sumevent.jetpy;

  // smear lepton e, or mu
  double slep_px = lep_px, slep_py = lep_py, slep_pz = lep_pz;
  if(abs(sumevent.in_lepton_pdgid) == 12) {
    smear_p(11, lep_px, lep_py, lep_pz, &slep_px, &slep_py, &slep_pz);
    lep_px = slep_px;
    lep_py = slep_py;
    lep_pz = slep_pz;
  }
  if(abs(sumevent.in_lepton_pdgid) == 14) {
    smear_p(13, lep_px, lep_py, lep_pz, &slep_px, &slep_py, &slep_pz);
    lep_px = slep_px;
    lep_py = slep_py;
    lep_pz = slep_pz;
  }

  // smear jet
  double sjet_px, sjet_py, sjet_pz;
  smear_p(211, sumevent.jetpx, sumevent.jetpy, sumevent.jetpz, &sjet_px, &sjet_py, &sjet_pz);
  sumevent.jetpx = sjet_px;
  sumevent.jetpy = sjet_py;
  sumevent.jetpz = sjet_pz;

  // recompute smeared total energy
  sumevent.vis_spx = lep_px+sumevent.jetpx;
  sumevent.vis_spy = lep_py+sumevent.jetpy;
  sumevent.vis_spz = lep_pz+sumevent.jetpz;
  
  //  sumevent.Evis = sqrt(sumevent.vis_spx*sumevent.vis_spx + sumevent.vis_spy*sumevent.vis_spy + sumevent.vis_spz*sumevent.vis_spz);
  sumevent.ptmiss = sqrt(sumevent.vis_spx*sumevent.vis_spx + sumevent.vis_spy*sumevent.vis_spy);
}

void extra_kinematics() {
  double lep_px, lep_py, lep_pz;
  double ptmissx, ptmissy;

  if(abs(sumevent.in_lepton_pdgid) == 16) {
    lep_px = sumevent.tauvis_px;
    lep_py = sumevent.tauvis_py;
    lep_pz = sumevent.tauvis_pz;
    ptmissx = -(sumevent.tauvis_px+sumevent.jetpx);
    ptmissy = -(sumevent.tauvis_py+sumevent.jetpy);
  } else {
    lep_px = sumevent.vis_spx-sumevent.jetpx;
    lep_py = sumevent.vis_spy-sumevent.jetpy;
    lep_pz = sumevent.vis_spz-sumevent.jetpz;
    ptmissx = -sumevent.vis_spx;
    ptmissy = -sumevent.vis_spy;
  }
  sumevent.cost = (lep_px*sumevent.jetpx + lep_py*sumevent.jetpy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(sumevent.jetpx*sumevent.jetpx+sumevent.jetpy*sumevent.jetpy));
  sumevent.cosf = (lep_px*ptmissx + lep_py*ptmissy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(ptmissx*ptmissx+ptmissy*ptmissy));

  sumevent.plep = sqrt(lep_px*lep_px+lep_py*lep_py+lep_pz*lep_pz);
  sumevent.ptlep = sqrt(lep_px*lep_px+lep_py*lep_py);
  // projection of lepton momentum on hadronic jet
  double lepproj = (lep_px*sumevent.jetpx + lep_py*sumevent.jetpy + lep_pz*sumevent.jetpz)/sqrt(sumevent.jetpx*sumevent.jetpx+sumevent.jetpy*sumevent.jetpy+sumevent.jetpz*sumevent.jetpz);
  double lep2 = lep_px*lep_px+lep_py*lep_py+lep_pz*lep_pz;
  sumevent.qtlep = sqrt(lep2-lepproj*lepproj);

  if (std::isnan(sumevent.cost) || std::isnan(sumevent.cosf) || std::isnan(sumevent.ptlep) || std::isnan(sumevent.qtlep)) {
    std::cout << "kinematics nan" << std::endl;
    std::cout << "ptlep " << sumevent.ptlep << std::endl;
    std::cout << "qtlep " << sumevent.qtlep << std::endl;
    std::cout << "lepproj " << lepproj << std::endl;
    std::cout << "lep2 " << lep2 << std::endl;
    std::cout << "lep " << lep_px << " " << lep_py << std::endl;
    std::cout << "jet " << sumevent.jetpx << " " << sumevent.jetpy << std::endl;
    std::cout << "ptmiss " << ptmissx << " " << ptmissy << std::endl;
    std::cout << "istau " << sumevent.in_lepton_pdgid << std::endl;
  }
}

void s(){

  stats.infid = stats.nueCC = stats.numuCC = stats.nutauCC = stats.NC = stats.ES = 0;
  
  open_tuples();

   Long64_t nentries = event_tree->GetEntries();

   for (Long64_t i=0; i<nentries;i++) {
     event_tree->GetEntry(i);

     // check vtx fiducial

     sumevent.infid = sumevent.prim_vz > -3000 && sumevent.prim_vz<-1900 && abs(sumevent.prim_vx)<140 && abs(sumevent.prim_vy)<200;
     //     if(!infid) continue;

     stats.infid++;
     
     // stats
     double jetp = sqrt(sumevent.jetpx*sumevent.jetpx+sumevent.jetpy*sumevent.jetpy+sumevent.jetpz*sumevent.jetpz);
     double tauvis = sqrt(sumevent.tauvis_px*sumevent.tauvis_px+sumevent.tauvis_py*sumevent.tauvis_py+sumevent.tauvis_pz*sumevent.tauvis_pz);

     // check for ELASTIC scattering on electrons
     if(sumevent.isES) {
       stats.ES++;
       continue;
     } else {
       if(sumevent.isCC) {
	 if(abs(sumevent.in_lepton_pdgid) == 12) {
	   stats.nueCC++;
	 }
	 if(abs(sumevent.in_lepton_pdgid) == 14) {
	   stats.numuCC++;
	 }
	 if(abs(sumevent.in_lepton_pdgid) == 16) {
	   if(tauvis > 0) {
	     stats.nutauCC++;
	   } else {
	     stats.nutauCClost++;
	     continue;
	   }
	 }
       } else {
	 stats.NC++;
	 //	 continue;
       };
     };
     if(sumevent.isCC) {
       smear_event();
       extra_kinematics();
     };
     fill_histos();
   }

   std::cout << "in fiducial = " << stats.infid;
   std::cout << " nueCC = " << stats.nueCC;
   std::cout << " numuCC = " << stats.numuCC;
   std::cout << " nutauCC = " << stats.nutauCC;
   std::cout << " nutauCClost = " << stats.nutauCClost;
   std::cout << " NC = " << stats.NC << std::endl;
   std::cout << " ES = " << stats.ES << std::endl;

   close_tuples();
   std::cout << "Done!" << std::endl;
}
