// create signal and background summaries
// A. Rubbia/June 2024
//

// summary of event
struct SUMEVENT {
   Bool_t          isCC;
   Bool_t          istau;
   Double_t        prim_vx;
   Double_t        prim_vy;
   Double_t        prim_vz;
   Int_t           tau_decaymode; // =1 e, =2 mu, =3 1-prong, =4 rho =5 3-prong, =6 other
   ULong_t         n_particles;
   Int_t           in_lepton_pdgid;
   Double_t        vis_spx;
   Double_t        vis_spy;
   Double_t        vis_spz;
   Double_t        jetpx;
   Double_t        jetpy;
   Double_t        jetpz;
   Double_t        tauvis_px;
   Double_t        tauvis_py;
   Double_t        tauvis_pz;
   Double_t        tautracklength;
   Double_t        Evis;
   Double_t        ptmiss;
} sumevent;

struct stats {
  int nueCC;
  int numuCC;
  int nutauCC;
  int NC;
} stats;

void open_tubles() {
  TChain *event_tree = new TChain("event_tree");
  
  event_tree->Add("event_data.root");
  
  std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

     // Set branch addresses.

   event_tree->SetBranchAddress("isCC",&sumevent.isCC);
   event_tree->SetBranchAddress("istau",&sumevent.istau);
   event_tree->SetBranchAddress("tau_decaymode",&sumevent.tau_decaymode);
   event_tree->SetBranchAddress("n_particles",&sumevent.n_particles);
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
}


void s(){
  open_tuples();

  
}
