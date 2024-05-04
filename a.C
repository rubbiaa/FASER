void a() {


  TChain *event_tree = new TChain("event_tree");
  
  event_tree->Add("event_data.root");

  std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

  //   tree->MakeCode("skel.h");

  //Declaration of leaves types
   Bool_t          isCC;
   Bool_t          istau;
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
   Double_t        Evis;
   Double_t        ptmiss;

   // Set branch addresses.
   event_tree->SetBranchAddress("isCC",&isCC);
   event_tree->SetBranchAddress("istau",&istau);
   event_tree->SetBranchAddress("tau_decaymode",&tau_decaymode);
   event_tree->SetBranchAddress("n_particles",&n_particles);
   event_tree->SetBranchAddress("in_lepton_pdgid",&in_lepton_pdgid);
   event_tree->SetBranchAddress("vis_spx",&vis_spx);
   event_tree->SetBranchAddress("vis_spy",&vis_spy);
   event_tree->SetBranchAddress("vis_spz",&vis_spz);
   event_tree->SetBranchAddress("jetpx",&jetpx);
   event_tree->SetBranchAddress("jetpy",&jetpy);
   event_tree->SetBranchAddress("jetpz",&jetpz);
   event_tree->SetBranchAddress("tauvis_px",&tauvis_px);
   event_tree->SetBranchAddress("tauvis_py",&tauvis_py);
   event_tree->SetBranchAddress("tauvis_pz",&tauvis_pz);
   event_tree->SetBranchAddress("Evis",&Evis);
   event_tree->SetBranchAddress("ptmiss",&ptmiss);

   // histograms
   TH1D* nueCC_Evis = new TH1D("nueCC_Evis","Evis nueCC",100,0,2000);
   TH1D* nueCC_ptmiss = new TH1D("nueCC_ptmiss","ptmiss nueCC",100,0,20);

   TH1D* nutaueCC_Evis = new TH1D("nutaueCC_Evis","Evis nutau->e CC",100,0,2000);
   TH1D* nutaueCC_ptmiss = new TH1D("nutaueCC_ptmiss","ptmiss nutau->e CC",100,0,20);
   
   Long64_t nentries = event_tree->GetEntries();

   for (Long64_t i=0; i<nentries;i++) {
     event_tree->GetEntry(i);

     //
     // tau->e channel
     //
     // BACKGROUND
     if(abs(in_lepton_pdgid) == 12) {
       nueCC_Evis->Fill(Evis);
       nueCC_ptmiss->Fill(ptmiss);
     }
     // signal
     if(abs(in_lepton_pdgid) == 16 && tau_decaymode == 1) {
       nutaueCC_Evis->Fill(Evis);
       nutaueCC_ptmiss->Fill(ptmiss);       
     }
     
   }
}
