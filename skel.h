{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Sat May  4 01:51:06 2024 by ROOT version6.30/04)
//   from TChain event_tree/
//////////////////////////////////////////////////////////


//Reset ROOT and connect tree file
   gROOT->Reset();

#ifdef SINGLE_TREE
   // The following code should be used if you want this code to access
   // a single tree instead of a chain
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("event_data.root");
   if (!f) {
      f = new TFile("event_data.root");
   }
    f->GetObject("event_tree",tree);

#else // SINGLE_TREE

   // The following code should be used if you want this code to access a chain
   // of trees.
   TChain *event_tree = new TChain("event_tree","");
   event_tree->Add("event_data.root/event_tree");
#endif // SINGLE_TREE

//Declaration of leaves types
   Bool_t          isCC;
   Bool_t          istau;
   ULong_t         n_particles;
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
   event_tree->SetBranchAddress("n_particles",&n_particles);
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

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
// event_tree->SetBranchStatus("*",0);  // disable all branches
// TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

   Long64_t nentries = event_tree->GetEntries();

   Long64_t nbytes = 0;
//   for (Long64_t i=0; i<nentries;i++) {
//      nbytes += event_tree->GetEntry(i);
//   }
}
