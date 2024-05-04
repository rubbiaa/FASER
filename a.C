#include <TPad.h>

struct cut {
  const char *name;
  int sig_passing, back_passing;
};

#define MAXCUTS 10
struct cuts {
  const char *name;
  int ncuts;
  struct cut cut[MAXCUTS];
};
    

void dump_cuts(struct cuts *cuts) {
  std::cout << "------ " << cuts->name << " -------------------------------------------------------------------" << std::endl;
  for (int ic = 0; ic<cuts->ncuts;ic++) {
    double rsig =  cuts->cut[ic].sig_passing*100.0/ cuts->cut[0].sig_passing;
    double rback =  cuts->cut[ic].back_passing*100.0/ cuts->cut[0].back_passing;
    std::cout << std::setw(10) << ic << " " << cuts->cut[ic].name << " sig=" << cuts->cut[ic].sig_passing << "(" << rsig << "%) ";
    std::cout << std::setw(10) << " bkng= " << cuts->cut[ic].back_passing << "(" << rback << "%) ";
    std::cout << std::endl;
  }
  std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
}

struct stats {
  int nueCC;
  int numuCC;
  int nutauCC;
  int NC;
} stats;

void a() {

  stats.nueCC = stats.numuCC = stats.nutauCC = stats.NC = 0;
  
  struct cuts tauecuts;
  struct cuts taumucuts;

  tauecuts.ncuts = 3;
  tauecuts.name = "taue analysis";
  struct cut tauecut_1 = {"All",0,0};
  struct cut tauecut_2 = {"Evis<2000",0,0};
  struct cut tauecut_3 = {"ptmiss>2",0,0};
  tauecuts.cut[0] = tauecut_1;
  tauecuts.cut[1] = tauecut_2;
  tauecuts.cut[2] = tauecut_3;

  taumucuts.ncuts = 4;
  taumucuts.name = "taumu analysis";
  struct cut taumucut_1 = {"All",0,0};
  struct cut taumucut_2 = {"Evis<2000",0,0};
  struct cut taumucut_3 = {"ptmiss>2",0,0};
  struct cut taumucut_4 = {"cost-cosf cut",0,0};
  taumucuts.cut[0] = taumucut_1;
  taumucuts.cut[1] = taumucut_2;
  taumucuts.cut[2] = taumucut_3;
  taumucuts.cut[3] = taumucut_4;

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
   TH1D* nueCC_Evis = new TH1D("nueCC_Evis","Evis nueCC",100,0,4000);
   TH1D* nueCC_ptmiss = new TH1D("nueCC_ptmiss","ptmiss nueCC",100,0,20);
   TH2D* nueCC_costcosf = new TH2D("nueCC_costcosf", "angles nueCC", 100, -1., 1., 100, -1, 1.);
     
   TH1D* nutaueCC_Evis = new TH1D("nutaueCC_Evis","Evis tau->e CC",100,0,4000);
   TH1D* nutaueCC_ptmiss = new TH1D("nutaueCC_ptmiss","ptmiss tau->e CC",100,0,20);
   TH2D* nutaueCC_costcosf = new TH2D("nutaueCC_costcosf", "angles tau->e CC", 100, -1., 1., 100, -1, 1.);
   
   TH1D* numuCC_Evis = new TH1D("numuCC_Evis","Evis numuCC",100,0,4000);
   TH1D* numuCC_ptmiss = new TH1D("numuCC_ptmiss","ptmiss numuCC",100,0,20);
   TH2D* numuCC_costcosf = new TH2D("numuCC_costcosf", "angles numuCC", 100, -1., 1., 100, -1, 1.);
   
   TH1D* nutaumuCC_Evis = new TH1D("nutaumuCC_Evis","Evis tau->mu CC",100,0,4000);
   TH1D* nutaumuCC_ptmiss = new TH1D("nutaumuCC_ptmiss","ptmiss tau->mu CC",100,0,20);
   TH2D* nutaumuCC_costcosf = new TH2D("nutaumuCC_costcosf", "angles tau-> numu CC", 100, -1., 1., 100, -1, 1.);
   
   Long64_t nentries = event_tree->GetEntries();

   for (Long64_t i=0; i<nentries;i++) {
     event_tree->GetEntry(i);

     // stats
     if(isCC) {
       if(abs(in_lepton_pdgid) == 12) {
	 stats.nueCC++;
       }
       if(abs(in_lepton_pdgid) == 14) {
	 stats.numuCC++;
       }
       if(abs(in_lepton_pdgid) == 16) {
	 stats.nutauCC++;
       }
     } else {
       stats.NC++;
     };
     
     double cost, cosf;
     // extra kinematics
     if(abs(in_lepton_pdgid) == 16) {
       cost = (tauvis_px*jetpx + tauvis_py*jetpy)/(sqrt(tauvis_px*tauvis_px+tauvis_py*tauvis_py)*sqrt(jetpx*jetpx+jetpy*jetpy));
       double ptmissx = -(tauvis_px+jetpx);
       double ptmissy = -(tauvis_py+jetpy);
       cosf = (tauvis_px*ptmissx + tauvis_py*ptmissy)/(sqrt(tauvis_px*tauvis_px+tauvis_py*tauvis_py)*sqrt(ptmissx*ptmissx+ptmissy*ptmissy));
     } else {
       double lep_px = vis_spx-jetpx;
       double lep_py = vis_spy-jetpy;
       cost = (lep_px*jetpx + lep_py*jetpy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(jetpx*jetpx+jetpy*jetpy));
       double ptmissx = -vis_spx;
       double ptmissy = -vis_spy;
       cosf = (lep_px*ptmissx + lep_py*ptmissy)/(sqrt(lep_px*lep_px+lep_py*lep_py)*sqrt(ptmissx*ptmissx+ptmissy*ptmissy));
     }

     //
     // tau->e channel
     //
     // BACKGROUND
     if(abs(in_lepton_pdgid) == 12 && isCC) {
       nueCC_Evis->Fill(Evis);
       nueCC_ptmiss->Fill(ptmiss);
       nueCC_costcosf->Fill(cost,cosf);
       tauecuts.cut[0].back_passing++;
       if(Evis < 2000){
	 tauecuts.cut[1].back_passing++;
	 if(ptmiss > 2){
	   tauecuts.cut[2].back_passing++;
	 }
       }
     }
     // signal
     if(istau && isCC && tau_decaymode == 1) {
       nutaueCC_Evis->Fill(Evis);
       nutaueCC_ptmiss->Fill(ptmiss);
       nutaueCC_costcosf->Fill(cost,cosf);
       tauecuts.cut[0].sig_passing++;
       if(Evis < 2000){
	 tauecuts.cut[1].sig_passing++;
	 if(ptmiss > 2){
	   tauecuts.cut[2].sig_passing++;
	 }
       }
     }

     //
     // tau->mu channel
     //
     // BACKGROUND
     if(abs(in_lepton_pdgid) == 14 && isCC) {
       numuCC_Evis->Fill(Evis);
       numuCC_ptmiss->Fill(ptmiss);
       numuCC_costcosf->Fill(cost,cosf);
       taumucuts.cut[0].back_passing++;
       if(Evis < 2000){
	 taumucuts.cut[1].back_passing++;
	 if(ptmiss > 2){
	   taumucuts.cut[2].back_passing++;
	   if(cost<-0.8 && cosf>0.8){
	     taumucuts.cut[3].back_passing++;
	   }
	 }
       }
     }
     // signal
     if(istau && isCC && tau_decaymode == 2) {
       nutaumuCC_Evis->Fill(Evis);
       nutaumuCC_ptmiss->Fill(ptmiss);
       nutaumuCC_costcosf->Fill(cost,cosf);
       taumucuts.cut[0].sig_passing++;
       if(Evis < 2000){
	 taumucuts.cut[1].sig_passing++;
	 if(ptmiss > 2){
	   taumucuts.cut[2].sig_passing++;
	   if(cost<-0.8 && cosf>0.8){
	     taumucuts.cut[3].sig_passing++;
	   }
	 }
       }
     }
     
   }
   
   TCanvas *c1 = new TCanvas("c1", "electron Evis", 800, 600);
   c1->Divide(1, 2);
   c1->cd(1);
   nueCC_Evis->Draw();
   c1->cd(2);
   nutaueCC_Evis->Draw();
   c1->SaveAs("c1.png");

   TCanvas *c2 = new TCanvas("c2", "electron ptmiss", 800, 600);
   c2->Divide(1, 2);
   c2->cd(1);
   gPad->SetLogy();
   nueCC_ptmiss->Draw();
   c2->cd(2);
   nutaueCC_ptmiss->Draw();
   c2->SaveAs("c2.png");

   TCanvas *c3 = new TCanvas("c3", "electron angles", 800, 600);
   c3->Divide(1, 2);
   c3->cd(1);
   nueCC_costcosf->Draw();
   c3->cd(2);
   nutaueCC_costcosf->Draw();
   c3->SaveAs("c3.png");

   TCanvas *c1m = new TCanvas("c1m", "muon Evis", 800, 600);
   c1m->Divide(1, 2);
   c1m->cd(1);
   numuCC_Evis->Draw();
   c1m->cd(2);
   nutaumuCC_Evis->Draw();
   c1m->SaveAs("c1m.png");

   TCanvas *c2m = new TCanvas("c2m", "tau ptmiss", 800, 600);
   c2m->Divide(1, 2);
   c2m->cd(1);
   gPad->SetLogy();
   numuCC_ptmiss->Draw();
   c2m->cd(2);
   nutaumuCC_ptmiss->Draw();
   c2m->SaveAs("c2m.png");

   TCanvas *c3m = new TCanvas("c3m", "muon angles", 800, 600);
   c3m->Divide(1, 2);
   c3m->cd(1);
   numuCC_costcosf->Draw();
   c3m->cd(2);
   nutaumuCC_costcosf->Draw();
   c3m->SaveAs("c3m.png");
   
   dump_cuts(&tauecuts);
   dump_cuts(&taumucuts);

   std::cout << "nueCC = " << stats.nueCC;
   std::cout << " numuCC = " << stats.numuCC;
   std::cout << " nutauCC = " << stats.nutauCC;
   std::cout << " NC = " << stats.NC << std::endl;

 }
