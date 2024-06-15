#include <TMVA/Factory.h>
#include <TMVA/DataLoader.h>
#include <TMVA/Tools.h>
#include <TMVA/MethodBDT.h>
#include <TMVA/Reader.h>
#include <TCut.h>

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


void trainBDT(const char *dataset, const char *tmvaoutput, const char *sigtuple, const char *sigtree, const char *bkgtuple, const char *bkgtree) {
    // Initialize TMVA
    TMVA::Tools::Instance();

    // Create a new ROOT output file
    TFile* outputFile = TFile::Open(tmvaoutput, "RECREATE");

    // Create the factory and dataloader
    TMVA::Factory factory("TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D");
    TMVA::DataLoader dataloader(dataset);

    // Load your signal and background files
    TFile signalFile(sigtuple);
    TTree* signalTree = (TTree*)signalFile.Get(sigtree);

    TFile backgroundFile(bkgtuple);
    TTree* backgroundTree = (TTree*)backgroundFile.Get(bkgtree);

    dataloader.AddVariable("evis", 'F');
    dataloader.AddVariable("ptmiss", 'F');
    dataloader.AddVariable("cost", 'F');
    dataloader.AddVariable("cosf", 'F');
    dataloader.AddVariable("plep", 'F');
    dataloader.AddVariable("ptlep", 'F');
    dataloader.AddVariable("qtlep", 'F');

    // Add signal and background trees
    dataloader.AddSignalTree(signalTree, 1.0);
    dataloader.AddBackgroundTree(backgroundTree, 1.0);

    // Apply cuts if needed
    TCut signalCut = "";     // Define any signal-specific cuts
    TCut backgroundCut = ""; // Define any background-specific cuts

    // Prepare training and test trees
    dataloader.PrepareTrainingAndTestTree(signalCut, backgroundCut, "nTrain_Signal=1000:nTrain_Background=10000:nTest_Signal=500:nTest_Background=5000:SplitMode=Random:NormMode=NumEvents:!V");

    // Set BDT parameters
    factory.BookMethod(&dataloader, TMVA::Types::kBDT, "BDT",
                       "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");

    // Train the BDT
    factory.TrainAllMethods();

    // Test the BDT
    factory.TestAllMethods();

    // Evaluate the BDT
    factory.EvaluateAllMethods();
    
    outputFile->Close();
}

// Plot the performance of the BDT
void PlotBDTPerformance(const char *tmvaoutput) {
    // Launch the TMVAGui
    if (!gSystem->AccessPathName(tmvaoutput)) {
        TMVA::TMVAGui(tmvaoutput);
    } else {
        std::cerr << "Error: Could not find TMVAOutput.root" << std::endl;
    }
}

// histograms
TH1D* nueCC_Evis = new TH1D("nueCC_Evis","Evis nueCC",100,0,4000);
TH1D* nueCC_ptmiss = new TH1D("nueCC_ptmiss","ptmiss nueCC",100,0,20);
TH2D* nueCC_costcosf = new TH2D("nueCC_costcosf", "angles nueCC", 100, -1., 1., 100, -1, 1.);

// const int nBins = 6;
//Double_t binEdges[nBins+1] = {-0.3, -0.1, 0.0, 0.1, 0.2, 0.4, 0.5};

#if 0
const int nBins = 83;
Double_t binEdges[] = {
  -0.3000, -0.2938, -0.2875, -0.2813, -0.2750, -0.2688, -0.2625, -0.2563, -0.2500, -0.2438,
  -0.2375, -0.2313, -0.2250, -0.2188, -0.2125, -0.2063, -0.2000, -0.1938, -0.1875, -0.1813,
  -0.1750, -0.1688, -0.1625, -0.1563, -0.1500, -0.1438, -0.1375, -0.1313, -0.1250, -0.1188,
  -0.1125, -0.1063, -0.1000, -0.0938, -0.0875, -0.0813, -0.0750, -0.0688, -0.0625, -0.0563,
  -0.0500, -0.0438, -0.0375, -0.0313, -0.0250, -0.0188, -0.0125, -0.0063, 0.0000,  0.0063,
  0.0125,  0.0188,  0.0250,  0.0313,  0.0375,  0.0438,  0.0500,  0.0563,  0.0625,  0.0688,
  0.0750,  0.0813,  0.0875,  0.0938,  0.1000,  0.1063,  0.1125,  0.1188,  0.1250,  0.1313,
  0.1375,  0.1438,  0.1500,  0.1563,  0.1625,  0.1688,  0.1750,  0.1813,  0.1875,  0.1938,
  0.2000,
  // Coarse binning from 0.2 to 0.5
  0.3000,  0.4000,  0.5000
};
#endif

// TH1D* nueCC_bdt = new TH1D("nueCC_bdt","bdt nueCC",nBins,binEdges);
TH1D* nueCC_bdt = new TH1D("nueCC_bdt","bdt tau->e", 100,-0.3,0.5);

TH1D* nutaueCC_Evis = new TH1D("nutaueCC_Evis","Evis tau->e CC",100,0,4000);
TH1D* nutaueCC_ptmiss = new TH1D("nutaueCC_ptmiss","ptmiss tau->e CC",100,0,20);
TH2D* nutaueCC_costcosf = new TH2D("nutaueCC_costcosf", "angles tau->e CC", 100, -1., 1., 100, -1, 1.);
// TH1D* nutaueCC_bdt = new TH1D("nutaueCC_bdt","bdt nutaueCC",nBins,binEdges);
TH1D* nutaueCC_bdt = new TH1D("nutaueCC_bdt","bdt nutaueCC",100,-0.3,0.5);

TH1D* numuCC_Evis = new TH1D("numuCC_Evis","Evis numuCC",100,0,4000);
TH1D* numuCC_ptmiss = new TH1D("numuCC_ptmiss","ptmiss numuCC",100,0,20);
TH2D* numuCC_costcosf = new TH2D("numuCC_costcosf", "angles numuCC", 100, -1., 1., 100, -1, 1.);
TH1D* numuCC_bdt = new TH1D("numuCC_bdt","bdt tau->mu", 100,-0.3,0.5);

TH1D* nutaumuCC_Evis = new TH1D("nutaumuCC_Evis","Evis tau->mu CC",100,0,4000);
TH1D* nutaumuCC_ptmiss = new TH1D("nutaumuCC_ptmiss","ptmiss tau->mu CC",100,0,20);
TH2D* nutaumuCC_costcosf = new TH2D("nutaumuCC_costcosf", "angles tau-> numu CC", 100, -1., 1., 100, -1, 1.);
TH1D* nutaumuCC_bdt = new TH1D("nutaumuCC_bdt","bdt nutaumuCC", 100,-0.3,0.5);

TH1D* nutaupi_sig_bdt = new TH1D("nutaupi_sig_bdt","bdt tau->pi", 100,-0.3,0.5);
TH1D* nutaupi_bkg_bdt = new TH1D("nutaupi_bkg_bdt","bdt tau->pi", 100,-0.3,0.5);

TH1D* nutaurho_sig_bdt = new TH1D("nutaurho_sig_bdt","bdt tau->pi", 100,-0.3,0.5);
TH1D* nutaurho_bkg_bdt = new TH1D("nutaurho_bkg_bdt","bdt tau->pi", 100,-0.3,0.5);

void dump_cuts(struct cuts *cuts) {
  std::cout << "------ " << cuts->name << " -------------------------------------------------------------------" << std::endl;
  for (int ic = 0; ic<cuts->ncuts;ic++) {
    double rsig =  cuts->cut[ic].sig_passing*100.0/ cuts->cut[0].sig_passing;
    double rback =  cuts->cut[ic].back_passing*100.0/ cuts->cut[0].back_passing;
    std::cout << std::setw(10) << ic << " " << cuts->cut[ic].name << " sig=" << cuts->cut[ic].sig_passing << "(" << rsig << "%) ";
    std::cout << std::setw(10) << " bkng= " << cuts->cut[ic].back_passing << "(" << rback << "%) ";
    if(cuts->cut[ic].back_passing>0) {
      std::cout << " S/√B = " << cuts->cut[ic].sig_passing/sqrt(cuts->cut[ic].back_passing);
    }
    std::cout << std::endl;
  }
  std::cout << "--------------------------------------------------------------------------------------------------------" << std::endl;
}

void rescale_cuts(struct cuts *cuts){
  std::cout << "Recaling " << cuts->name << std::endl;
  double factor = 0.68/10.0;
  for (int ic = 0; ic<cuts->ncuts;ic++) {
    cuts->cut[ic].sig_passing *= factor;
    cuts->cut[ic].back_passing *= factor;
  }
};
    
struct SELTREE {
  Bool_t   infid;
  Float_t  Evis;
  Float_t  ptmiss;
  Float_t  cost;
  Float_t  cosf;
  Float_t  plep;
  Float_t  ptlep;
  Float_t  qtlep;
} seltree;

void process_file(int sig, const char *dataset, const char *fname, const char *chain, const char *tmvaoutput, struct cuts *cuts, TH1D* evis, TH1D *ptmiss, TH2D *costcosf, TH1D *bdt) {
  TChain *event_tree = new TChain(chain);
  event_tree->Add(fname);
  std::cout << "Number of entries " << event_tree->GetEntries() << std::endl;

  event_tree->SetBranchAddress("infid",&seltree.infid);
  event_tree->SetBranchAddress("evis",&seltree.Evis);
  event_tree->SetBranchAddress("ptmiss",&seltree.ptmiss);
  event_tree->SetBranchAddress("cost",&seltree.cost);
  event_tree->SetBranchAddress("cosf",&seltree.cosf);
  event_tree->SetBranchAddress("plep",&seltree.plep);
  event_tree->SetBranchAddress("ptlep",&seltree.ptlep);
  event_tree->SetBranchAddress("qtlep",&seltree.qtlep);

  // Initialize TMVA Reader
  TMVA::Reader reader("!Color:!Silent");
  
  // Define the variables to be read
  reader.AddVariable("evis", &seltree.Evis);
  reader.AddVariable("ptmiss", &seltree.ptmiss);
  reader.AddVariable("cost", &seltree.cost);
  reader.AddVariable("cosf", &seltree.cosf);
  reader.AddVariable("plep", &seltree.plep);
  reader.AddVariable("ptlep", &seltree.ptlep);
  reader.AddVariable("qtlep", &seltree.qtlep);

  TFile* weightsFile = TFile::Open(tmvaoutput);
  // Get the path to the weights
  std::string s;
  s = string(dataset) + string("/weights/TMVAClassification_BDT.weights.xml");
  const char* weightsFilePath = s.c_str();
  
  // Book the BDT method using the weights file path
  reader.BookMVA("BDT", weightsFilePath);
  
  Long64_t nentries = event_tree->GetEntries();
  
  for (Long64_t i=0; i<nentries;i++) {
    event_tree->GetEntry(i);

    // keep only events in fiducial
    if(!seltree.infid) continue;

    if(evis != nullptr) {
      evis->Fill(seltree.Evis);
      ptmiss->Fill(seltree.ptmiss);
      costcosf->Fill(seltree.cost,seltree.cosf);
    }

    double bdtResponse = reader.EvaluateMVA("BDT");
    if(bdt != nullptr) {
      bdt->Fill(bdtResponse);
    };

    double cutval[] = {-999, 0, 0.15, 0.2, 0.25, 0.3};
    for (int ic=0; ic<6; ic++) {
      if(bdtResponse>cutval[ic]) {
	if(sig) {
	  cuts->cut[ic].sig_passing++;
	} else {
	  cuts->cut[ic].back_passing++;
	}
      }
    }
  }
  weightsFile->Close();
}

void a() {

  //  trainBDT("taueset", "nueCC_TMVAOutput.root", "nuecc_signal.root", "nuecc_signal", "nuecc_bkg.root", "nuecc_bkg");
  //  PlotBDTPerformance("nueCC_TMVAOutput.root");
  
  //    trainBDT("taumuset", "numuCC_TMVAOutput.root", "numucc_signal.root", "numucc_signal", "numucc_bkg.root", "numucc_bkg");
  //  PlotBDTPerformance("numuCC_TMVAOutput.root");

  //  trainBDT("taupiset", "nutaupi_TMVAOutput.root", "nutaupi_signal.root", "nutaupi_signal", "nutaupi_bkg.root", "nutaupi_bkg");
  // PlotBDTPerformance("nutaupi_TMVAOutput.root");
  
  //  trainBDT("taurhoset", "nutaurho_TMVAOutput.root", "nutaurho_signal.root", "nutaurho_signal", "nutaurho_bkg.root", "nutaurho_bkg");
  //  PlotBDTPerformance("nutaurho_TMVAOutput.root");

  struct cuts tauecuts;
  struct cuts taumucuts;
  struct cuts taupicuts;
  struct cuts taurhocuts;

  tauecuts.ncuts = 6;
  tauecuts.name = "taue analysis";
  struct cut tauecut_1 = {"All",0,0};
  struct cut tauecut_2 = {"bdt>0",0,0};
  struct cut tauecut_3 = {"bdt>0.15",0,0};
  struct cut tauecut_4 = {"bdt>0.2",0,0};
  struct cut tauecut_5 = {"bdt>0.25",0,0};
  struct cut tauecut_6 = {"bdt>0.3",0,0};
  tauecuts.cut[0] = tauecut_1;
  tauecuts.cut[1] = tauecut_2;
  tauecuts.cut[2] = tauecut_3;
  tauecuts.cut[3] = tauecut_4;
  tauecuts.cut[4] = tauecut_5;
  tauecuts.cut[5] = tauecut_6;

  taumucuts.ncuts = 6;
  taumucuts.name = "taumu analysis";
  taumucuts.cut[0] = tauecut_1;
  taumucuts.cut[1] = tauecut_2;
  taumucuts.cut[2] = tauecut_3;
  taumucuts.cut[3] = tauecut_4;
  taumucuts.cut[4] = tauecut_5;
  taumucuts.cut[5] = tauecut_6;

  taupicuts.ncuts = 6;
  taupicuts.name = "taupi analysis";
  taupicuts.cut[0] = tauecut_1;
  taupicuts.cut[1] = tauecut_2;
  taupicuts.cut[2] = tauecut_3;
  taupicuts.cut[3] = tauecut_4;
  taupicuts.cut[4] = tauecut_5;
  taupicuts.cut[5] = tauecut_6;

  taurhocuts.ncuts = 6;
  taurhocuts.name = "taurho analysis";
  taurhocuts.cut[0] = tauecut_1;
  taurhocuts.cut[1] = tauecut_2;
  taurhocuts.cut[2] = tauecut_3;
  taurhocuts.cut[3] = tauecut_4;
  taurhocuts.cut[4] = tauecut_5;
  taurhocuts.cut[5] = tauecut_6;

  process_file(0, "taueset", "nuecc_bkg.root","nuecc_bkg", "nueCC_TMVAOutput.root", &tauecuts, nueCC_Evis, nueCC_ptmiss, nueCC_costcosf, nueCC_bdt);
  process_file(1, "taueset", "nuecc_signal.root","nuecc_signal", "nueCC_TMVAOutput.root", &tauecuts, nutaueCC_Evis, nutaueCC_ptmiss, nutaueCC_costcosf, nutaueCC_bdt);

  process_file(0, "taumuset", "numucc_bkg.root","numucc_bkg", "numuCC_TMVAOutput.root", &taumucuts, numuCC_Evis, numuCC_ptmiss, numuCC_costcosf, numuCC_bdt);
  process_file(1, "taumuset", "numucc_signal.root","numucc_signal", "numuCC_TMVAOutput.root", &taumucuts, nutaumuCC_Evis, nutaumuCC_ptmiss, nutaumuCC_costcosf, nutaumuCC_bdt);

  process_file(0, "taupiset", "nutaupi_bkg.root","nutaupi_bkg", "nutaupi_TMVAOutput.root", &taupicuts, nullptr, nullptr, nullptr, nutaupi_bkg_bdt);
  process_file(1, "taupiset", "nutaupi_signal.root","nutaupi_signal", "nutaupi_TMVAOutput.root", &taupicuts, 
	       nullptr, nullptr, nullptr, nutaupi_sig_bdt);

  process_file(0, "taurhoset", "nutaurho_bkg.root","nutaurho_bkg", "nutaurho_TMVAOutput.root", &taurhocuts, nullptr, nullptr, nullptr, nutaurho_bkg_bdt);
  process_file(1, "taurhoset", "nutaurho_signal.root","nutaurho_signal", "nutaurho_TMVAOutput.root", &taurhocuts, 
	       nullptr, nullptr, nullptr, nutaurho_sig_bdt);
  
  #if 0
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
   //   gPad->SetLogy();
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
   #endif

   TCanvas *c4 = new TCanvas("c4", "bdt", 800, 600);
   gPad->SetLogy();
   nueCC_bdt->SetFillColor(kYellow);  
   nueCC_bdt->Draw("HIST");
   //   nutaueCC_bdt->Draw("same");
   TH1F* sum_bdt = (TH1F*)nueCC_bdt->Clone("sum_bdt");
   sum_bdt->Add(nutaueCC_bdt);
   sum_bdt->SetLineColor(kRed);  // Set color for the sum histogram
   sum_bdt->SetLineWidth(2);     // Set line width for better visibility
   sum_bdt->SetTitle("Sum of nueCC and nutaueCC BDT");
   sum_bdt->Draw("e SAME");
   c4->SaveAs("c4.png");

   #if 0
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
   //   gPad->SetLogy();
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
#endif

   TCanvas *c4pi = new TCanvas("c4pi", "bdt", 800, 600);
   gPad->SetLogy();
   nutaupi_bkg_bdt->SetFillColor(kYellow);  
   nutaupi_bkg_bdt->Draw("HIST");
   //   nutaueCC_bdt->Draw("same");
   TH1F* sum_pi_bdt = (TH1F*)nutaupi_bkg_bdt->Clone("sum_bdt");
   sum_pi_bdt->Add(nutaumuCC_bdt);
   sum_pi_bdt->SetLineColor(kRed);  // Set color for the sum histogram
   sum_pi_bdt->SetLineWidth(2);     // Set line width for better visibility
   sum_pi_bdt->SetTitle("Sum of NC and nutaupiCC BDT");
   sum_pi_bdt->Draw("e SAME");
   c4pi->SaveAs("c4pi.png");

   TCanvas *c4rho = new TCanvas("c4rho", "bdt", 800, 600);
   gPad->SetLogy();
   nutaurho_bkg_bdt->SetFillColor(kYellow);  
   nutaurho_bkg_bdt->Draw("HIST");
   //   nutaueCC_bdt->Draw("same");
   TH1F* sum_rho_bdt = (TH1F*)nutaurho_bkg_bdt->Clone("sum_bdt");
   sum_rho_bdt->Add(nutaumuCC_bdt);
   sum_rho_bdt->SetLineColor(kRed);  // Set color for the sum histogram
   sum_rho_bdt->SetLineWidth(2);     // Set line width for better visibility
   sum_rho_bdt->SetTitle("Sum of NC and nutaurhoCC BDT");
   sum_rho_bdt->Draw("e SAME");
   c4rho->SaveAs("c4rho.png");

   dump_cuts(&tauecuts);
   dump_cuts(&taumucuts);
   dump_cuts(&taupicuts);
   dump_cuts(&taurhocuts);

   rescale_cuts(&tauecuts);
   rescale_cuts(&taumucuts);
   rescale_cuts(&taupicuts);
   rescale_cuts(&taurhocuts);
   
   dump_cuts(&tauecuts);
   dump_cuts(&taumucuts);
   dump_cuts(&taupicuts);
   dump_cuts(&taurhocuts);

 }
