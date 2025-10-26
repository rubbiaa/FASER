#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

void plot_muon_genfit_taubin(const char* inpath = "build/scifi_hits_all.root", const char* outdir_override = "") {
  gStyle->SetOptStat(0);

  // Try primary path; fallback to local directory
  std::string path = inpath;
  if (gSystem->AccessPathName(path.c_str())) {
    if (!gSystem->AccessPathName("scifi_hits_all.root")) path = "scifi_hits_all.root";
  }

  // Determine output directory: if override provided, use it; otherwise infer from filename suffix *_<TAG>.root where TAG ends with GeV
  auto inferTagFrom = [](const std::string& p) -> std::string {
    // strip directories
    size_t s = p.find_last_of("/\\");
    std::string base = (s==std::string::npos) ? p : p.substr(s+1);
    // strip extension
    size_t dot = base.rfind('.');
    std::string stem = (dot==std::string::npos) ? base : base.substr(0, dot);
    // take last underscore suffix
    size_t us = stem.rfind('_');
    if (us == std::string::npos) return std::string();
    std::string tag = stem.substr(us+1);
    if (tag.find("GeV") != std::string::npos) return tag; // use tag like 50GeV, 100GeV
    return std::string();
  };

  std::string tag;
  std::string outDir = "build";
  if (outdir_override && std::string(outdir_override).size() > 0) {
    outDir = std::string(outdir_override);
  } else {
    tag = inferTagFrom(path);
    if (!tag.empty()) outDir = std::string("build/") + tag;
  }
  // ensure output directory exists
  gSystem->mkdir(outDir.c_str(), true);

  TFile* fin = TFile::Open(path.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << path << std::endl;
    return;
  }
  TTree* t = (TTree*)fin->Get("Hits");
  if (!t) { std::cerr << "Tree 'Hits' not found in file." << std::endl; return; }

  // Branch pointers (vectors per event)
  std::vector<double> *muon_truth_p = nullptr;
  std::vector<double> *muon_truth_px = nullptr, *muon_truth_py = nullptr, *muon_truth_pz = nullptr;
  std::vector<double> *genfit_p = nullptr, *genfit_px = nullptr, *genfit_py = nullptr, *genfit_pz = nullptr;
  std::vector<double> *genfit_chi2 = nullptr, *genfit_ndf = nullptr, *genfit_pval = nullptr;
  std::vector<bool>   *genfit_fit_success = nullptr;
  std::vector<double> *taubin_p = nullptr, *taubin_p_err = nullptr, *taubin_px = nullptr, *taubin_py = nullptr, *taubin_pz = nullptr;
  std::vector<double> *taubin_chi2 = nullptr, *taubin_ndf = nullptr, *taubin_pval = nullptr;
  std::vector<int>    *taubin_charge = nullptr;
  std::vector<bool>   *taubin_fit_success = nullptr;
  std::vector<int>    *muon_nhits = nullptr, *muon_nstations = nullptr;

  t->SetBranchAddress("muon_truth_p", &muon_truth_p);
  t->SetBranchAddress("muon_truth_px", &muon_truth_px);
  t->SetBranchAddress("muon_truth_py", &muon_truth_py);
  t->SetBranchAddress("muon_truth_pz", &muon_truth_pz);

  t->SetBranchAddress("genfit_p", &genfit_p);
  t->SetBranchAddress("genfit_px", &genfit_px);
  t->SetBranchAddress("genfit_py", &genfit_py);
  t->SetBranchAddress("genfit_pz", &genfit_pz);
  t->SetBranchAddress("genfit_chi2", &genfit_chi2);
  t->SetBranchAddress("genfit_ndf", &genfit_ndf);
  t->SetBranchAddress("genfit_pval", &genfit_pval);
  t->SetBranchAddress("genfit_fit_success", &genfit_fit_success);

  t->SetBranchAddress("taubin_p", &taubin_p);
  t->SetBranchAddress("taubin_p_err", &taubin_p_err);
  t->SetBranchAddress("taubin_px", &taubin_px);
  t->SetBranchAddress("taubin_py", &taubin_py);
  t->SetBranchAddress("taubin_pz", &taubin_pz);
  t->SetBranchAddress("taubin_chi2", &taubin_chi2);
  t->SetBranchAddress("taubin_ndf", &taubin_ndf);
  t->SetBranchAddress("taubin_pval", &taubin_pval);
  t->SetBranchAddress("taubin_charge", &taubin_charge);
  t->SetBranchAddress("taubin_fit_success", &taubin_fit_success);
  // extra axes for resolution vs track complexity
  t->SetBranchAddress("muon_nhits", &muon_nhits);
  t->SetBranchAddress("muon_nstations", &muon_nstations);

  // Histograms
  TH2D* h_p_truth_vs_genfit = new TH2D("h_p_truth_vs_genfit", "Truth p vs GenFit p;Truth p [GeV/c];GenFit p [GeV/c]", 100, 0, 200, 100, 0, 200);
  TH2D* h_p_truth_vs_taubin = new TH2D("h_p_truth_vs_taubin", "Truth p vs Taubin p;Truth p [GeV/c];Taubin p [GeV/c]", 100, 0, 200, 100, 0, 200);

  TH1D* h_p_res_genfit = new TH1D("h_p_res_genfit", "(p_{reco}-p_{truth})/p_{truth} (GenFit);Rel. error;Tracks", 100, -1, 1);
  TH1D* h_p_res_taubin = new TH1D("h_p_res_taubin", "(p_{reco}-p_{truth})/p_{truth} (Taubin);Rel. error;Tracks", 100, -1, 1);

  // Widen x-ranges to better capture long tails observed in residuals
  TH1D* h_px_res_genfit = new TH1D("h_px_res_genfit", "(p_{x,reco}-p_{x,truth});#Delta p_{x} [GeV/c];Tracks", 200, -50, 50);
  TH1D* h_px_res_taubin = new TH1D("h_px_res_taubin", "(p_{x,reco}-p_{x,truth});#Delta p_{x} [GeV/c];Tracks", 200, -50, 50);
  TH1D* h_py_res_genfit = new TH1D("h_py_res_genfit", "(p_{y,reco}-p_{y,truth});#Delta p_{y} [GeV/c];Tracks", 200, -50, 50);
  TH1D* h_py_res_taubin = new TH1D("h_py_res_taubin", "(p_{y,reco}-p_{y,truth});#Delta p_{y} [GeV/c];Tracks", 200, -50, 50);
  TH1D* h_pz_res_genfit = new TH1D("h_pz_res_genfit", "(p_{z,reco}-p_{z,truth});#Delta p_{z} [GeV/c];Tracks", 200, -100, 100);
  TH1D* h_pz_res_taubin = new TH1D("h_pz_res_taubin", "(p_{z,reco}-p_{z,truth});#Delta p_{z} [GeV/c];Tracks", 200, -100, 100);

  // 2D: relative p resolution vs number of stations / hits
  // x binning: stations typically up to ~11; hits equals stations after averaging in our pipeline
  TH2D* h2_prel_genfit_vs_st = new TH2D("h2_prel_genfit_vs_st", "(p_{reco}-p_{truth})/p_{truth} vs nstations (GenFit);nstations;Rel. error", 13, -0.5, 12.5, 200, -1.0, 1.0);
  TH2D* h2_prel_taubin_vs_st = new TH2D("h2_prel_taubin_vs_st", "(p_{reco}-p_{truth})/p_{truth} vs nstations (Taubin);nstations;Rel. error", 13, -0.5, 12.5, 200, -1.0, 1.0);
  TH2D* h2_prel_genfit_vs_nh = new TH2D("h2_prel_genfit_vs_nh", "(p_{reco}-p_{truth})/p_{truth} vs nhits (GenFit);nhits;Rel. error", 45, -0.5, 44.5, 200, -1.0, 1.0);
  TH2D* h2_prel_taubin_vs_nh = new TH2D("h2_prel_taubin_vs_nh", "(p_{reco}-p_{truth})/p_{truth} vs nhits (Taubin);nhits;Rel. error", 45, -0.5, 44.5, 200, -1.0, 1.0);

  Long64_t nent = t->GetEntries();
  Long64_t used_genfit = 0, used_taubin = 0;

  for (Long64_t ie = 0; ie < nent; ++ie) {
    t->GetEntry(ie);
    if (!muon_truth_p) continue;
    size_t ntrk = muon_truth_p->size();
    for (size_t i = 0; i < ntrk; ++i) {
      double p_truth  = muon_truth_p->at(i);
      double px_truth = muon_truth_px->at(i);
      double py_truth = muon_truth_py->at(i);
      double pz_truth = muon_truth_pz->at(i);
      int nst = (muon_nstations && i < muon_nstations->size()) ? muon_nstations->at(i) : -1;
      int nh  = (muon_nhits && i < muon_nhits->size()) ? muon_nhits->at(i) : -1;

      // GenFit
      if (genfit_p && genfit_fit_success && i < genfit_p->size() && i < genfit_fit_success->size() && genfit_fit_success->at(i)) {
        double p_reco = genfit_p->at(i);
        double px_reco = genfit_px->at(i);
        double py_reco = genfit_py->at(i);
        double pz_reco = genfit_pz->at(i);
        if (p_truth > 1e-6) {
          double r = (p_reco - p_truth)/p_truth;
          h_p_res_genfit->Fill(r);
          if (nst >= 0) h2_prel_genfit_vs_st->Fill(nst, r);
          if (nh  >= 0) h2_prel_genfit_vs_nh->Fill(nh,  r);
        }
        h_p_truth_vs_genfit->Fill(p_truth, p_reco);
        h_px_res_genfit->Fill(px_reco - px_truth);
        h_py_res_genfit->Fill(py_reco - py_truth);
        h_pz_res_genfit->Fill(pz_reco - pz_truth);
        ++used_genfit;
      }

      // Taubin
      if (taubin_p && taubin_fit_success && i < taubin_p->size() && i < taubin_fit_success->size() && taubin_fit_success->at(i)) {
        double p_reco = taubin_p->at(i);
        double px_reco = taubin_px->at(i);
        double py_reco = taubin_py->at(i);
        double pz_reco = taubin_pz->at(i);
        if (p_truth > 1e-6) {
          double r = (p_reco - p_truth)/p_truth;
          h_p_res_taubin->Fill(r);
          if (nst >= 0) h2_prel_taubin_vs_st->Fill(nst, r);
          if (nh  >= 0) h2_prel_taubin_vs_nh->Fill(nh,  r);
        }
        h_p_truth_vs_taubin->Fill(p_truth, p_reco);
        h_px_res_taubin->Fill(px_reco - px_truth);
        h_py_res_taubin->Fill(py_reco - py_truth);
        h_pz_res_taubin->Fill(pz_reco - pz_truth);
        ++used_taubin;
      }
    }
  }

  // Create output canvases
  TCanvas* c1 = new TCanvas("c1","p truth vs reco", 1200, 500);
  c1->Divide(2,1);
  c1->cd(1); h_p_truth_vs_genfit->Draw("COLZ");
  c1->cd(2); h_p_truth_vs_taubin->Draw("COLZ");
  c1->SaveAs((outDir+"/muon_p_truth_vs_reco_genfit_taubin.png").c_str());
  c1->SaveAs((outDir+"/muon_p_truth_vs_reco_genfit_taubin.pdf").c_str());

  TCanvas* c2 = new TCanvas("c2","p relative resolution", 1200, 500);
  c2->Divide(2,1);
  c2->cd(1); h_p_res_genfit->SetLineColor(kBlue+1); h_p_res_genfit->Draw();
  c2->cd(2); h_p_res_taubin->SetLineColor(kRed+1); h_p_res_taubin->Draw();
  c2->SaveAs((outDir+"/muon_p_resolution_genfit_taubin.png").c_str());
  c2->SaveAs((outDir+"/muon_p_resolution_genfit_taubin.pdf").c_str());

  TCanvas* c3 = new TCanvas("c3","px/py/pz absolute residuals", 1500, 600);
  c3->Divide(3,2);
  c3->cd(1); h_px_res_genfit->SetLineColor(kBlue+1); h_px_res_genfit->Draw();
  c3->cd(2); h_px_res_taubin->SetLineColor(kRed+1); h_px_res_taubin->Draw();
  c3->cd(3); h_py_res_genfit->SetLineColor(kBlue+1); h_py_res_genfit->Draw();
  c3->cd(4); h_py_res_taubin->SetLineColor(kRed+1); h_py_res_taubin->Draw();
  c3->cd(5); h_pz_res_genfit->SetLineColor(kBlue+1); h_pz_res_genfit->Draw();
  c3->cd(6); h_pz_res_taubin->SetLineColor(kRed+1); h_pz_res_taubin->Draw();
  c3->SaveAs((outDir+"/muon_components_resolution_genfit_taubin.png").c_str());
  c3->SaveAs((outDir+"/muon_components_resolution_genfit_taubin.pdf").c_str());

  // Profiles vs stations
  TProfile* pr_pg_vs_st = h2_prel_genfit_vs_st->ProfileX("pr_pg_vs_st");
  TProfile* pr_pt_vs_st = h2_prel_taubin_vs_st->ProfileX("pr_pt_vs_st");
  pr_pg_vs_st->SetLineColor(kBlue+1); pr_pg_vs_st->SetMarkerColor(kBlue+1); pr_pg_vs_st->SetMarkerStyle(20);
  pr_pt_vs_st->SetLineColor(kRed+1);  pr_pt_vs_st->SetMarkerColor(kRed+1);  pr_pt_vs_st->SetMarkerStyle(21);
  pr_pg_vs_st->SetTitle("Mean rel. p error vs nstations;nstations;Mean (p_{reco}-p_{truth})/p_{truth}");
  pr_pg_vs_st->SetErrorOption("s"); // show RMS as error bars
  pr_pt_vs_st->SetErrorOption("s");
  TCanvas* c4 = new TCanvas("c4","p rel vs stations", 800, 600);
  pr_pg_vs_st->Draw("E1");
  pr_pt_vs_st->Draw("E1 SAME");
  auto leg4 = new TLegend(0.60,0.75,0.88,0.90); leg4->AddEntry(pr_pg_vs_st,"GenFit","lep"); leg4->AddEntry(pr_pt_vs_st,"Taubin","lep"); leg4->Draw();
  c4->SaveAs((outDir+"/muon_prel_vs_nstations_genfit_taubin.png").c_str());
  c4->SaveAs((outDir+"/muon_prel_vs_nstations_genfit_taubin.pdf").c_str());

  // Profiles vs hits
  TProfile* pr_pg_vs_nh = h2_prel_genfit_vs_nh->ProfileX("pr_pg_vs_nh");
  TProfile* pr_pt_vs_nh = h2_prel_taubin_vs_nh->ProfileX("pr_pt_vs_nh");
  pr_pg_vs_nh->SetLineColor(kBlue+1); pr_pg_vs_nh->SetMarkerColor(kBlue+1); pr_pg_vs_nh->SetMarkerStyle(20);
  pr_pt_vs_nh->SetLineColor(kRed+1);  pr_pt_vs_nh->SetMarkerColor(kRed+1);  pr_pt_vs_nh->SetMarkerStyle(21);
  pr_pg_vs_nh->SetTitle("Mean rel. p error vs nhits;nhits;Mean (p_{reco}-p_{truth})/p_{truth}");
  pr_pg_vs_nh->SetErrorOption("s");
  pr_pt_vs_nh->SetErrorOption("s");
  TCanvas* c5 = new TCanvas("c5","p rel vs hits", 800, 600);
  pr_pg_vs_nh->Draw("E1");
  pr_pt_vs_nh->Draw("E1 SAME");
  auto leg5 = new TLegend(0.60,0.75,0.88,0.90); leg5->AddEntry(pr_pg_vs_nh,"GenFit","lep"); leg5->AddEntry(pr_pt_vs_nh,"Taubin","lep"); leg5->Draw();
  c5->SaveAs((outDir+"/muon_prel_vs_nhits_genfit_taubin.png").c_str());
  c5->SaveAs((outDir+"/muon_prel_vs_nhits_genfit_taubin.pdf").c_str());

  // CSV summaries for profiles (bin center, entries, mean, RMS as error)
  auto dump_profile_csv = [](TProfile* p, const char* path){
    std::ofstream f(path);
    f << "bin, x_center, entries, mean, rms\n";
    for (int i=1;i<=p->GetNbinsX();++i){
      f << i << "," << p->GetBinCenter(i) << "," << p->GetBinEntries(i) << "," << p->GetBinContent(i) << "," << p->GetBinError(i) << "\n";
    }
    f.close();
  };
  dump_profile_csv(pr_pg_vs_st, (outDir+"/muon_prel_vs_nstations_genfit.csv").c_str());
  dump_profile_csv(pr_pt_vs_st, (outDir+"/muon_prel_vs_nstations_taubin.csv").c_str());
  dump_profile_csv(pr_pg_vs_nh, (outDir+"/muon_prel_vs_nhits_genfit.csv").c_str());
  dump_profile_csv(pr_pt_vs_nh, (outDir+"/muon_prel_vs_nhits_taubin.csv").c_str());

  // Resolution table (CSV)
  auto mean_sigma = [](TH1D* h){
    struct S { double n, mean, sigma; } s; s.n=h->GetEntries(); s.mean=h->GetMean(); s.sigma=h->GetRMS(); return s; };
  auto s_pg  = mean_sigma(h_p_res_genfit);
  auto s_pt  = mean_sigma(h_p_res_taubin);
  auto s_pxg = mean_sigma(h_px_res_genfit);
  auto s_pxt = mean_sigma(h_px_res_taubin);
  auto s_pyg = mean_sigma(h_py_res_genfit);
  auto s_pyt = mean_sigma(h_py_res_taubin);
  auto s_pzg = mean_sigma(h_pz_res_genfit);
  auto s_pzt = mean_sigma(h_pz_res_taubin);

  std::ofstream csv((outDir+"/muon_resolution_table.csv").c_str());
  csv << "method,param,entries,mean_bias,sigma\n";
  csv << "genfit,p_rel," << s_pg.n  << "," << s_pg.mean  << "," << s_pg.sigma  << "\n";
  csv << "taubin,p_rel," << s_pt.n  << "," << s_pt.mean  << "," << s_pt.sigma  << "\n";
  csv << "genfit,px_abs,"<< s_pxg.n << "," << s_pxg.mean << "," << s_pxg.sigma << "\n";
  csv << "taubin,px_abs,"<< s_pxt.n << "," << s_pxt.mean << "," << s_pxt.sigma << "\n";
  csv << "genfit,py_abs,"<< s_pyg.n << "," << s_pyg.mean << "," << s_pyg.sigma << "\n";
  csv << "taubin,py_abs,"<< s_pyt.n << "," << s_pyt.mean << "," << s_pyt.sigma << "\n";
  csv << "genfit,pz_abs,"<< s_pzg.n << "," << s_pzg.mean << "," << s_pzg.sigma << "\n";
  csv << "taubin,pz_abs,"<< s_pzt.n << "," << s_pzt.mean << "," << s_pzt.sigma << "\n";
  csv.close();

  std::cout << "Saved plots to " << outDir << "/ and table to " << outDir << "/muon_resolution_table.csv\n";
  std::cout << "Saved profiles: prel vs nstations/nhits (PNG/PDF) and CSVs in " << outDir << "/\n";
  std::cout << "Tracks used: GenFit=" << used_genfit << ", Taubin=" << used_taubin << std::endl;

  fin->Close();
}
