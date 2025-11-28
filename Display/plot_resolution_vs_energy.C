#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TStyle.h>
#include <TPaveText.h>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

struct ResPoint {
  double energyGeV{0.0};
  // genfit
  double gf_entries{0.0};
  double gf_mean{0.0};
  double gf_sigma{0.0};
  // taubin
  double tb_entries{0.0};
  double tb_mean{0.0};
  double tb_sigma{0.0};
};

static bool endsWith(const std::string &s, const std::string &suffix) {
  return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string trim(const std::string &s) {
  size_t i = 0, j = s.size();
  while (i < j && std::isspace((unsigned char)s[i])) ++i;
  while (j > i && std::isspace((unsigned char)s[j-1])) --j;
  return s.substr(i, j-i);
}

static double parseEnergyTag(const std::string &dirName) {
  // Expect something like "100GeV" -> 100
  if (!endsWith(dirName, "GeV")) return -1.0;
  std::string num = dirName.substr(0, dirName.size()-3);
  // use strtod which does not throw
  const char* c = num.c_str();
  char* endptr = nullptr;
  double val = std::strtod(c, &endptr);
  if (endptr == c) return -1.0; // no conversion
  return val;
}

static bool readSummaryCSV(const std::string &csvPath, ResPoint &out) {
  std::ifstream in(csvPath);
  if (!in) return false;
  std::string line;
  // header
  if (!std::getline(in, line)) return false;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty()) continue;
    // format: method,param,entries,mean_bias,sigma
    std::istringstream ss(line);
    std::string method, param, entries, mean, sigma;
    if (!std::getline(ss, method, ',')) continue;
    if (!std::getline(ss, param, ',')) continue;
    if (!std::getline(ss, entries, ',')) continue;
    if (!std::getline(ss, mean, ',')) continue;
    if (!std::getline(ss, sigma, ',')) continue;
    method = trim(method); param = trim(param);
    if (param != "p_rel") continue; // we summarize relative p resolution
    double n = atof(entries.c_str());
    double m = atof(mean.c_str());
    double s = atof(sigma.c_str());
    if (method == "genfit") {
      out.gf_entries = n; out.gf_mean = m; out.gf_sigma = s;
    } else if (method == "taubin") {
      out.tb_entries = n; out.tb_mean = m; out.tb_sigma = s;
    }
  }
  return true;
}

void plot_resolution_vs_energy(const char* baseDir = "build", const char* outDir = "build/summary") {
  gStyle->SetOptStat(0);
  // ensure output directory exists
  gSystem->mkdir(outDir, true);

  // List subdirectories under baseDir
  std::vector<ResPoint> points;
  TSystemDirectory based("base", baseDir);
  TList* list = based.GetListOfFiles();
  if (!list) {
    std::cerr << "No files in baseDir: " << baseDir << std::endl;
  } else {
    TIter next(list);
    while (TObject* obj = next()) {
      TSystemFile* f = dynamic_cast<TSystemFile*>(obj);
      if (!f) continue;
      const char* n = f->GetName();
      if (!n || n[0]=='.') continue;
      if (!f->IsDirectory()) continue;
      std::string dname = n;
      if (!endsWith(dname, "GeV")) continue;
      double E = parseEnergyTag(dname);
      if (E <= 0) continue;
      std::string csvPath = std::string(baseDir) + "/" + dname + "/muon_resolution_table.csv";
      ResPoint rp; rp.energyGeV = E;
      if (readSummaryCSV(csvPath, rp)) {
        points.push_back(rp);
      }
    }
  }

  if (points.empty()) {
    std::cerr << "No energy-tagged resolution tables found under " << baseDir << std::endl;
    return;
  }

  // sort by energy
  std::sort(points.begin(), points.end(), [](const ResPoint&a, const ResPoint&b){return a.energyGeV < b.energyGeV;});

  // Prepare graphs (sigma of relative p error)
  std::vector<double> x, y_gf, y_tb, yb_gf, yb_tb;
  std::vector<double> ye_gf, ye_tb, ybe_gf, ybe_tb; // error bars
  x.reserve(points.size()); y_gf.reserve(points.size()); y_tb.reserve(points.size());
  yb_gf.reserve(points.size()); yb_tb.reserve(points.size());
  ye_gf.reserve(points.size()); ye_tb.reserve(points.size());
  ybe_gf.reserve(points.size()); ybe_tb.reserve(points.size());

  for (const auto &rp : points) {
    x.push_back(rp.energyGeV);
    y_gf.push_back(rp.gf_sigma);
    y_tb.push_back(rp.tb_sigma);
    yb_gf.push_back(rp.gf_mean);
    yb_tb.push_back(rp.tb_mean);
    // Errors:
    // - sigma uncertainty ~ sigma / sqrt(2*(N-1)) for normal distribution
    double gf_se_sigma = (rp.gf_entries > 1 ? rp.gf_sigma / std::sqrt(2.0 * (rp.gf_entries - 1.0)) : 0.0);
    double tb_se_sigma = (rp.tb_entries > 1 ? rp.tb_sigma / std::sqrt(2.0 * (rp.tb_entries - 1.0)) : 0.0);
    ye_gf.push_back(gf_se_sigma);
    ye_tb.push_back(tb_se_sigma);
    // - mean uncertainty ~ sigma / sqrt(N)
    double gf_se_mean = (rp.gf_entries > 0 ? rp.gf_sigma / std::sqrt(rp.gf_entries) : 0.0);
    double tb_se_mean = (rp.tb_entries > 0 ? rp.tb_sigma / std::sqrt(rp.tb_entries) : 0.0);
    ybe_gf.push_back(gf_se_mean);
    ybe_tb.push_back(tb_se_mean);
  }

  // Use TGraphErrors with symmetric y-errors; no x-errors
  TGraphErrors* gr_gf = new TGraphErrors((int)x.size());
  TGraphErrors* gr_tb = new TGraphErrors((int)x.size());
  gr_gf->SetName("gr_genfit_sigma");
  gr_tb->SetName("gr_taubin_sigma");
  for (int i=0;i<(int)x.size();++i){
    gr_gf->SetPoint(i, x[i], y_gf[i]); gr_gf->SetPointError(i, 0.0, ye_gf[i]);
    gr_tb->SetPoint(i, x[i], y_tb[i]); gr_tb->SetPointError(i, 0.0, ye_tb[i]);
  }
  gr_gf->SetLineColor(kBlue+1); gr_gf->SetMarkerColor(kBlue+1); gr_gf->SetMarkerStyle(20);
  gr_tb->SetLineColor(kRed+1);  gr_tb->SetMarkerColor(kRed+1);  gr_tb->SetMarkerStyle(21);

  // Plot sigma vs energy
  TCanvas* c1 = new TCanvas("c_res_vs_E", "Relative p resolution vs energy", 900, 600);
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(gr_gf, "LP");
  mg->Add(gr_tb, "LP");
  mg->SetTitle("Muon relative p resolution vs energy;Energy [GeV];#sigma((p_{reco}-p_{truth})/p_{truth})");
  mg->Draw("A");
  auto leg = new TLegend(0.60,0.75,0.88,0.90);
  leg->AddEntry(gr_gf, "GenFit", "lp");
  leg->AddEntry(gr_tb, "Taubin", "lp");
  leg->Draw();
  c1->SaveAs((std::string(outDir)+"/muon_prel_sigma_vs_energy.png").c_str());
  c1->SaveAs((std::string(outDir)+"/muon_prel_sigma_vs_energy.pdf").c_str());

  // Optional: bias vs energy (mean of relative error)
  TGraphErrors* gr_gf_bias = new TGraphErrors((int)x.size());
  TGraphErrors* gr_tb_bias = new TGraphErrors((int)x.size());
  gr_gf_bias->SetName("gr_genfit_bias");
  gr_tb_bias->SetName("gr_taubin_bias");
  for (int i=0;i<(int)x.size();++i){
    gr_gf_bias->SetPoint(i, x[i], yb_gf[i]); gr_gf_bias->SetPointError(i, 0.0, ybe_gf[i]);
    gr_tb_bias->SetPoint(i, x[i], yb_tb[i]); gr_tb_bias->SetPointError(i, 0.0, ybe_tb[i]);
  }
  gr_gf_bias->SetLineColor(kBlue+1); gr_gf_bias->SetMarkerColor(kBlue+1); gr_gf_bias->SetMarkerStyle(20);
  gr_tb_bias->SetLineColor(kRed+1);  gr_tb_bias->SetMarkerColor(kRed+1);  gr_tb_bias->SetMarkerStyle(21);

  TCanvas* c2 = new TCanvas("c_bias_vs_E", "Relative p bias vs energy", 900, 600);
  TMultiGraph* mgb = new TMultiGraph();
  mgb->Add(gr_gf_bias, "LP");
  mgb->Add(gr_tb_bias, "LP");
  mgb->SetTitle("Muon relative p bias vs energy;Energy [GeV];Mean((p_{reco}-p_{truth})/p_{truth})");
  mgb->Draw("A");
  auto legb = new TLegend(0.60,0.75,0.88,0.90);
  legb->AddEntry(gr_gf_bias, "GenFit", "lp");
  legb->AddEntry(gr_tb_bias, "Taubin", "lp");
  legb->Draw();
  c2->SaveAs((std::string(outDir)+"/muon_prel_bias_vs_energy.png").c_str());
  c2->SaveAs((std::string(outDir)+"/muon_prel_bias_vs_energy.pdf").c_str());

  // CSV summary
  std::ofstream csv((std::string(outDir)+"/muon_prel_vs_energy_summary.csv").c_str());
  csv << "energy_GeV,genfit_entries,genfit_mean,genfit_sigma,genfit_mean_err,genfit_sigma_err,taubin_entries,taubin_mean,taubin_sigma,taubin_mean_err,taubin_sigma_err\n";
  for (const auto &rp : points) {
    double gf_se_sigma = (rp.gf_entries > 1 ? rp.gf_sigma / std::sqrt(2.0 * (rp.gf_entries - 1.0)) : 0.0);
    double tb_se_sigma = (rp.tb_entries > 1 ? rp.tb_sigma / std::sqrt(2.0 * (rp.tb_entries - 1.0)) : 0.0);
    double gf_se_mean  = (rp.gf_entries > 0 ? rp.gf_sigma / std::sqrt(rp.gf_entries) : 0.0);
    double tb_se_mean  = (rp.tb_entries > 0 ? rp.tb_sigma / std::sqrt(rp.tb_entries) : 0.0);
    csv << rp.energyGeV << ","
        << rp.gf_entries << "," << rp.gf_mean << "," << rp.gf_sigma << "," << gf_se_mean  << "," << gf_se_sigma << ","
        << rp.tb_entries << "," << rp.tb_mean << "," << rp.tb_sigma << "," << tb_se_mean  << "," << tb_se_sigma << "\n";
  }
  csv.close();

  std::cout << "Saved energy summary plots to " << outDir << "/ and CSV to " << outDir << "/muon_prel_vs_energy_summary.csv\n";
}
