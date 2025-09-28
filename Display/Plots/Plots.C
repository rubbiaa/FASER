
#include "Riostream.h"

#include <algorithm>
#include <limits.h>
#include <list>
#include <iterator>


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "TROOT.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TText.h"
#include "TString.h"
#include "TSystem.h"


#ifdef __MAKECINT__
#pragma link off all class;
#pragma link C++ class FaserCalRecoPlot+;
#pragma link C++ class TreeManager+;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;
#endif

using namespace std;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Loading TreeManager in order to access TTree
// You just need to give the correct input file
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//root v6
R__LOAD_LIBRARY(TreeManager_C.so);
#include "TreeManager.h"

TreeManager* filereader = nullptr;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ReadFile()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();
  cout << "Reading all data " << entries << endl;
  for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
      data->GetEvent(ientry);
    } // end_of_for_ientry_loop
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void overlay_plot(TH1* h1, TH1* h2, const TString& title, const TString& xlabel, TCanvas* c, const TString& pdf_output)
{
  // Normalize both histograms if they have non-zero integral
  if (h1->Integral() > 0) h1->Scale(1.0 / h1->Integral());
  if (h2->Integral() > 0) h2->Scale(1.0 / h2->Integral());
  
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetStats(0);
  
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);
  h2->SetStats(0);

  h1->SetTitle(title);
  h1->GetXaxis()->SetTitle(xlabel);
  h1->GetYaxis()->SetTitle("Normalized entries");

  // Determine max y value and set axis range
  double max1 = h1->GetMaximum();
  double max2 = h2->GetMaximum();
  double maxY = std::max(max1, max2);
  h1->SetMaximum(1.2 * maxY);  // Add 20% padding

  h1->Draw("HIST");
  h2->Draw("HIST SAME");

  TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
  legend->AddEntry(h1, "Charm", "l");
  legend->AddEntry(h2, "Non-Charm", "l");
  legend->Draw();

  c->Print(pdf_output);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotCharmComparison() 
{
  FaserCalRecoPlot *data = filereader->tmCD();
  Long64_t entries = data->GetInputTree()->GetEntries();
  
  TH1F* hCharm = new TH1F("hCharm", "TotalMultiplicity: Charm vs Non-Charm;Multiplicity;Events", 30, 0, 30);
  TH1F* hNonCharm = new TH1F("hNonCharm", "TotalMultiplicity: Charm vs Non-Charm;Multiplicity;Events", 30, 0, 30);
  
  TH1F* h_VisibleEnergy_charm = new TH1F("h_VisibleEnergy_charm", "Visible Energy;E_{vis} [GeV];Events", 100, 0, 400);
  TH1F* h_VisibleEnergy_nocharm = new TH1F("h_VisibleEnergy_nocharm", "", 100, 0, 400);
  TH1F* h_VisibleEnergy_charm_multimu = new TH1F("h_VisibleEnergy_charm_multimu", "", 100, 0, 400);
  TH1F* h_VisibleEnergy_nocharm_multimu = new TH1F("h_VisibleEnergy_nocharm_multimu", "", 100, 0, 400);

  TH1F* h_TotalEnergy_charm = new TH1F("h_TotalEnergy_charm", "Total Energy;E_{vis} [GeV];Events", 100, 0, 10);
  TH1F* h_TotalEnergy_nocharm = new TH1F("h_TotalEnergy_nocharm", "", 100, 0, 10);
  TH1F* h_TotalEnergy_charm_multimu = new TH1F("h_TotalEnergy_charm_multimu", "", 100, 0, 10);
  TH1F* h_TotalEnergy_nocharm_multimu = new TH1F("h_TotalEnergy_nocharm_multimu", "", 100, 0, 10);

  // Create histograms for charm and non-charm MuTag variables
  TH1F* h_mom_charm = new TH1F("h_mom_charm", "MuTag Momentum;Momentum [GeV];Entries", 100, 0, 3000);
  TH1F* h_mom_nocharm = new TH1F("h_mom_nocharm", "", 100, 0, 3000);
  TH1F* h_mom_charm_multimu = new TH1F("h_mom_charm_multimu", "", 100, 0, 3000);
  TH1F* h_mom_nocharm_multimu = new TH1F("h_mom_nocharm_multimu", "", 100, 0, 3000);

  TH1F* h_ene_charm = new TH1F("h_ene_charm", "MuTag Energy;Energy [GeV];Entries", 100, 0, 3000);
  TH1F* h_ene_nocharm = new TH1F("h_ene_nocharm", "", 100, 0, 3000);
  TH1F* h_ene_charm_multimu = new TH1F("h_ene_charm_multimu", "", 100, 0, 3000);
  TH1F* h_ene_nocharm_multimu = new TH1F("h_ene_nocharm_multimu", "", 100, 0, 3000);
  
  TH1F* h_dist_charm = new TH1F("h_dist_charm", "MuTag Distance;Distance [mm];Entries", 100, 0, 1000);
  TH1F* h_dist_nocharm = new TH1F("h_dist_nocharm", "", 100, 0, 1000);
  TH1F* h_dist_charm_multimu = new TH1F("h_dist_charm_multimu", "", 100, 0, 1000);
  TH1F* h_dist_nocharm_multimu = new TH1F("h_dist_nocharm_multimu", "", 100, 0, 1000);

  TH1F* h_muon_charm = new TH1F("h_muon_charm", "MuTag Muon Count;# Muons;Entries", 10, 0, 10);
  TH1F* h_muon_nocharm = new TH1F("h_muon_nocharm", "", 10, 0, 10);

  // Rear Calorimeter energies
  TH1F* h_rearEcal_charm = new TH1F("h_rearEcal_charm", "Rear ECal Energy;Energy [GeV];Entries", 100, 0, 300);
  TH1F* h_rearEcal_nocharm = new TH1F("h_rearEcal_nocharm", "", 100, 0, 300);
  TH1F* h_rearEcal_charm_multimu = new TH1F("h_rearEcal_charm_multimu", "", 100, 0, 300);
  TH1F* h_rearEcal_nocharm_multimu = new TH1F("h_rearEcal_nocharm_multimu", "", 100, 0, 300);
  
  TH1F* h_rearHcal_charm = new TH1F("h_rearHcal_charm", "Rear HCal Energy;Energy [GeV];Entries", 100, 0, 50);
  TH1F* h_rearHcal_nocharm = new TH1F("h_rearHcal_nocharm", "", 100, 0, 50);
  TH1F* h_rearHcal_charm_multimu = new TH1F("h_rearHcal_charm_multimu", "", 100, 0, 50);
  TH1F* h_rearHcal_nocharm_multimu = new TH1F("h_rearHcal_nocharm_multimu", "", 100, 0, 50);
  
  TH1F* h_rearMucal_charm = new TH1F("h_rearMucal_charm", "Rear MuCal Energy;Energy [GeV];Entries", 100, 0, 1000);
  TH1F* h_rearMucal_nocharm = new TH1F("h_rearMucal_nocharm", "", 100, 0, 1000);
  TH1F* h_rearMucal_charm_multimu = new TH1F("h_rearMucal_charm_multimu", "", 100, 0, 1000);
  TH1F* h_rearMucal_nocharm_multimu = new TH1F("h_rearMucal_nocharm_multimu", "", 100, 0, 1000);
  
  for (Long64_t i = 0; i < entries; ++i) {
    data->GetEvent(i);
    bool is_charm = (data->CharmParticleType() != 0); // charm event
    bool is_multimu = (data->MuonTagStatus() >= 2); // multimuon

    if (is_charm)  // charm event
      {
	std::cout << "CharmEvent" << std::endl;
	hCharm->Fill(data->Multiplicity());
	h_VisibleEnergy_charm->Fill(data->VisibleE());
	h_TotalEnergy_charm->Fill(data->TotalE());
   
	h_rearEcal_charm->Fill(data->RearECalE());
	h_rearHcal_charm->Fill(data->RearHCalE());
	h_rearMucal_charm->Fill(data->RearMuCalE());
	if (is_multimu)
	  {
	    h_rearEcal_charm_multimu->Fill(data->RearECalE());
	    h_rearHcal_charm_multimu->Fill(data->RearHCalE());
	    h_rearMucal_charm_multimu->Fill(data->RearMuCalE());
	    h_VisibleEnergy_charm_multimu->Fill(data->VisibleE());
	    h_TotalEnergy_charm_multimu->Fill(data->TotalE());
	  }
      }
    else
      {
	// non-charm
	hNonCharm->Fill(data->Multiplicity());
	h_VisibleEnergy_nocharm->Fill(data->VisibleE());
	h_TotalEnergy_nocharm->Fill(data->TotalE());
   
	h_rearEcal_nocharm->Fill(data->RearECalE());
	h_rearHcal_nocharm->Fill(data->RearHCalE());
	h_rearMucal_nocharm->Fill(data->RearMuCalE());      
	if (is_multimu)
	  {
	    h_rearEcal_nocharm_multimu->Fill(data->RearECalE());
	    h_rearHcal_nocharm_multimu->Fill(data->RearHCalE());
	    h_rearMucal_nocharm_multimu->Fill(data->RearMuCalE());
	    h_VisibleEnergy_nocharm_multimu->Fill(data->VisibleE());
	    h_TotalEnergy_nocharm_multimu->Fill(data->TotalE());

	  }
       }
    /////////////////////////////////////////////
    const auto& mom = data->MuTagMom();
    const auto& ene = data->MuTagEne();
    const auto& dist = data->MuTagDist();
    
    for (size_t i = 0; i < mom.size(); ++i) {
      if (is_charm) {
	h_mom_charm->Fill(mom[i]/1000); // convert to MeV
	h_ene_charm->Fill(ene[i]/1000);
	h_dist_charm->Fill(dist[i]);
	if (is_multimu)
	  {
	    h_mom_charm_multimu->Fill(mom[i]/1000);
	    h_ene_charm_multimu->Fill(ene[i]/1000);
	    h_dist_charm_multimu->Fill(dist[i]);
	  }
      } else {
	h_mom_nocharm->Fill(mom[i]/1000);
	h_ene_nocharm->Fill(ene[i]/1000);
	h_dist_nocharm->Fill(dist[i]);
	if (is_multimu)
	  {
	    h_mom_nocharm_multimu->Fill(mom[i]/1000);
	    h_ene_nocharm_multimu->Fill(ene[i]/1000);
	    h_dist_nocharm_multimu->Fill(dist[i]);
	  }
      }
      // MuTag_muon
      if (is_charm) h_muon_charm->Fill(data->MuonTagStatus());
      else h_muon_nocharm->Fill(data->MuonTagStatus());
      /////////////////////////////////////////////
    }    
  } 
  
  TCanvas* c = new TCanvas("c", "Overlay", 800, 600);
  TString pdf = "output.pdf";
  
  c->Print(pdf + "["); // open PDF file to store all the plots




  overlay_plot(h_muon_charm, h_muon_nocharm, "MuTag Muon Count", "# Muons", c, pdf);

  overlay_plot(h_mom_charm, h_mom_nocharm, "MuTag Momentum", "p [GeV]", c, pdf);
  overlay_plot(h_mom_charm_multimu, h_mom_nocharm_multimu, "MuTag Momentum [multi-Muon]", "p [GeV]", c, pdf);
 
  overlay_plot(h_ene_charm, h_ene_nocharm, "MuTag Energy", "E [GeV]", c, pdf);
  overlay_plot(h_ene_charm_multimu, h_ene_nocharm_multimu, "MuTag Energy [multi-Muon]", "E [GeV]", c, pdf);
  
  overlay_plot(h_dist_charm, h_dist_nocharm, "MuTag Distance", "Distance [mm]", c, pdf);
  overlay_plot(h_dist_charm_multimu, h_dist_nocharm_multimu, "MuTag Distance [multi-Muon]", "Distance [cm]", c, pdf);
  
  overlay_plot(h_rearEcal_charm, h_rearEcal_nocharm, "Rear ECal Energy", "Energy [GeV]", c, pdf);
  overlay_plot(h_rearEcal_charm_multimu, h_rearEcal_nocharm_multimu, "Rear ECal Energy [multi-Muon]", "Energy [GeV]", c, pdf);

  overlay_plot(h_rearHcal_charm, h_rearHcal_nocharm, "Rear HCal Energy", "Energy [GeV]", c, pdf);
  overlay_plot(h_rearHcal_charm_multimu, h_rearHcal_nocharm_multimu, "Rear HCal Energy [multi-Muon]", "Energy [GeV]", c, pdf);

  overlay_plot(h_rearMucal_charm, h_rearMucal_nocharm, "Rear MuCal Energy", "Energy [GeV]", c, pdf);
  overlay_plot(h_rearMucal_charm_multimu, h_rearMucal_nocharm_multimu, "Rear MuCal Energy [multi-Muon]", "Energy [GeV]", c, pdf);
 
  overlay_plot (h_VisibleEnergy_charm, h_VisibleEnergy_nocharm, "Visible Energy", "Energy [GeV]", c, pdf);
  overlay_plot (h_VisibleEnergy_charm_multimu, h_VisibleEnergy_nocharm_multimu, "Visible Energy [multi-Muon]", "Energy [GeV]", c, pdf);

  overlay_plot (h_TotalEnergy_charm, h_TotalEnergy_nocharm, "Total Energy", "Energy [GeV]", c, pdf);
  overlay_plot (h_TotalEnergy_charm_multimu, h_TotalEnergy_nocharm_multimu, "Total Energy [multi-Muon]", "Energy [GeV]", c, pdf);
  

  c->Print(pdf + "]"); // close PDF

}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetSummaryTable()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();
  // Maps for counting
  std::map<int, int> cc_nc_count;
  std::map<int, int> neutrino_type_count;
  std::map<std::string, int> charm_particle_count;
  std::map<std::string, int> decay_mode_count;
  
  for (Long64_t ientry = 0; ientry < entries; ++ientry)
    {
      data->GetEvent(ientry);
      
      bool is_charm = (data->CharmParticleType() != 0); // charm event

      // CC / NC (0 = NC, 1 = CC)
      cc_nc_count[data->CCNCType()]++;
      
      // Neutrino type (PDG)
      neutrino_type_count[data->NuType()]++;
      
      if (is_charm)
	{
	  charm_particle_count[data->CharmParticleName()]++;
	  decay_mode_count[data->CharmDecayM()]++;
	}
    }
  
  // Print tables
  cout << "\n=== CC / NC Summary ===" << endl;
  for (auto& [type, count] : cc_nc_count)
    cout << (type == 0 ? "NC" : "CC") << ": " << count << endl;
  
  cout << "\n=== Neutrino Type Summary (PDG) ===" << endl;
  for (auto& [pdg, count] : neutrino_type_count)
    cout << "PDG " << pdg << ": " << count << endl;
  
  cout << "\n=== Charm Particle Types ===" << endl;
  for (auto& [name, count] : charm_particle_count)
    cout << name << ": " << count << endl;
  
  cout << "\n=== Charm Decay Modes ===" << endl;
  for (auto& [mode, count] : decay_mode_count)
    cout << mode << ": " << count << endl;
  
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetDetailedSummaryTable()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();

  // Nested maps: [neutrino PDG][CC/NC or CharmType] = count
  std::map<int, std::map<int, int>> cc_nc_by_nu;         // [nu_pdg][0 or 1]
  std::map<int, std::map<std::string, int>> charm_by_nu; // [nu_pdg][charm name]
  std::map<int, std::map<std::string, int>> decay_by_nu; // [nu_pdg][decay mode]

  for (Long64_t ientry = 0; ientry < entries; ++ientry)
  {
    data->GetEvent(ientry);

    int nu_pdg = data->NuType();
    int cc_type = data->CCNCType(); // 0 = NC, 1 = CC
    bool is_charm = (data->CharmParticleType() != 0);

    // Count CC/NC by neutrino type
    cc_nc_by_nu[nu_pdg][cc_type]++;

    // Count charm-related breakdowns
    if (is_charm)
    {
      charm_by_nu[nu_pdg][data->CharmParticleName()]++;
      decay_by_nu[nu_pdg][data->CharmDecayM()]++;
    }
  }

  // Helper lambda to map PDG to readable label
  auto pdg_label = [](int pdg) -> std::string {
    switch (pdg) {
      case 12: return "nu_e";
      case -12: return "nu_e_bar";
      case 14: return "nu_mu";
      case -14: return "nu_mu_bar";
      case 16: return "nu_tau";
      case -16: return "nu_tau_bar";
      default: return "other";
    }
  };

  // CC/NC Breakdown
  std::cout << "\n=== CC / NC Breakdown by Neutrino Type ===\n";
  for (const auto& [pdg, type_map] : cc_nc_by_nu)
  {
    std::cout << "Neutrino: " << pdg_label(pdg) << " (PDG " << pdg << ")\n";
    for (const auto& [ccnc, count] : type_map)
      std::cout << "  " << (ccnc == 1 ? "CC" : "NC") << ": " << count << "\n";
  }

  // Charm Particle Breakdown
  std::cout << "\n=== Charm Particle Types by Neutrino Type ===\n";
  for (const auto& [pdg, name_map] : charm_by_nu)
  {
    std::cout << "Neutrino: " << pdg_label(pdg) << " (PDG " << pdg << ")\n";
    for (const auto& [name, count] : name_map)
      std::cout << "  " << name << ": " << count << "\n";
  }

  // Decay Mode Breakdown
  std::cout << "\n=== Charm Decay Modes by Neutrino Type ===\n";
  for (const auto& [pdg, decay_map] : decay_by_nu)
  {
    std::cout << "Neutrino: " << pdg_label(pdg) << " (PDG " << pdg << ")\n";
    for (const auto& [mode, count] : decay_map)
      std::cout << "  " << mode << ": " << count << "\n";
  }
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetSummaryTablePDF()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();

  // Count maps: category ? type ? count
  std::map<std::string, std::map<std::string, int>> summary;

  auto classify = [](int pdg) {
    if (pdg == 12) return "nu_e";
    if (pdg == -12) return "nu_e_bar";
    if (pdg == 14) return "nu_mu";
    if (pdg == -14) return "nu_mu_bar";
    if (pdg == 16) return "nu_tau";
    if (pdg == -16) return "nu_tau_bar";
    return "other";
  };

  for (Long64_t ientry = 0; ientry < entries; ++ientry)
  {
    data->GetEvent(ientry);
    std::string nutype = classify(data->NuType());
    if (nutype == "other") continue;

    bool is_charm = data->CharmParticleType() != 0;
    std::string cc_label = (data->CCNCType() == 1 ? "CC" : "NC");

    summary[cc_label][nutype]++;
    if (is_charm)
      summary["Charm"][nutype]++;
  }

  // ===> Prepare for display
  auto get_value = [&](const std::map<std::string, int>& m, const std::string& k) {
    return m.count(k) ? m.at(k) : 0;
  };

  // ===> Create Canvas
  TCanvas* c = new TCanvas("c", "Summary Table", 800, 600);
  TPaveText* table = new TPaveText(0.01, 0.1, 0.99, 0.9);
  table->SetTextFont(42);
  table->SetTextAlign(22);
  table->SetTextSize(0.03);
  table->SetFillColor(0);

  // Header
  table->AddText("Summary Table: nu_l | nu_l_bar | nu_l + nu_l_bar");
  table->AddText("--------------------------------------------------");
  table->AddText("Category         nu      nu_bar      Total");

  // Utility to format each row
  auto add_row = [&](const std::string& label, const std::string& x, const std::string& x_bar) {
    int a = get_value(summary[label], x);
    int b = get_value(summary[label], x_bar);
    int sum = a + b;
    TString line = Form("%-15s %6d  %6d  %6d", (label + " " + x.substr(3)).c_str(), a, b, sum);
    table->AddText(line);
  };

  // Add rows
  for (const std::string& label : {"CC", "NC", "Charm"})
  {
    add_row(label, "nu_e", "nu_e_bar");
    add_row(label, "nu_mu", "nu_mu_bar");
    add_row(label, "nu_tau", "nu_tau_bar");
  }

  table->AddText("--------------------------------------------------");

  c->cd();
  table->Draw();
  c->Print("Summary4ColumnTable.pdf");
  std::cout << "\n Saved 4-column summary table to SummaryTable.pdf\n";
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetAllSummaryTablePDF()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();
  
  // === Classification functions
  auto classify = [](int pdg) {
    if (pdg == 12) return "nu_e";
    if (pdg == -12) return "nu_e_bar";
    if (pdg == 14) return "nu_mu";
    if (pdg == -14) return "nu_mu_bar";
    if (pdg == 16) return "nu_tau";
    if (pdg == -16) return "nu_tau_bar";
    return "other";
  };
  
  // === Summary Containers
  std::map<std::string, std::map<std::string, int>> summary;                // CC/NC/Charm
  std::map<std::string, std::map<std::string, int>> charm_particle_by_flavor; // Particle
  std::map<std::string, std::map<std::string, int>> decay_mode_by_flavor;     // Decay
  
  // === Loop over events
  for (Long64_t ientry = 0; ientry < entries; ++ientry)
    {
      data->GetEvent(ientry);
      std::string nutype = classify(data->NuType());
      if (nutype == "other") continue;
      
      bool is_charm = data->CharmParticleType() != 0;
      std::string cc_label = (data->CCNCType() == 1 ? "CC" : "NC");
      
      summary[cc_label][nutype]++;
      if (is_charm) {
	summary["Charm"][nutype]++;
	charm_particle_by_flavor[data->CharmParticleName()][nutype]++;
	decay_mode_by_flavor[data->CharmDecayM()][nutype]++;
      }
    }
  //////
  auto get_val = [](const std::map<std::string, int>& m, const std::string& key) {
    return m.count(key) ? m.at(key) : 0;
  };
  //////
  auto add_table_page = [&](const std::string& title,
			    const std::map<std::string, std::map<std::string, int>>& table_data,
			    const std::vector<std::string>& flav_order,
			    const std::string& filename,
			    bool first = false,
			    bool last = false) {
    auto get_val_local = get_val;


    TCanvas* c = new TCanvas(TString::Format("c_%s", title.c_str()), title.c_str(), 1000, 800);
    TPaveText* text = new TPaveText(0.01, 0.05, 0.99, 0.95);
    text->SetTextFont(42);
    text->SetTextAlign(12);
    text->SetTextSize(0.02);
    text->SetFillColor(0);
    
    TString head = Form("%-25s %10s %10s %10s %10s %10s %10s %10s",
			title.c_str(),
			"nu_e", "nubar_e", "nu_mu", "nubar_mu",
			"nu_tau", "nubar_tau", "Total");
    text->AddText(head);
    text->AddText("------------------------------------------------------------------------------------------");

    // DEBUG: Print to console
    std::cout << "\n==== Table: " << title << " ====\n";
    std::cout << head << std::endl;
    std::cout << "------------------------------------------------------------------------------------------\n";
    
    for (const auto& [label, counts] : table_data) {
      int total = 0;
      TString row = Form("%-25s", label.c_str());
      // DEBUG
      std::cout << label;

      for (const auto& flav : flav_order) {
	int n = get_val(counts, flav);
	total += n;
	row += Form(" %10d", n);
      }
      row += Form(" %10d", total);
      text->AddText(row);
      // DEBUG
    std::cout << " -> " << row << std::endl; // Console output
    }
    
    c->cd();
    text->Draw();
    gPad->Update();
    gSystem->ProcessEvents();
    c->Update();
    //TString opt = first ? "(" : (last ? ")" : "");
    //c->Print(filename.c_str(), opt);
    std::string opt = "pdf";
    if (first)  opt = "pdf(";
    if (last)   opt = "pdf)";
    c->Print(filename.c_str(), opt.c_str());
    delete c;
  };
  
  // === Neutrino flavours
  std::vector<std::string> flavors = {
    "nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar"
  };
  
  // === PDF Filename
  TString outname = "SummaryAllTables.pdf";
  
  // === Page 1
  {
    std::cout << " I am here " << std::endl;
    TCanvas* c1 = new TCanvas("c1", "Summary", 800, 600);
    TPaveText* table = new TPaveText(0.01, 0.1, 0.99, 0.9);
    table->SetTextFont(42);
    table->SetTextAlign(22);
    table->SetTextSize(0.03);
    table->SetFillColor(0);
    
    table->AddText("Summary Table: nu_x | nubar_x | nu_x + nubar_x");
    table->AddText("--------------------------------------------------");
    table->AddText("Category         nu       nubar      Total");
    
    auto get = [&](const std::string& cat, const std::string& flav1, const std::string& flav2) {
      int a = get_val(summary[cat], flav1);
      int b = get_val(summary[cat], flav2);
      return std::tuple<int,int,int>(a, b, a + b);
    };
    
    for (const std::string& cat : {"CC", "NC", "Charm"}) {
      std::vector<std::pair<std::string, std::string>> nu_pairs = {
	{"nu_e", "nu_e_bar"}, {"nu_mu", "nu_mu_bar"}, {"nu_tau", "nu_tau_bar"}
      };
      for (const auto& [nu, nubar] : nu_pairs)
	{
	  auto [a, b, sum] = get(cat, nu, nubar);
	  TString row = Form("%-15s %6d  %6d  %6d", (cat + " " + nu.substr(3)).c_str(), a, b, sum);
	  table->AddText(row);
	}
    }
    
    table->AddText("--------------------------------------------------");
    c1->cd();
    table->Draw();
    c1->Update();
    // DEBUG: Dump text content to console
    std::cout << "\n==== Summary Table (Screen Dump) ====\n";
    for (int i = 0; i < table->GetListOfLines()->GetEntries(); ++i) {
      auto* line = (TText*)table->GetListOfLines()->At(i);
      std::cout << line->GetTitle() << std::endl;
    }
    c1->Print(outname.Data(), "pdf("); 
    delete c1;
  }
  
  // === Page 2: Charm particles
  add_table_page("Charm Particle Type", charm_particle_by_flavor, flavors, std::string(outname.Data()), false, false);
  
  // === Page 3: Charm decay modes
  add_table_page("Charm Decay Modes", decay_mode_by_flavor, flavors, std::string(outname.Data()), false, true); // last page
  
  std::cout << "\n 3-page summary saved in " << outname << "\n";
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void AnalyzeMultiMuonCharmEvents() {
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();

  int total_multimu = 0;
  int numu_cc_charm_muonic = 0;
  int numubar_cc_charm_muonic = 0;
  int from_nue_nuebar = 0;
  int from_nc = 0;
  int other_cases = 0;

  for (Long64_t ientry = 0; ientry < entries; ++ientry) {
    data->GetEvent(ientry);

    // Multimuon tag requirement
    if (data->MuonTagStatus() < 2)
      continue;

    total_multimu++;

    int nu = data->NuType();
    int ccnc = data->CCNCType(); // 0 = NC, 1 = CC
    bool is_charm = (data->CharmParticleType() != 0);
    bool is_muonic_decay = data->CharmDecayM().find("Muonic") != std::string::npos;

    if (ccnc == 1 && is_charm && is_muonic_decay) {
      if (nu == 14)       numu_cc_charm_muonic++;
      else if (nu == -14) numubar_cc_charm_muonic++;
    }
    else if (nu == 12 || nu == -12) {
      from_nue_nuebar++;
    }
    else if (ccnc == 0) {
      from_nc++;
    }
    else {
      other_cases++;
    }
  }

  std::cout << "\n==== Multi-Muon Charm Analysis ====\n";
  std::cout << "Total multi-muon events:             " << total_multimu << std::endl;
  std::cout << "numu CC + charm + muonic decay:      " << numu_cc_charm_muonic << std::endl;
  std::cout << "numubar CC + charm + muonic decay:   " << numubar_cc_charm_muonic << std::endl;
  std::cout << "From nue or nuebar (any type):       " << from_nue_nuebar << std::endl;
  std::cout << "From NC interactions:                " << from_nc << std::endl;
  std::cout << "Other multi-muon events:             " << other_cases << std::endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Breakdown of "Other" Multi-Muon Events
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void InvestigateOtherMultiMuonEvents()
{
  FaserCalRecoPlot *data = filereader->tmCD(); 
  Long64_t entries = data->GetInputTree()->GetEntries();

  std::map<std::string, int> other_breakdown;
  std::map<std::string, int> decay_mode_by_nue_charm;

  // For 2D histogram: Muon count vs Neutrino type
  TH2F* h_muon_vs_flavor = new TH2F("h_muon_vs_flavor", "Muon Count vs Neutrino Flavor",
                                    6, 0, 6, 10, 0, 10); // x: 6 flavors, y: 0-10 muons
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(1, "nu_e");
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(2, "nubar_e");
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(3, "nu_mu");
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(4, "nubar_mu");
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(5, "nu_tau");
  h_muon_vs_flavor->GetXaxis()->SetBinLabel(6, "nubar_tau");
  h_muon_vs_flavor->GetYaxis()->SetTitle("# muons (MuTag)");

  auto pdg_label = [](int pdg) {
    if (pdg == 12) return "nu_e";
    if (pdg == -12) return "nubar_e";
    if (pdg == 14) return "nu_mu";
    if (pdg == -14) return "nubar_mu";
    if (pdg == 16) return "nu_tau";
    if (pdg == -16) return "nubar_tau";
    return "other";
  };

  auto flavor_index = [](int pdg) {
    if (pdg == 12) return 0;
    if (pdg == -12) return 1;
    if (pdg == 14) return 2;
    if (pdg == -14) return 3;
    if (pdg == 16) return 4;
    if (pdg == -16) return 5;
    return -1;
  };

  std::cout << "\n==== Investigating Other Multi-Muon Events ====\n";

  int dumped = 0;

  for (Long64_t ientry = 0; ientry < entries; ++ientry)
  {
    data->GetEvent(ientry);

    int nmu = data->MuonTagStatus();
    if (nmu < 2) continue;

    int pdg = data->NuType();
    int ccnc = data->CCNCType();
    bool is_charm = (data->CharmParticleType() != 0);
    std::string decay = data->CharmDecayM();
    std::string nu_str = pdg_label(pdg);

    // Fill histogram
    int binx = flavor_index(pdg) + 1;
    if (binx > 0 && binx <= 6) h_muon_vs_flavor->Fill(binx - 1, nmu);

    bool is_muonic = (decay.find("Muonic") != std::string::npos);

    bool is_nue = (pdg == 12 || pdg == -12);
    bool is_numucc = (pdg == 14 && ccnc == 1 && is_charm && is_muonic);
    bool is_numubcc = (pdg == -14 && ccnc == 1 && is_charm && is_muonic);
    bool is_nc = (ccnc == 0);

    if (is_nue)
      {
	std::string label = is_charm ? "nue(nubar)_with_charm" : "nue(nubar)_no_charm";
	other_breakdown[label]++;
	if (is_charm)
	  {
	    std::string mode = data->CharmDecayM();
	    decay_mode_by_nue_charm[mode]++;
	  }      
      }

    if (!is_numucc && !is_numubcc && !is_nue && !is_nc)
    {
      // Build breakdown key
      std::string key = nu_str;
      key += (ccnc == 1 ? " CC" : " NC");
      key += (is_charm ? " +charm" : " no-charm");
      if (is_charm)
        key += " decay: " + decay;
      other_breakdown[key]++;

      if (dumped < 10)
      {
        std::cout << "OtherEvent: Run " << data->RunNumber()
                  << ", Event " << data->EventNumber()
                  << ", nu=" << nu_str
                  << ", CCNC=" << ccnc
                  << ", Charm=" << is_charm
                  << ", Decay=" << decay
                  << ", Muons=" << nmu
                  << std::endl;
        dumped++;
      }
    }
  }

  // Print full breakdown
  std::cout << "\n==== Other Multi-Muon Breakdown ====\n";
  for (const auto& [key, count] : other_breakdown)
  {
    std::cout << key << ": " << count << std::endl;
  }
  // === Breakdown of decay modes for nue/nuebar + charm
  std::cout << "\n==== Charm Decay Modes in nue/nuebar Multi-Muon Events ====\n";
  for (const auto& [mode, count] : decay_mode_by_nue_charm)
    {
      std::cout << "  " << mode << ": " << count << "\n";
    }
  // Plot the 2D histogram
  TCanvas* c = new TCanvas("c2d", "Muon Count vs Flavor", 800, 600);
  h_muon_vs_flavor->Draw("COLZ TEXT");
  c->SaveAs("MuonCount_vs_Flavor.pdf");

  std::cout << "\n2D histogram saved as MuonCount_vs_Flavor.pdf\n";
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SaveCSV(const string& filename, const vector<string>& headers, const vector<vector<string>>& rows) {
  ofstream fout(filename);
  for (const auto& h : headers)
    fout << h << ",";
  fout << "\n";
  for (const auto& row : rows) {
    for (const auto& cell : row)
      fout << cell << ",";
    fout << "\n";
  }
  fout.close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
string classify(int pdg) {
  if (pdg == 12) return "nu_e";
  if (pdg == -12) return "nu_e_bar";
  if (pdg == 14) return "nu_mu";
  if (pdg == -14) return "nu_mu_bar";
  if (pdg == 16) return "nu_tau";
  if (pdg == -16) return "nu_tau_bar";
  return "other";
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GenerateSummaryWithMuSignsPDFandCSV(FaserCalRecoPlot* data, const string& outprefix = "Summary") {
  Long64_t entries = data->GetInputTree()->GetEntries();

  map<string, map<string, int>> summary, charm_particles, decay_modes;
  map<string, map<string, int>> mutag_signs, multimuon_counts;
  vector<string> flavors = {"nu_e", "nu_e_bar", "nu_mu", "nu_mu_bar", "nu_tau", "nu_tau_bar"};

  for (Long64_t i = 0; i < entries; ++i) {
    data->GetEvent(i);
    string nu = classify(data->NuType());
    if (nu == "other") continue;

    string cc = data->CCNCType() == 1 ? "CC" : "NC";
    bool charm = data->CharmParticleType() != 0;

    summary[cc][nu]++;
    if (charm) {
      summary["Charm"][nu]++;
      charm_particles[data->CharmParticleName()][nu]++;
      decay_modes[data->CharmDecayM()][nu]++;
    }

    // MuTag counts
    vector<int> signs = data->MuTagMuSign();
    for (int s : signs) {
      if (s > 0) mutag_signs["mu+"][nu]++;
      else if (s < 0) mutag_signs["mu-"][nu]++;
    }

    if (data->MuonTagStatus() >= 2)
      multimuon_counts["multi-muon"][nu]++;
  }

  TString pdffile = outprefix + ".pdf";
  TCanvas* c = new TCanvas("c", "Summary", 1000, 800);
  c->Print(pdffile + "[");

  auto make_table = [&](const string& title, map<string, map<string, int>>& table, const string& csvname) {
    TPaveText* pt = new TPaveText(0.01, 0.05, 0.99, 0.95);
    pt->SetTextFont(42);
    pt->SetTextSize(0.025);
    pt->SetFillColor(0);

    pt->AddText(title.c_str());
    string header = "Label";
    for (auto& f : flavors) header += " | " + f;
    header += " | Total";
    pt->AddText(header.c_str());

    vector<vector<string>> rows;
    vector<string> headers = {"Label"};
    headers.insert(headers.end(), flavors.begin(), flavors.end());
    headers.push_back("Total");

    for (auto& [label, countmap] : table) {
      vector<string> row = {label};
      int total = 0;
      for (auto& f : flavors) {
        int n = countmap[f];
        total += n;
        row.push_back(to_string(n));
      }
      row.push_back(to_string(total));
      pt->AddText((label + " : " + to_string(total)).c_str());
      rows.push_back(row);
    }

    c->cd();
    pt->Draw();
    c->Print(pdffile);
    SaveCSV(csvname, headers, rows);
    delete pt;
  };

  make_table("[1] CC/NC/Charm Summary", summary, outprefix + "_summary.csv");
  make_table("[2] Charm Particles", charm_particles, outprefix + "_charm_particles.csv");
  make_table("[3] Charm Decay Modes", decay_modes, outprefix + "_decay_modes.csv");
  make_table("[4] MuTag Signs (mu+/mu-)", mutag_signs, outprefix + "_mutag_signs.csv");
  make_table("[5] Multi-Muon Events", multimuon_counts, outprefix + "_multimuon.csv");

  c->Print(pdffile + "]");
  std::cout << "Summary tables saved in PDF and CSV format with prefix: " << outprefix << std::endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Multi-Muon Charm Summary with Muon Sign Breakdown
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetMultiMuonCharmSummaryWithSigns() {
  FaserCalRecoPlot* data = filereader->tmCD();
  Long64_t entries = data->GetInputTree()->GetEntries();

  int total_multimuon = 0;
  int numu_cc_charm_muonic = 0;
  int numubar_cc_charm_muonic = 0;
  int from_nue_or_nuebar = 0;
  int from_nc = 0;
  int other = 0;

  // Muon sign counters for each category
  std::map<std::string, std::pair<int, int>> muon_signs; // key -> <mu+, mu->
  auto update_muon_signs = [&](const std::string& label, const std::vector<int>& signs) {
    for (int s : signs) {
      if (s > 0) muon_signs[label].first++;
      else if (s < 0) muon_signs[label].second++;
    }
  };

  for (Long64_t i = 0; i < entries; ++i) {
    data->GetEvent(i);

    if (data->MuonTagStatus() < 2) continue;

    ++total_multimuon;
    const int nu = data->NuType();
    const int cc = data->CCNCType();
    const int charm = data->CharmParticleType();
    const std::string decay = data->CharmDecayM();
    const auto& signs = data->MuTagMuSign();

    if (nu == 14 && cc == 1 && charm != 0 && decay.find("Muonic") != std::string::npos) {
      ++numu_cc_charm_muonic;
      update_muon_signs("numu CC + charm + muonic", signs);
    }
    else if (nu == -14 && cc == 1 && charm != 0 && decay.find("Muonic") != std::string::npos) {
      ++numubar_cc_charm_muonic;
      update_muon_signs("numubar CC + charm + muonic", signs);
    }
    else if (std::abs(nu) == 12) {
      ++from_nue_or_nuebar;
      update_muon_signs("From nue or nuebar", signs);
    }
    else if (cc == 0) {
      ++from_nc;
      update_muon_signs("From NC interactions", signs);
    }
    else {
      ++other;
      update_muon_signs("Other multi-muon events", signs);
    }
  }

  // Print the summary table
  std::cout << "\n==== Multi-Muon Charm Analysis with Muon Signs ====\n";
  std::cout << "Total multi-muon events:             " << total_multimuon << std::endl;

  auto print_row = [&](const std::string& label, int count) {
    int mup = muon_signs[label].first;
    int mum = muon_signs[label].second;
    std::cout << std::setw(40) << std::left << label << ": "
              << std::setw(6) << count
              << "   (mu+: " << mup << ", mu-: " << mum << ")" << std::endl;
  };

  print_row("numu CC + charm + muonic", numu_cc_charm_muonic);
  print_row("numubar CC + charm + muonic", numubar_cc_charm_muonic);
  print_row("From nue or nuebar", from_nue_or_nuebar);
  print_row("From NC interactions", from_nc);
  print_row("Other multi-muon events", other);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
int main(int argc, char* argv[])
{
  if ( argc != 2 )
    {
      cout << "Usage: " << argv[0] << " filename.root" << endl;
      return -1;
    }
  
  string InputFile = argv[1];
  filereader = new TreeManager(InputFile);
  ReadFile();
  PlotCharmComparison();
  GetAllSummaryTablePDF();
  AnalyzeMultiMuonCharmEvents();
  InvestigateOtherMultiMuonEvents();
  GetMultiMuonCharmSummaryWithSigns();
  delete filereader;
  return 0;
}
