// Plot neutrino interaction vertices by detector region
// Extracts detector boundaries and tilt angle from GDML file
//
// Usage: root -l -b -q 'plotInteractionsByDetector_gdml.C("input.gfaser.root", "geometry.gdml")'
//

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <THStack.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <map>

// Structure to hold detector geometry information
struct DetectorGeometry {
  std::string name;
  double zmin;
  double zmax;
  int color;
};

// Parse GDML file to extract detector dimensions and tilt angle
bool parseGDML(const std::string& gdmlFile, std::vector<DetectorGeometry>& detectors, double& tilt_angle) {
  std::ifstream file(gdmlFile);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open GDML file " << gdmlFile << std::endl;
    return false;
  }
  
  std::string line;
  std::map<std::string, double> box_sizes;
  std::map<std::string, double> box_positions;
  std::map<std::string, std::string> volume_to_solid;
  
  // Read file line by line
  while (std::getline(file, line)) {
    // Look for box definitions with z dimension
    // Pattern: <box lunit="mm" name="NAME" x="..." y="..." z="VALUE"/>
    if (line.find("<box") != std::string::npos && line.find("lunit=\"mm\"") != std::string::npos) {
      size_t name_start = line.find("name=\"");
      size_t z_start = line.find("z=\"");
      if (name_start != std::string::npos && z_start != std::string::npos) {
        name_start += 6; // skip 'name="'
        size_t name_end = line.find("\"", name_start);
        std::string name = line.substr(name_start, name_end - name_start);
        
        z_start += 3; // skip 'z="'
        size_t z_end = line.find("\"", z_start);
        std::string z_str = line.substr(z_start, z_end - z_start);
        double z_size = std::stod(z_str);
        box_sizes[name] = z_size;
      }
    }
    
    // Look for physical volume positions
    // Pattern: <position name="NAME_pos" unit="mm" x="..." y="..." z="VALUE"/>
    if (line.find("<position") != std::string::npos && line.find("_pos\"") != std::string::npos) {
      size_t name_start = line.find("name=\"");
      size_t z_start = line.find("z=\"");
      if (name_start != std::string::npos && z_start != std::string::npos) {
        name_start += 6; // skip 'name="'
        size_t name_end = line.find("_pos\"", name_start);
        if (name_end != std::string::npos) {
          std::string name = line.substr(name_start, name_end - name_start);
          
          z_start += 3; // skip 'z="'
          size_t z_end = line.find("\"", z_start);
          std::string z_str = line.substr(z_start, z_end - z_start);
          double z_pos = std::stod(z_str);
          box_positions[name] = z_pos;
        }
      }
    }
    
    // Look for detector assembly rotation
    // Pattern: <rotation name="DetectorAssemblyPV..._rot" unit="deg" x="0" y="VALUE" z="0"/>
    if (line.find("<rotation") != std::string::npos && 
        line.find("DetectorAssemblyPV") != std::string::npos && 
        line.find("_rot\"") != std::string::npos) {
      size_t y_start = line.find("y=\"");
      if (y_start != std::string::npos) {
        y_start += 3; // skip 'y="'
        size_t y_end = line.find("\"", y_start);
        std::string y_str = line.substr(y_start, y_end - y_start);
        tilt_angle = std::stod(y_str);
      }
    }
  }
  
  file.close();
  
  // Calculate detector boundaries
  // 3DCAL: ContainerBox centered at origin
  double z_3dcal = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerBox") != std::string::npos) {
      z_3dcal = kv.second;
      break;
    }
  }
  
  // ECAL: ContainerEcal at specific position
  double z_ecal_size = -999;
  double z_ecal_pos = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerEcal") != std::string::npos) {
      z_ecal_size = kv.second;
      break;
    }
  }
  for (const auto& kv : box_positions) {
    if (kv.first.find("rearCal") != std::string::npos && 
        kv.first.find("Abs") == std::string::npos && 
        kv.first.find("Scint") == std::string::npos) {
      z_ecal_pos = kv.second;
      break;
    }
  }
  
  // AHCAL: ContainerHcal at specific position
  double z_hcal_size = -999;
  double z_hcal_pos = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("ContainerHcal") != std::string::npos) {
      z_hcal_size = kv.second;
      break;
    }
  }
  for (const auto& kv : box_positions) {
    if (kv.first.find("rearHCal") != std::string::npos) {
      z_hcal_pos = kv.second;
      break;
    }
  }
  
  // Muon Spectrometer
  double z_muon_size = -999;
  double z_muon_pos = -999;
  for (const auto& kv : box_sizes) {
    if (kv.first.find("MuonSpectrometer") != std::string::npos) {
      z_muon_size = kv.second;
      break;
    }
  }
  for (const auto& kv : box_positions) {
    if (kv.first.find("MuonSpectrometer") != std::string::npos) {
      z_muon_pos = kv.second;
      break;
    }
  }
  
  // Validate and construct detector list
  if (z_3dcal > 0) {
    DetectorGeometry det;
    det.name = "3DCAL";
    det.zmin = -z_3dcal / 2.0;
    det.zmax = z_3dcal / 2.0;
    det.color = kBlue;
    detectors.push_back(det);
    std::cout << "Found 3DCAL: " << det.zmin << " to " << det.zmax << " mm" << std::endl;
  }
  
  if (z_ecal_size > 0 && z_ecal_pos != -999) {
    DetectorGeometry det;
    det.name = "ECAL";
    det.zmin = z_ecal_pos - z_ecal_size / 2.0;
    det.zmax = z_ecal_pos + z_ecal_size / 2.0;
    det.color = kRed;
    detectors.push_back(det);
    std::cout << "Found ECAL: " << det.zmin << " to " << det.zmax << " mm" << std::endl;
  }
  
  if (z_hcal_size > 0 && z_hcal_pos != -999) {
    DetectorGeometry det;
    det.name = "AHCAL";
    det.zmin = z_hcal_pos - z_hcal_size / 2.0;
    det.zmax = z_hcal_pos + z_hcal_size / 2.0;
    det.color = kGreen+2;
    detectors.push_back(det);
    std::cout << "Found AHCAL: " << det.zmin << " to " << det.zmax << " mm" << std::endl;
  }
  
  if (z_muon_size > 0 && z_muon_pos != -999) {
    DetectorGeometry det;
    det.name = "MuonSpec";
    det.zmin = z_muon_pos - z_muon_size / 2.0;
    det.zmax = z_muon_pos + z_muon_size / 2.0;
    det.color = kMagenta;
    detectors.push_back(det);
    std::cout << "Found MuonSpec: " << det.zmin << " to " << det.zmax << " mm" << std::endl;
  }
  
  std::cout << "Tilt angle: " << tilt_angle << " degrees" << std::endl;
  
  return detectors.size() > 0;
}

// Apply tilt correction to vertex coordinates
void apply_tilt_correction(double &x, double &z, double tilt_deg) {
  double tilt_rad = tilt_deg * M_PI / 180.0;
  double cos_tilt = cos(tilt_rad);
  double sin_tilt = sin(tilt_rad);
  
  // Rotate coordinates: tilt around Y-axis
  double x_new = x * cos_tilt + z * sin_tilt;
  double z_new = -x * sin_tilt + z * cos_tilt;
  
  x = x_new;
  z = z_new;
}

void plotInteractionsByDetector_gdml(const char* inputFile, const char* gdmlFile) {
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  
  // Parse GDML file to get detector geometry
  std::vector<DetectorGeometry> detectors;
  double tilt_angle = 0.0;
  
  std::cout << "Parsing GDML file: " << gdmlFile << std::endl;
  if (!parseGDML(gdmlFile, detectors, tilt_angle)) {
    std::cerr << "Error: Failed to parse GDML file or no detectors found" << std::endl;
    return;
  }
  
  // Open input file
  TFile *file = TFile::Open(inputFile);
  if (!file || !file->IsOpen()) {
    std::cerr << "Error: Cannot open file " << inputFile << std::endl;
    return;
  }
  
  TTree *tree = (TTree*)file->Get("gFaser");
  if (!tree) {
    std::cerr << "Error: Cannot find tree gFaser" << std::endl;
    return;
  }
  
  // Set up branches
  Double_t vx, vy, vz;
  tree->SetBranchAddress("vx", &vx);
  tree->SetBranchAddress("vy", &vy);
  tree->SetBranchAddress("vz", &vz);
  
  // Create histograms
  TH1D *h_z_all = new TH1D("h_z_all", "All Interactions;Z [mm];Events", 400, -2000, 5000);
  TH1D *h_x_all = new TH1D("h_x_all", "All Interactions;X [mm];Events", 200, -1000, 1000);
  TH1D *h_y_all = new TH1D("h_y_all", "All Interactions;Y [mm];Events", 200, -1000, 1000);
  TH2D *h_xy_all = new TH2D("h_xy_all", "All Interactions;X [mm];Y [mm]", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *h_xz_all = new TH2D("h_xz_all", "All Interactions;Z [mm];X [mm]", 200, -2000, 5000, 100, -1000, 1000);
  TH2D *h_yz_all = new TH2D("h_yz_all", "All Interactions;Z [mm];Y [mm]", 200, -2000, 5000, 100, -1000, 1000);
  
  // Histograms for each detector
  std::map<std::string, TH1D*> h_z_detector;
  std::map<std::string, TH1D*> h_x_detector;
  std::map<std::string, TH1D*> h_y_detector;
  std::map<std::string, TH2D*> h_xy_detector;
  std::map<std::string, TH2D*> h_xz_detector;
  std::map<std::string, TH2D*> h_yz_detector;
  
  for (const auto& det : detectors) {
    h_z_detector[det.name] = new TH1D(Form("h_z_%s", det.name.c_str()), 
                                       Form("%s Interactions;Z [mm];Events", det.name.c_str()), 
                                       400, -2000, 5000);
    h_x_detector[det.name] = new TH1D(Form("h_x_%s", det.name.c_str()), 
                                       Form("%s Interactions;X [mm];Events", det.name.c_str()), 
                                       200, -1000, 1000);
    h_y_detector[det.name] = new TH1D(Form("h_y_%s", det.name.c_str()), 
                                       Form("%s Interactions;Y [mm];Events", det.name.c_str()), 
                                       200, -1000, 1000);
    h_xy_detector[det.name] = new TH2D(Form("h_xy_%s", det.name.c_str()), 
                                        Form("%s Interactions;X [mm];Y [mm]", det.name.c_str()), 
                                        100, -1000, 1000, 100, -1000, 1000);
    h_xz_detector[det.name] = new TH2D(Form("h_xz_%s", det.name.c_str()), 
                                        Form("%s Interactions;Z [mm];X [mm]", det.name.c_str()), 
                                        200, -2000, 5000, 100, -1000, 1000);
    h_yz_detector[det.name] = new TH2D(Form("h_yz_%s", det.name.c_str()), 
                                        Form("%s Interactions;Z [mm];Y [mm]", det.name.c_str()), 
                                        200, -2000, 5000, 100, -1000, 1000);
    h_z_detector[det.name]->SetLineColor(det.color);
    h_z_detector[det.name]->SetLineWidth(2);
  }
  
  // Event counters
  Long64_t total_events = 0;
  std::map<std::string, Long64_t> events_per_detector;
  for (const auto& det : detectors) {
    events_per_detector[det.name] = 0;
  }
  
  // Loop over events
  Long64_t nentries = tree->GetEntries();
  std::cout << "\nProcessing " << nentries << " events..." << std::endl;
  
  for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    
    if (i % 10000 == 0) {
      std::cout << "Processing event " << i << " / " << nentries << std::endl;
    }
    
    // Convert to mm and apply tilt correction
    double vtx_x = vx * 1e3;
    double vtx_y = vy * 1e3;
    double vtx_z = vz * 1e3;
    
    if (tilt_angle != 0.0) {
      apply_tilt_correction(vtx_x, vtx_z, tilt_angle);
    }
    
    // Fill all interactions
    h_z_all->Fill(vtx_z);
    h_x_all->Fill(vtx_x);
    h_y_all->Fill(vtx_y);
    h_xy_all->Fill(vtx_x, vtx_y);
    h_xz_all->Fill(vtx_z, vtx_x);
    h_yz_all->Fill(vtx_z, vtx_y);
    total_events++;
    
    // Check which detector it belongs to
    for (const auto& det : detectors) {
      if (vtx_z >= det.zmin && vtx_z <= det.zmax) {
        h_z_detector[det.name]->Fill(vtx_z);
        h_x_detector[det.name]->Fill(vtx_x);
        h_y_detector[det.name]->Fill(vtx_y);
        h_xy_detector[det.name]->Fill(vtx_x, vtx_y);
        h_xz_detector[det.name]->Fill(vtx_z, vtx_x);
        h_yz_detector[det.name]->Fill(vtx_z, vtx_y);
        events_per_detector[det.name]++;
        break; // Each event can only be in one detector
      }
    }
  }
  
  // Print statistics
  std::cout << "\n=== Event Statistics ===" << std::endl;
  std::cout << "Total events: " << total_events << std::endl;
  for (const auto& det : detectors) {
    double fraction = 100.0 * events_per_detector[det.name] / total_events;
    std::cout << det.name << ": " << events_per_detector[det.name] 
              << " (" << fraction << "%)" << std::endl;
  }
  
  // Create output file
  TFile *outfile = new TFile("interaction_plots.root", "RECREATE");
  
  // Plot 1: Z distribution - All vs individual detectors
  TCanvas *c1 = new TCanvas("c1", "Z Distribution", 1200, 800);
  h_z_all->SetLineColor(kBlack);
  h_z_all->SetLineWidth(2);
  h_z_all->Draw();
  
  TLegend *leg1 = new TLegend(0.7, 0.5, 0.88, 0.88);
  leg1->AddEntry(h_z_all, Form("All (%lld)", total_events), "l");
  for (const auto& det : detectors) {
    h_z_detector[det.name]->Draw("same");
    leg1->AddEntry(h_z_detector[det.name], 
                   Form("%s (%lld)", det.name.c_str(), events_per_detector[det.name]), "l");
  }
  leg1->Draw();
  
  // Add text box with event statistics
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  double y_pos = 0.88;
  for (const auto& det : detectors) {
    latex.DrawLatex(0.12, y_pos, Form("%s: %lld events", det.name.c_str(), events_per_detector[det.name]));
    y_pos -= 0.03;
  }
  latex.DrawLatex(0.12, y_pos, Form("Total: %lld events", total_events));
  
  c1->Write();
  c1->SaveAs("z_distribution_all.png");
  
  // Plot 2: Z distribution - stacked by detector
  TCanvas *c2 = new TCanvas("c2", "Z Distribution by Detector", 1200, 800);
  THStack *hs_z = new THStack("hs_z", "Interaction Vertices by Detector;Z [mm];Events");
  for (const auto& det : detectors) {
    h_z_detector[det.name]->SetFillColor(det.color);
    h_z_detector[det.name]->SetFillStyle(1001);
    hs_z->Add(h_z_detector[det.name]);
  }
  hs_z->Draw();
  
  TLegend *leg2 = new TLegend(0.7, 0.5, 0.88, 0.88);
  for (const auto& det : detectors) {
    leg2->AddEntry(h_z_detector[det.name], 
                   Form("%s (%lld)", det.name.c_str(), events_per_detector[det.name]), "f");
  }
  leg2->Draw();
  c2->Write();
  c2->SaveAs("z_distribution_stacked.png");
  
  // Plot 3: 2D plots - All interactions
  TCanvas *c3 = new TCanvas("c3", "2D Distributions - All", 1800, 600);
  c3->Divide(3, 1);
  c3->cd(1);
  h_xy_all->Draw("colz");
  c3->cd(2);
  h_xz_all->Draw("colz");
  c3->cd(3);
  h_yz_all->Draw("colz");
  c3->Write();
  c3->SaveAs("2d_distributions_all.png");
  
  // Plot 4+: Individual detector plots (only for 3DCAL, ECAL, AHCAL)
  for (const auto& det : detectors) {
    if (det.name == "MuonSpec") continue; // Skip muon spec for individual plots
    
    TCanvas *c_det = new TCanvas(Form("c_%s", det.name.c_str()), 
                                  Form("%s Distributions", det.name.c_str()), 
                                  1800, 1200);
    c_det->Divide(3, 2);
    
    c_det->cd(1);
    h_x_detector[det.name]->Draw();
    
    c_det->cd(2);
    h_y_detector[det.name]->Draw();
    
    c_det->cd(3);
    h_z_detector[det.name]->Draw();
    
    c_det->cd(4);
    h_xy_detector[det.name]->Draw("colz");
    
    c_det->cd(5);
    h_xz_detector[det.name]->Draw("colz");
    
    c_det->cd(6);
    h_yz_detector[det.name]->Draw("colz");
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    c_det->cd(1);
    latex.DrawLatex(0.15, 0.85, Form("Z range: %.0f to %.0f mm", det.zmin, det.zmax));
    latex.DrawLatex(0.15, 0.80, Form("Events: %lld", events_per_detector[det.name]));
    
    c_det->Write();
    c_det->SaveAs(Form("%s_distributions.png", det.name.c_str()));
  }
  
  // Write all histograms
  h_z_all->Write();
  h_x_all->Write();
  h_y_all->Write();
  h_xy_all->Write();
  h_xz_all->Write();
  h_yz_all->Write();
  
  for (const auto& det : detectors) {
    h_z_detector[det.name]->Write();
    h_x_detector[det.name]->Write();
    h_y_detector[det.name]->Write();
    h_xy_detector[det.name]->Write();
    h_xz_detector[det.name]->Write();
    h_yz_detector[det.name]->Write();
  }
  
  outfile->Close();
  file->Close();
  
  std::cout << "\nPlots saved to interaction_plots.root" << std::endl;
  std::cout << "PNG files created for each detector" << std::endl;
}
