#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TGraph2D.h>
#include <TDatabasePDG.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

#include "../CoreUtils/TPOEvent.hh"
// Directly include the implementation to avoid unresolved symbols if dictionary/library not preloaded.
#include "../CoreUtils/TPOEvent.cc"

// Auto-detect a TTree that contains a TPOEvent branch (named either "event" or something similar)
static TTree* FindPOEventTree(TFile* f, std::string& branchName) {
    if(!f) return nullptr;
    TIter keyIter(f->GetListOfKeys());
    while(auto *key = (TKey*)keyIter()) {
        if(strcmp(key->GetClassName(), "TTree")!=0) continue;
        TTree* t = (TTree*)key->ReadObj();
        // Try common branch names
        std::vector<std::string> candidates = {"event", "poevent", "POEvent", "TPOEvent"};
        for(const auto& b : candidates) {
            if(t->GetBranch(b.c_str())) {
                TBranch* br = t->GetBranch(b.c_str());
                if(br && br->GetClassName() && std::string(br->GetClassName()).find("TPOEvent")!=std::string::npos) {
                    branchName = b; return t;
                }
            }
        }
        // Fallback: scan all branches for class type containing TPOEvent
        TObjArray* bl = t->GetListOfBranches();
        for(int i=0;i<bl->GetEntries();++i) {
            TBranch* br = (TBranch*)bl->At(i);
            if(br && br->GetClassName() && std::string(br->GetClassName()).find("TPOEvent")!=std::string::npos) {
                branchName = br->GetName(); return t;
            }
        }
    }
    return nullptr;
}

void plotPO(const char* filename = "FASERMC-PO-Run500-0_99261.root") {
    // Load required ROOT libs
    gSystem->Load("libCore");
    gSystem->Load("libTree");

    // Implementation already included above; no need for ACLiC compile step.

    TFile* f = TFile::Open(filename, "READ");
    if(!f || f->IsZombie()) {
        std::cerr << "Cannot open ROOT file: " << filename << std::endl;
        return;
    }

    std::string branchName;
    TTree* tree = FindPOEventTree(f, branchName);
    if(!tree) {
        std::cerr << "No TTree with a TPOEvent branch found in file." << std::endl;
        f->Close();
        return;
    }
    std::cout << "Using tree: " << tree->GetName() << " branch: " << branchName << std::endl;

    TPOEvent* event = nullptr;
    tree->SetBranchAddress(branchName.c_str(), &event);

    TH1D* h_nParticles   = new TH1D("h_nParticles", "Number of generator particles;N_{PO};Events", 100, 0, 200);
    TH1D* h_visE         = new TH1D("h_visE", "Visible energy;E_{vis} (GeV);Events", 100, 0, 5000);
    TH1D* h_leptonP      = new TH1D("h_leptonP", "Outgoing lepton momentum;|p_{l}| (GeV);Events", 100, 0, 5000);
    TH1D* h_ptmiss       = new TH1D("h_ptmiss", "Missing transverse momentum;p_{T}^{miss} (GeV);Events", 100, 0, 30);
    TH1D* h_leptonPDG    = new TH1D("h_leptonPDG", "Outgoing lepton PDG;PDG;Events", 61, -30.5, 30.5);
    TH2D* h_leptonAngleP = new TH2D("h_leptonAngleP", "Lepton polar angle vs momentum;|p_{l}| (GeV);#theta (rad)", 100, 0, 5000, 90, 0, 0.9);

    // Vertex position histograms (in mm)
    TH1D* h_vx = new TH1D("h_vx", "Primary vertex X;X (mm);Events", 120, -600, 600);
    TH1D* h_vy = new TH1D("h_vy", "Primary vertex Y;Y (mm);Events", 120, -600, 600);
    TH1D* h_vz = new TH1D("h_vz", "Primary vertex Z;Z (mm);Events", 120, -1500, 1500);
    TH2D* h_vxy = new TH2D("h_vxy", "Primary vertex X-Y;X (mm);Y (mm)", 120, -600, 600, 120, -600, 600);
    TH2D* h_vxz = new TH2D("h_vxz", "Primary vertex X-Z;X (mm);Z (mm)", 120, -600, 600, 120, -1500, 1500);
    TH2D* h_vyz = new TH2D("h_vyz", "Primary vertex Y-Z;Y (mm);Z (mm)", 120, -600, 600, 120, -1500, 1500);
    //h_vx->SetLineWidth(2);
    //h_vy->SetLineWidth(2);
    //h_vz->SetLineWidth(2);
    h_vxy->SetStats(0);
    h_vxz->SetStats(0);
    h_vyz->SetStats(0);
    //TGraph2D* g_vxyz = new TGraph2D();
    //g_vxyz->SetName("g_vxyz");
    //g_vxyz->SetTitle("Primary vertex in 3D;X (mm);Y (mm);Z (mm)");

    // Incoming neutrino energy (all reactions)
    TH1D* h_nuE_all = new TH1D("h_nuE_all", "Incoming neutrino energy (all);E_{#nu} (GeV);Events", 100, 0, 5000);
    // Per-interaction type neutrino energy histograms
    std::vector<std::string> reactionTypes = {
        "nueCC","numuCC","nutauCC","nuNC","nuES",
        "antinueCC","antinumuCC","antinutauCC"
    };
    std::map<std::string, TH1D*> h_nuE_type;
    for(const auto& r : reactionTypes) {
        std::string name = std::string("h_nuE_") + r;
        std::string title = std::string("Incoming neutrino energy ") + r + ";E_{#nu} (GeV);Events";
        h_nuE_type[r] = new TH1D(name.c_str(), title.c_str(), 100, 0, 5000);
    }

    Long64_t nEntries = tree->GetEntries();
    std::cout << "Entries: " << nEntries << std::endl;

    for(Long64_t i=0;i<nEntries;++i) {
        tree->GetEntry(i);
        if(!event) continue;
        event->kinematics_event();
    // Fill vertex positions
    double vx = event->prim_vx.x();
    double vy = event->prim_vx.y();
    double vz = event->prim_vx.z();
    h_vx->Fill(vx);
    h_vy->Fill(vy);
    h_vz->Fill(vz);
    h_vxy->Fill(vx, vy);
    h_vxz->Fill(vx, vz);
    h_vyz->Fill(vy, vz);
    //g_vxyz->SetPoint(g_vxyz->GetN(), vx, vy, vz);

        h_nParticles->Fill(event->n_particles());
        h_visE->Fill(event->Evis);
        h_ptmiss->Fill(event->ptmiss);

        // Incoming neutrino energy (convert MeV to GeV if stored that way)
        double Enu = event->in_neutrino.m_energy;
        std::cout << "Event " << i << " Neutrino energy: " << Enu << " GeV" << std::endl;
        h_nuE_all->Fill(Enu);
        std::string reac = event->isES() ? std::string("nuES") : std::string(event->reaction_desc());
        auto it = h_nuE_type.find(reac);
        if(it != h_nuE_type.end()) it->second->Fill(Enu);

        double plep = sqrt(event->out_lepton.m_px*event->out_lepton.m_px +
                           event->out_lepton.m_py*event->out_lepton.m_py +
                           event->out_lepton.m_pz*event->out_lepton.m_pz);
        h_leptonP->Fill(plep);
        h_leptonPDG->Fill(event->out_lepton.m_pdg_id);
        double theta = atan2(sqrt(event->out_lepton.m_px*event->out_lepton.m_px + event->out_lepton.m_py*event->out_lepton.m_py), fabs(event->out_lepton.m_pz));
        h_leptonAngleP->Fill(plep, theta);
    }

    TCanvas* c1 = new TCanvas("c1", "POEvent Summary", 1800, 1000);
    c1->Divide(4,2); // expand grid to include neutrino energy summary
    c1->cd(1); h_nParticles->Draw();
    c1->cd(2); h_visE->Draw();
    c1->cd(3); h_leptonP->Draw();
    c1->cd(4); h_ptmiss->Draw();
    c1->cd(5); h_leptonPDG->Draw();
    c1->cd(6); h_leptonAngleP->Draw("COLZ");
    c1->cd(7); h_nuE_all->Draw();
    // leave pad 8 for legend of interaction types (optional cumulative)
    c1->cd(8);
    // Find maximum bin content among per-type histograms to set a common y-scale
    double maxY = 0.0;
    for(const auto& r : reactionTypes) {
        double m = h_nuE_type[r]->GetMaximum();
        if(m > maxY) maxY = m;
    }
    if(maxY <= 0) maxY = 1.0;
    for(const auto& r : reactionTypes) h_nuE_type[r]->SetMaximum(maxY*1.10);

    int idx=0;  
    for(const auto& r : reactionTypes) {
        h_nuE_type[r]->SetLineColor(idx+1);
        if(idx==0){
            h_nuE_type[r]->Draw();
            h_nuE_type[r]->SetTitle("Incoming neutrino energy by interaction type;E_{#nu} (GeV);Events");
        }
        else
        h_nuE_type[r]->Draw("same");
        
        idx++;
    }
    c1->SaveAs("POEventSummary.pdf");

    // Separate canvas for per-reaction neutrino energy plots
    TCanvas* cNu = new TCanvas("cNu", "Neutrino Energy by Reaction", 1800, 1000);
    cNu->Divide(3,3);
    int pad=1;
    for(const auto& r : reactionTypes) {
        cNu->cd(pad);
        h_nuE_type[r]->SetLineColor(pad+1);
        h_nuE_type[r]->Draw();
        pad++;
    }
    cNu->SaveAs("POEventNuE.pdf");

    // Vertex summary canvas
    TCanvas* cV = new TCanvas("cV", "Primary Vertex", 1800, 1000);
    cV->Divide(3,2);
    cV->cd(1); h_vx->Draw();
    cV->cd(2); h_vy->Draw();
    cV->cd(3); h_vz->Draw();
    cV->cd(4); h_vxy->Draw("COLZ");
    cV->cd(5); h_vxz->Draw("COLZ");
    cV->cd(6); h_vyz->Draw("COLZ");
    cV->SaveAs("POEventVertex.pdf");

    // 3D scatter of vertex positions
    /*
    TCanvas* cV3 = new TCanvas("cV3", "Primary Vertex 3D", 1200, 900);
    g_vxyz->Draw("P0");
    cV3->SaveAs("POEventVertex3D.pdf");
    */
    TFile out("POEventSummary.root", "RECREATE");
    h_nParticles->Write(); h_visE->Write(); h_leptonP->Write(); h_ptmiss->Write(); h_leptonPDG->Write(); h_leptonAngleP->Write();
    h_nuE_all->Write();
    for(const auto& r : reactionTypes) h_nuE_type[r]->Write();
    h_vx->Write(); h_vy->Write(); h_vz->Write();
    h_vxy->Write(); h_vxz->Write(); h_vyz->Write();
    //g_vxyz->Write();
    out.Close();

    f->Close();
    std::cout << "Finished. Output: POEventSummary.pdf, POEventNuE.pdf, POEventSummary.root" << std::endl;
}
