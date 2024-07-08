#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <TGeoManager.h>

#include "TcalEvent.hh"
#include "TPORecoEvent.hh"

void load_geometry() {
    // Load the GDML geometry
    TGeoManager::Import("../GeomGDML/geometry.gdml");

    // Draw the geometry
//    gGeoManager->GetTopVolume()->Draw("ogl");
}

int main(int argc, char** argv) {

    load_geometry();

    std::string base_path = "input/tcalevent_";

    TFile *m_rootFile = new TFile("batchreco.root", "RECREATE", "", 0); // last is the compression level
    if (!m_rootFile || !m_rootFile->IsOpen())
    {
        throw std::runtime_error("Could not create ROOT file");
    }
    m_rootFile->cd();

    TH1D e_em_energy = TH1D("e_em_energy", "electrons: energy fraction", 200, -1., 1.);
    TH2D e_em_energy2 = TH2D("e_em_energy2", "electrons: energy fraction", 100, 0., 1000., 200, -1., 1.);
    TH1D pi_em_energy = TH1D("pi_em_energy", "pions: energy fraction", 100, -1., 1.);
    TH2D pi_em_energy2 = TH2D("pi_em_energy2", "pions : energy fraction vs E", 100, 0.,200.,100,-1.,1.);
    TH1D p_em_energy = TH1D("p_em_energy", "protons: energy fraction", 100, -1., 1.);
    TH2D p_em_energy2 = TH2D("p_em_energy2", "protons : energy fraction vs E", 100, 0.,200.,100,-1.,1.);

    int ievent = 0;
    int error = 0;

    while (error == 0 && ievent<99999999) {

        // Create an instance of TcalEvent and TPOEvent
        TcalEvent *fTcalEvent = new TcalEvent();
        TPOEvent *POevent = new TPOEvent();

        error = fTcalEvent -> Load_event(base_path, ievent++, POevent);
        if(error != 0) break;
    
        std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
        std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fTotalLength << " mm " << std::endl;
	    std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
        std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;

        std::cout << " copied digitized tracks " << fTcalEvent->fTracks.size() << std::endl;

        fTcalEvent -> fTPOEvent -> dump_event();

        TPORecoEvent* fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
        fPORecoEvent -> Reconstruct();
        fPORecoEvent -> Dump(); 

        for(auto it : fPORecoEvent -> fPORecs) {
            int POID = it->POID;
            struct PO* aPO = &fTcalEvent->fTPOEvent->POs[POID];
            int PDG = aPO->m_pdg_id;
            if(abs(PDG) == 11) {
                double POEne = aPO->m_energy;
                double RecoEne = it->fTotal.Ecompensated;
                double f = (RecoEne-POEne)/POEne;
                e_em_energy.Fill(f);
                e_em_energy2.Fill(POEne, f);
            } else if(abs(PDG) == 111 || abs(PDG) == 211) {
                double POEne = aPO->m_energy;
                double RecoEne = it->fTotal.Ecompensated;
                double f = (RecoEne-POEne)/POEne;
                pi_em_energy.Fill(f);
                pi_em_energy2.Fill(POEne, f);
            } else if(abs(PDG) == 2112) {
                double POEne = aPO->m_energy;
                double RecoEne = it->fTotal.Ecompensated;
                double f = (RecoEne-POEne)/POEne;
                p_em_energy.Fill(f);
                p_em_energy2.Fill(POEne, f);
            }
        };

        delete fPORecoEvent;
        delete POevent;
        delete fTcalEvent;       
    }

    m_rootFile->Write();
    m_rootFile->Close();

    return 0;
}
