#ifndef __FaserCal_display_FaserCalDisplay__
#define __FaserCal_display_FaserCalDisplay__

#include <TVector3.h>
#include <RQ_OBJECT.h>
#include <TQObject.h>
#include <string>
#include <TROOT.h>
#include <TCanvas.h>
#include <TRootEmbeddedCanvas.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH1.h>
#include <TH2I.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>
#include <TObject.h> 
// libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveProjections.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveElement.h>
#include <TEveRGBAPalette.h>
#include <TEveRGBAPaletteOverlay.h>
#include <TEveTrans.h>
#include <TEveBoxSet.h>
#include <TEveArrow.h>
#include <TEveText.h>
#include <TEveSelection.h>
#include <TEveCompound.h>
#include <TEveDigitSet.h>
#include <TEveQuadSet.h>
#include <TEveLine.h>

//
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeometry.h>
#include <TGDMLParse.h>
#include <TGeoMatrix.h>
#include "TGeoSphere.h"
#include "TGeoBBox.h"
#include "TGeoMedium.h"
#include <TView3D.h>

// libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
//#include <TGRadioButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGFileDialog.h>
#include <TGLFontManager.h>
#include <TGNumberEntry.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
#include <TGLViewer.h>
#include <TGLOrthoCamera.h>
#include <TGLCamera.h>
#include <TGProgressBar.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>
#include <TGStatusBar.h>
#include <TGLLightSet.h>
//libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>

#include <TPolyLine3D.h>
#include <TLegend.h>
#include <TGHtml.h>
#include "TPaveText.h"
#include <TText.h>

#include "../CoreUtils/TcalEvent.hh" 
#include "../CoreUtils/TPORecoEvent.hh"
#include "../CoreUtils/TPOEvent.hh"
//#include "fastjet/ClusterSequence.hh"

#include <vector>
#include <string>
#include <map>
#include <set>


struct HitInfo : public TObject
{
  enum EType {
    kCalVoxel,
    kPixelHit,
    kMuTagHit,
    kMuonSpectHit
  };

  EType     fType;
  int       fPDG      = 0;
  int       fTrackID  = -1;
  int       fParentID = -1;
  int       fPrimaryID= -1;
  long      fChannelID= -1;
  double    fEDep     = 0.0;
  double    fX = 0.0, fY = 0.0, fZ = 0.0;

  HitInfo() = default;
};

namespace display 
{
    class FaserCalDisplay
    {
        RQ_OBJECT("FaserCalDisplay")
    public:
        FaserCalDisplay();
      virtual ~FaserCalDisplay();
      // Declare a map to store the association between TGeoVolume and track information
      std::unordered_map<TGeoVolume*, DigitizedTrack*> volumeTrackMap;

      //const std::vector<fastjet::PseudoJet>& GetJets() const { return fJets; }
      //fastjet::PseudoJet computeMomentumFromVoxel(ROOT::Math::XYZVector position, double totalEnergy);

    std::vector<std::unique_ptr<HitInfo>> fHitInfos;  // Store hit information for selection

      void EventDisplay();
      void GetEventDisplay();
      void GetDetector();
      ROOT::Math::XYZVector DetToWorld(const ROOT::Math::XYZVector &detPos) const;
      
      void InitGeometryPointers();

      void DoExit();
      
      void LoadEvent(int iRun, int iEvent, int iMask );
      void DrawEvent(int iRun, int iEvent, int iMask);
      void ZoomEvent(int iEvent);
      
      Bool_t CheckEventInFile(size_t ievent);
      void SetEventNumber();
      void ShowEvent();
      void ZoomingEvent();
      void NextEvent();
      void PreviousEvent();
      void Isolate();
      void DumpEvent();
      void DoSlider(Int_t position);
      TGStatusBar *fStatusBar;

      void ShowAxis();
      void LoadAllEvents();
      void ShowPrimary();
      void ShowPixelHits();
      void ShowMuonParticles();
      void ShowProtonParticles();
      void ShowPionParticles();
      void ShowKaonParticles();
      void ShowSecondaryShowers();
      void ShowSecondaryHadShowers();
      void ShowPixelRecoTrack();
      void ShowShortLivedParticle();
      void JetReconstructions();
      void ShowJetHits();
      void DrawMCTruthVertexPoint();

      void ShowClusterHits();
      void ShowRecoVoxelHits();

      void GetClusterInfo();
      void Get3DClusterInfo();
      void GetReconstructedVoxels(TGeoShape *bigbox, TGeoMedium *air, TGeoShape *box);
      void ShowOnlyGhostHits();
      void GetRearECAL();

      Double_t fVisibleEnergy = 0.0;
      Double_t fTotalEnergy = 0.0;
      Double_t fRearECalEnergy = 0.0;
      Double_t fRearHCalEnergy = 0.0;
      Double_t fRearMuCalEnergy = 0.0;

      void CountHitsInCube();
      void PlotParticleFractionsAndHitsInCubes();

      void SetMyStyle();
      void BackgroundColor();

      void EnablePicking();
      void OnPick(TGLPhysicalShape *shape, TObject *obj, Int_t event);
      void OnShapeSelected(TEveElement* selectedElement);
      
      void PlotHistogramsFromRootFile();
      void CleanViewer();
      void CleanCanvas();
      static TCanvas * CreateCanvas(const char *plot_name, int tb);
      TCanvas *CreateTabs(const char* name);
      void AddCustomNucleusParticles();
      TGCheckButton *eGCBAnimation;
      TGNumberEntry *fNumberEntry;
      TGNumberEntryField *fNumberEntryRun;
      TGCheckButton *fIsolate;
      TGHProgressBar* gProgress;
      TGHtml *fgHtml = 0;
      TGLabel *fExt_l;
      TGCheckButton *fExt_b;
      TGLabel *fLabelEve;
      TGRadioButton *fb_sel;
      std::string HandleEventTypeSelection(Int_t id);

      void GetRecoTracks();
      void GetMuTagInfo();
      void GetMuonSpectrometerInfo(TGeoShape *bigbox, TGeoMedium *air, TGeoShape *box);

      void AnalyzeScintVoxelsAndLayerOccupancy(bool drawPlots=false,bool clampOutOfRange=false, int expectedNz=20);

      TObjArray *get_selected(int printsel = 0);
      TObjArray *selected = new TObjArray;
      void print_selected();
      void OnSelectionChanged();

      TGCheckButton *fPrimary;
      TGCheckButton *fMuons;
      TGCheckButton *fProtons;
      TGCheckButton *fPions;
      TGCheckButton *fKaons;
      TGCheckButton *fSecondary;
      TGCheckButton *fEMShowers;
      TGCheckButton *fHadronShowers;
      TGCheckButton *fPixelTracker;
      TGCheckButton *fPixelRecoTrack;
      TGCheckButton *fRear;
      TGCheckButton *fMuTag;
      TGCheckButton *fGhostHits;
      TGCheckButton *fRecoVoxHits;
     

      TGCheckButton *fShortLivedParticle;

      TGCheckButton *fPandora;
      TGCheckButton *fZoom;
      TGComboBox* fEventTypeComboBox;
      TGHSlider *fHSlider_ext;
      TGNumberEntry *fRawDigitCut;
      Double_t timePitch = 0;
      Double_t extPitch = 0;
      Double_t extension = 0;
      TGNumberEntryField *fExt_NE;
      TGNumberEntry *fTSteps;
      
      void DumpMCTruth();
      int CharmParentID();
      bool CharmedEvent();
      void ShowRearCalos();

      TEveElement * GetSelectedElement () { TEveElement::List_i it= gEve->GetSelection()->EndChildren(); return *(--it);}
      //////
      bool IsCharmed(int pdg) {
        int abs_pdg = std::abs(pdg);	
        // Check if it matches charmed meson or baryon
        return (abs_pdg / 100 == 4) ||       // Charmed mesons: 4xx
        (abs_pdg / 1000 == 4) ||      // Charmed baryons: 4xxx
        (abs_pdg / 1000 == 104) ||    // Excited mesons: 104xx
        (abs_pdg / 10000 == 104);     // Excited baryons: 104xxx
      }
      ///////
    bool IsShortLivedParticle(int pdg_id, int& particleType);
    void TauDecayMode();
    void IdentifyTauDecayMode(const std::vector<std::pair<int, int>>& daughterTracks);
       bool ShortLivedParticleEvent();

std::vector<int> fSLPParentIDs;
    std::vector<std::string> fSLPNames;
    std::vector<int> fSLPTypes; // 1=charm, 2=tau, 3=other

      void GetPryMuon();
      void CharmDecayMode();
      Int_t fCharmCharge;
      Double_t fCharmEnergy;
      Int_t fTauCharge;
      Double_t fTauEnergy;

      Int_t fEventNumber = 0;
      Int_t fRunNumber = 0;
      Int_t fMaskNumber = 0;
      Int_t fCharmParentID = -1;
      Int_t fCharmDaughterID = -1;
      Int_t fTauParentID = -1;
      Int_t fTauDaughterID = -1;


      Int_t fPryMuonID = -1;
      Int_t fCurrentEventNumber;
      Int_t fnumLayers = 0;
      Int_t fnumChargedDaughters=-1;
      Int_t fnumNeutralDaughters=-1;
      Int_t fnumTauChargedDaughters=-1;
      Int_t fnumTauNeutralDaughters=-1;

      Int_t fCharm;
      Int_t fTau;
      std::string fdecayMode;
      std::string fcharmname=" ";
      std::string ftauDecayMode;  
      std::string ftauname=" ";

      void IdentifyCharmDecayMode(const std::vector<std::pair<int, int>> &daughterTracks);
      Double_t fdecay_vx, fdecay_vy, fdecay_vz;
      Double_t ftau_vx, ftau_vy, ftau_vz;

      Int_t fdecay_module;
      Double_t fdecayFlightLength;
      Int_t ftau_decay_module;
      Double_t ftauDecayFlightLength;

      std::vector<double> fmuTag_p;
      std::vector<double> fmuTag_px;
      std::vector<double> fmuTag_py;
      std::vector<double> fmuTag_pz;
      std::vector<double> fmuTag_E;
      std::vector<int> fmuTag_id;
      std::vector<double> fmuTag_dist;
      std::vector<int> fmuTag_alltrk;
      std::vector<int> fmuTag_muonSign;
      int fmuPry_id;


      Int_t fmuTag_muon;

      TcalEvent* fTcalEvent;
      TPOEvent *POevent;
      TPORecoEvent* fPORecoEvent;

    // FASERCal full container
    TGeoNode* fFaserCalNode = nullptr;
    // Rear ECAL has 24 modules → store them in a vector
    std::vector<TGeoNode*> fRearCalNodes;
    // Rear HCAL container
    TGeoNode* fRearHCalNode = nullptr;  
    // Muon Spectrometer container
    TGeoNode* fMuonSpectNode = nullptr;
 
    private:
      Bool_t fVisibleVol[10000] = { kFALSE };
      //std::map<TEveGeoShape*, TcalEvent->getfTracks()[0]> shapeTrackMap; // Auto-detect the type
      std::map<TEveGeoShape*, DigitizedTrack*> shapeTrackMap; 
  //Fastjet
  //std::vector<fastjet::PseudoJet> fJets;
      Bool_t ApplyIsolation = kFALSE;
      Bool_t ShowTruth = kFALSE;
      Bool_t ShowCRT = kFALSE;
      Bool_t ShowPMT = kFALSE;
      Bool_t ShowTPC = kFALSE;
      Bool_t ShowPandora = kFALSE;
      
      Bool_t fVerbose = kFALSE;
      
      //FaserCalData faserCalData;
      int ievent;

      TText *runText = nullptr;
      TText *eventypeText = nullptr;
      TText *energyText = nullptr;
      
      
      TEveSelection* fSelMgr = nullptr; // EVE selection manager




      TGeoVolume *primary;
      TGeoVolume *muonhit;
      TGeoVolume *pionhit;
      TGeoVolume *protonhit;
      TGeoVolume *kaonhit;
      TGeoVolume *secondary;
      TGeoVolume *secondary_em;
      TGeoVolume *secondary_had;
      TGeoVolume *si_tracker;
      TGeoVolume *ShortLivedP;
      TGeoVolume *rearecal;
      TGeoVolume *rearhcal;
      TGeoVolume *rearmucal;
      TGeoVolume *voxVol;
      TGeoVolume *muonSpectVol;

      TEveElementList* fDetectorElements;
      
      TEveElementList* fHitElements;
      TEveElementList* fZoomHitElements;
      TEveElementList* fPixelHitElements;
      TEveElementList* fPrimaryElements;
      TEveElementList* fSecondaryShowerElements;
      TEveElementList* fSecondaryHadShowerElements;
      TEveElementList* fShortLivedParticleHitElements;
      TEveElementList* fPixelRecoTrackElements;
      TEveElementList* fMuonHitElements;
      TEveElementList* fProtonHitElements;
      TEveElementList* fPionHitElements;
      TEveElementList* fKaonHitElements;

      TEveElementList* fRearECALElements;
      TEveElementList* fRearHCALElements;
      TEveElementList* fRearMuCALElements;

      TEveElementList* fMuTagHitElements;
      TEveElementList* fMuSpectHitElements;
      TEveElementList* fMuonSpectFitElements;

      TEveElementList* fClusterHitElements;
      TEveElementList* fVoxHitElements;
      TEveElementList* fVoxGhostElements;

      TEvePointSet* fTruthVertex = nullptr; 

      TCanvas* fEventInfoCanvas;
      ClassDef(FaserCalDisplay, 1)
    };
  
  
 
    
    
    class FaserCalDisplayCamera
    {
        RQ_OBJECT("FaserCalDisplayCamera")
    public:
        FaserCalDisplayCamera();
        virtual ~FaserCalDisplayCamera();
        enum { X, Y, Z, B };

        void SetCamera(int projection = Y);
        void Animation();
        void SetAnimationMode();
        void StartAnimation();
        void StopAnimation();
        void SaveAnimation(char *filename = NULL, int n = 100, int interval = 10);
        void SavePictures();
        void Snapshot(char *filename = NULL);
        void SetProjection();

    private:
        TGLViewer *eViewer = nullptr;
        TTimer *eTimer = nullptr;
        int eAnimationMode = 0;
        float eTheta = 3.14 / 2;
        float eAnimationAmplitude = 0.002;

        ClassDef(FaserCalDisplayCamera, 0)
    };

}  // namespace display

#endif  // __FaserCal_display_FaserCalDisplay__
