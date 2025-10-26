#include "FaserCalDisplay.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <TEveWindow.h>
#include <TEveViewer.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TGFrame.h>
#include <TGFont.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TEveText.h>
#include <TEveTrans.h> 
#include <TGDMLParse.h>
#include <TString.h>

#include <TGComboBox.h>

#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include <filesystem>
#include <vector>
#include <TLegend.h>

#include <TSystem.h>
#include <unordered_map>

namespace display
{
  FaserCalDisplay::FaserCalDisplay() {
    fNumberEntry = nullptr;
    fStatusBar = nullptr;
    fDetectorElements = new TEveElementList("Detector Elements");
    fHitElements = new TEveElementList("Hit Elements");
    fPrimaryElements = new TEveElementList("Primary Elements");  
    fSecondaryShowerElements = new TEveElementList("Secondary Shower Elements");  
    fSecondaryHadShowerElements = new TEveElementList("Secondary Hadron Shower Elements");  
    fPixelHitElements = new TEveElementList("Pixel Hit Elements");
    fClusterHitElements = new TEveElementList("Cluster Hit Elements");
    fShortLivedParticleHitElements =  new TEveElementList("Short Lived Particle Elements");
    fMuonHitElements = new TEveElementList("Muon Hit Elements");
    fProtonHitElements = new TEveElementList("Proton Hit Elements");
    fPionHitElements = new TEveElementList("Pion Hit Elements");
    fKaonHitElements = new TEveElementList("Kaon Hit Elements");
    fRearECALElements = new TEveElementList("RearECAL Elements");  
    fRearHCALElements = new TEveElementList("RearHCAL Elements");
    fRearMuCALElements = new TEveElementList("RearMuCAL Elements");
    fMuTagHitElements = new TEveElementList("MuTagHit Elements");
    fVoxHitElements = new TEveElementList("RecoVoxel Elements");
    fVoxGhostElements = new TEveElementList("GhostVoxel Elements");
    fMuSpectHitElements = new TEveElementList("MuonSpectHit Elements");
    fMuonSpectFitElements = new TEveElementList("MuonSpectFit Elements");
  }
  //////////////////////////////////////////////////////////
  FaserCalDisplay::~FaserCalDisplay() {}
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::GetDetector()
  {
    std::cout << "Starting GetDetector()" << std::endl;
    // Load the GDML file using TGeoManager::Import
    //TGeoManager::Import("/data/sw/FASERCAL/FASER/GeomGDML/geometry.gdml");
    // use for Run120 and v5.0
    //TGeoManager::Import("/home/hyperk/sw/FASERCAL/FASER_March2025/GeomGDML/geometry_v5.gdml");
    TGeoManager::Import("../../GeomGDML/geometry.gdml");
    //
    if (!gGeoManager) {
      std::cerr << "Failed to import GDML file." << std::endl;
      return;
    }
    std::cout << "GDML file imported successfully." << std::endl;
    //
    TGeoVolume* gdmlTop = gGeoManager->GetTopVolume();
    if (!gdmlTop) {
      std::cerr << "Failed to get top volume." << std::endl;
      return;
    }
    //std::cout << "Top volume set: " << gdmlTop->GetName() << std::endl;
    TGeoIterator nextNode(gdmlTop);
    TGeoNode* curNode;
    //
    while ((curNode = nextNode())) {
      TGeoVolume* vol = curNode->GetVolume();
      if (!vol) {
	std::cerr << "Volume is null for node: " << curNode->GetName() << std::endl;
	continue;
      }
      //
      TGeoShape* shape = vol->GetShape();
      if (!shape) {
	std::cerr << "Shape is null for volume: " << vol->GetName() << std::endl;
	continue;
      }
      //
      TString nodeName(curNode->GetName());
      TString nodePath;
      nextNode.GetPath(nodePath);
      //std::cout << "Processing node: " << nodeName << ", path: " << nodePath << std::endl;
      // Get the transformation matrix of the current node
      const TGeoMatrix* matrix = nextNode.GetCurrentMatrix();
      if (!matrix) {
	std::cerr << "Matrix is null for node: " << curNode->GetName() << std::endl;
	continue;
      }
      //
      const Double_t* trans = matrix->GetTranslation();
      const Double_t* rotMatrix = matrix->GetRotationMatrix();
      if (!trans || !rotMatrix) {
	std::cerr << "Transformation matrix is null for node: " << curNode->GetName() << std::endl;
	continue;
      }
      //
      TGeoRotation rotation;
      rotation.SetMatrix(rotMatrix);
      TGeoCombiTrans transform(trans[0], trans[1], trans[2], &rotation);
      // Create the TEveGeoShape for visualization
      TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
      eveShape->SetShape(shape);
      eveShape->SetMainTransparency(90); // Set transparency
      eveShape->SetTransMatrix(transform);
      //gEve->AddGlobalElement(eveShape);
      fDetectorElements->AddElement(eveShape);
    }
    gEve->AddGlobalElement(fDetectorElements);
    gEve->Redraw3D(kTRUE);
    std::cout << "GetDetector() completed." << std::endl;
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::GetEventDisplay()
  {
    gSystem->Setenv("ROOT_OPENGL_FALLBACK", "1");
    gROOT->SetBatch(kFALSE);
    
    // Force software rendering
    if (gSystem->Getenv("DISPLAY") == nullptr) {
        gSystem->Setenv("DISPLAY", ":0.0");
    }
    
    TEveManager* gEve = TEveManager::Create(kTRUE, "V");
    if (gROOT->IsBatch())
      gROOT->SetBatch(kFALSE);
    
    TEveRGBAPalette* pal = new TEveRGBAPalette(0, 1000);
    TEveViewer* ev = gEve->GetDefaultViewer();
    //ev->GetPickEvent()->Connect(this, "OnPick(TGLPhysicalShape*, TObject*, Int_t)");
    
    TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    TEveWindowPack* pack = slot->MakePack();
    pack->SetElementName("Multi View");
    pack->SetHorizontal();
    pack->SetShowTitleBar(kFALSE);
    
    pack->NewSlot()->MakeCurrent();
    TEveViewer* T3DView = gEve->SpawnNewViewer("Y-Z View", "");
    T3DView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
    T3DView->AddScene(gEve->GetGlobalScene());
    T3DView->AddScene(gEve->GetEventScene());
    
    pack = pack->NewSlot()->MakePack();
    pack->SetShowTitleBar(kFALSE);
    pack->NewSlot()->MakeCurrent();
    TEveViewer* TRPhiView = gEve->SpawnNewViewer("X-Y View", "");
    TRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    TRPhiView->AddScene(gEve->GetGlobalScene());
    TRPhiView->AddScene(gEve->GetEventScene());

    pack->NewSlot()->MakeCurrent();
    TEveViewer* TRhoZView = gEve->SpawnNewViewer("X-Z View", "");
    TRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
    TRhoZView->AddScene(gEve->GetGlobalScene());
    TRhoZView->AddScene(gEve->GetEventScene());

    TEveBrowser* browser = gEve->GetBrowser();
    browser->StartEmbedding(TRootBrowser::kLeft);

    TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
    frmMain->SetWindowName("XX GUI");
    frmMain->SetCleanup(kDeepCleanup);

    TGVerticalFrame* hf = new TGVerticalFrame(frmMain);
    TGGroupFrame* fGroupFrame2 = new TGGroupFrame(hf, "Event Display");
    fGroupFrame2->SetLayoutBroken(kTRUE);
    /////////////////////////////////////////////////////////////////////
    // Set Event Number Button
    int posy = 20;
        fNumberEntryRun = new TGNumberEntryField(fGroupFrame2, -1, 200025, TGNumberFormat::kNESInteger);
    fGroupFrame2->AddFrame(fNumberEntryRun, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 62, 2));
    fNumberEntryRun->MoveResize(20, posy, 90, 18);
    
    fNumberEntry = new TGNumberEntry(fGroupFrame2, 0, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 1000000);
    fGroupFrame2->AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 62, 2));
    fNumberEntry->MoveResize(120, posy, 90, 18);
    //MaskEntry = new TGNumberEntry(fGroupFrame2, 0, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 10000);
    //fGroupFrame2->AddFrame(fNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 62, 2));
    //fNumberEntry->MoveResize(20, posy, 90, 18);
    posy += 28;
    
    fEventTypeComboBox = new TGComboBox(fGroupFrame2, "Event Type");
    fEventTypeComboBox->AddEntry("NoMask", 0);
    fEventTypeComboBox->AddEntry("nueCC", 1);
    fEventTypeComboBox->AddEntry("numuCC", 2);
    fEventTypeComboBox->AddEntry("nutauCC", 3);
    fEventTypeComboBox->AddEntry("nuNC", 4);
    fEventTypeComboBox->AddEntry("nuES", 5);
    
    fEventTypeComboBox->Resize(150, 20);
    fGroupFrame2->AddFrame(fEventTypeComboBox, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fEventTypeComboBox->MoveResize(20, posy, 90, 20);

    TGTextButton* b = new TGTextButton(fGroupFrame2, "Set Event#");
    b->Connect("Clicked()", "display::FaserCalDisplay", this, "SetEventNumber()");
    fGroupFrame2->AddFrame(b, new TGLayoutHints(kLHintsExpandX));
    b->MoveResize(120, posy, 90, 18);
    posy += 28;
    // Show Event Button
    TGTextButton* fb = new TGTextButton(fGroupFrame2, "ShowEvent");
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowEvent()");
    fb->MoveResize(10, posy, 100, 18);
    fb->SetToolTipText("Show event");
    // Dump Event Button
    fb = new TGTextButton(fGroupFrame2, "DumpEvent");
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "DumpEvent()");
    fb->MoveResize(110, posy, 100, 18);
    fb->SetToolTipText("Dump event");

    posy += 25;
    // Next Event Button
    fb = new TGTextButton(fGroupFrame2, "Next");
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "NextEvent()");
    fb->MoveResize(10, posy, 100, 18);
    fb->SetToolTipText("Show next event");
    // Previous Event Button
    fb = new TGTextButton(fGroupFrame2, "Previous");
    fb->MoveResize(110, posy, 100, 18);
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "PreviousEvent()");
    fb->SetToolTipText("Show previous event");
    // Hide Detector Geometry
    posy += 25;
    fIsolate = new TGCheckButton(fGroupFrame2, "Hide Geometry");
    fIsolate->MoveResize(10, posy, 0, 18);
    fIsolate->Connect("Clicked()", "display::FaserCalDisplay", this, "Isolate()");

    fRear = new TGCheckButton(fGroupFrame2, "RearE/H/MuCal");
    fRear->MoveResize(115, posy, 0, 18);
    fRear->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowRearCalos()");

    ////////////
    posy += 25;
    fPrimary = new TGCheckButton(fGroupFrame2, "1ryParticles");
    fPrimary->MoveResize(10, posy, 0, 18);
    fPrimary->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowPrimary()");

    fShortLivedParticle = new TGCheckButton(fGroupFrame2, "SLP");
    fShortLivedParticle->MoveResize(110, posy, 0, 18);
    fShortLivedParticle->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowShortLivedParticle()");

    fMuons = new TGCheckButton(fGroupFrame2, "Muons");
    fMuons->MoveResize(160, posy, 0, 18);
    fMuons->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowMuonParticles()");

    posy += 25;
    fKaons = new TGCheckButton(fGroupFrame2, "Kaons");
    fKaons->MoveResize(10, posy, 0, 18);
    fKaons->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowKaonParticles()");
    fPions = new TGCheckButton(fGroupFrame2, "Pions");
    fPions->MoveResize(80, posy, 0, 18);
    fPions->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowPionParticles()");
    fProtons = new TGCheckButton(fGroupFrame2, "Protons");
    fProtons->MoveResize(150, posy, 0, 18);
    fProtons->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowProtonParticles()");  

    posy += 25;
    fEMShowers = new TGCheckButton(fGroupFrame2, "EM Showers");
    fEMShowers->MoveResize(10, posy, 0, 18);
    fEMShowers->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowSecondaryShowers()");

    fHadronShowers = new TGCheckButton(fGroupFrame2, "Hadron Showers");
    fHadronShowers->MoveResize(110, posy, 0, 18);
    fHadronShowers->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowSecondaryHadShowers()");

    posy += 25;
    fPixelTracker = new TGCheckButton(fGroupFrame2, "Pixel Tracker");
    fPixelTracker->MoveResize(10, posy, 0, 18);
    fPixelTracker->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowPixelHits()");

    fPixelRecoTrack = new TGCheckButton(fGroupFrame2, "RecoTracks");
    fPixelRecoTrack->MoveResize(110, posy, 0, 18);
    fPixelRecoTrack->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowPixelRecoTrack()");

    posy += 25;

    fRecoVoxHits = new TGCheckButton(fGroupFrame2, "RecoVox");
    fRecoVoxHits->MoveResize(10, posy, 0, 18);
    fRecoVoxHits->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowRecoVoxelHits()");

    fGhostHits = new TGCheckButton(fGroupFrame2, "Ghost");
    fGhostHits->MoveResize(110, posy, 0, 18);
    fGhostHits->Connect("Clicked()", "display::FaserCalDisplay", this, "ShowOnlyGhostHits()");

    posy += 15;


    // Connect the selection event to the callback
    //gEve->GetSelection()->Connect("Selected(TEveElement*)", 
    //    "display::FaserCalDisplay", this, "OnShapeSelected(TEveElement*)");
    // std::cout << "Selection callback connected!" << std::endl;
    // here is the selection
    //gEve->GetSelection()->Connect("Selected(TEveElement*)", 
    //"display::FaserCalDisplay", this, "OnShapeSelected(TEveElement*)");

    hf->AddFrame(fGroupFrame2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame2->MoveResize(0, 200, 260, 270);


    TGGroupFrame* fGroupFrame30 = new TGGroupFrame(hf, "Analysis");
    fGroupFrame30->SetLayoutBroken(kTRUE);
    fb_sel = new TGRadioButton(fGroupFrame30,"Selection");
    fb_sel->MoveResize(10,10,0,18);
    fb_sel->Connect("Clicked()","display::FaserCalDisplay",this,"EnablePicking()");
    fb_sel->SetState((EButtonState)1);
    hf->AddFrame(fGroupFrame30, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame30->MoveResize(0, 220, 250, 50);


    TGGroupFrame* fGroupFrame3 = new TGGroupFrame(hf, "Plots");
    fGroupFrame3->SetLayoutBroken(kTRUE);
    fb = new TGTextButton(fGroupFrame3, "HitsInCubes");
    fb->Connect("Clicked()","display::FaserCalDisplay",this,"CountHitsInCube()");
    fb->MoveResize(10, 20, 110, 18);
    
    fb = new TGTextButton(fGroupFrame3, "AllEvents");
    fb->MoveResize(10, 40, 110, 18);
    fb->Connect("Clicked()","display::FaserCalDisplay",this,"PlotHistogramsFromRootFile()");

      
    fb->SetToolTipText("CRT plots");
    fb = new TGTextButton(fGroupFrame3, "EnergyDepositions");
    fb->MoveResize(125, 20, 110, 18);
    fb = new TGTextButton(fGroupFrame3, "JetReconstructions");
    fb->Connect("Clicked()","display::FaserCalDisplay",this,"JetReconstructions()");
    fb->MoveResize(125, 40, 110, 18);

    fb = new TGTextButton(fGroupFrame3, "SummaryTable");
    fb->MoveResize(10, 60, 110, 18);

    hf->AddFrame(fGroupFrame3, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame3->MoveResize(0, 220, 240, 100);


    TGGroupFrame* fGroupFrame4 = new TGGroupFrame(hf, "Zooming");
    fGroupFrame4->SetLayoutBroken(kTRUE);
    posy = 0;
    int posx = 10;
    int dx = 20;

    //fExt_b = new TGCheckButton(fGroupFrame4, "");
    //fGroupFrame4->AddFrame(fExt_b, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    //fExt_b->MoveResize(10, posy, 110, 18);

    //TGLabel* lTStep = new TGLabel(fGroupFrame4, "TSteps");
    //lTStep->MoveResize(30, posy, 0, 18);

    //fTSteps = new TGNumberEntry(fGroupFrame4, 50, 6, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, 0, 10000);
    //fTSteps->MoveResize(80, posy, 0, 18);

    //fExt_l = new TGLabel(fGroupFrame4, " ns ");
    //fExt_l->MoveResize(140, posy, 20, 18);
    //fExt_NE = new TGNumberEntryField(fGroupFrame4, -1, 0, TGNumberFormat::kNESInteger);
    //fExt_NE->MoveResize(155, posy, 75, 18);

    posy += 20;
    fHSlider_ext = new TGHSlider(fGroupFrame4, 220, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
    fHSlider_ext->SetRange(0, 50);
    fHSlider_ext->SetPosition(0);
    fHSlider_ext->SetPosition(extPitch);
    fHSlider_ext->MoveResize(5, posy, 225, 15);
    fHSlider_ext->Connect("PositionChanged(Int_t)", "display::FaserCalDisplay", this, "DoSlider(Int_t)");

    hf->AddFrame(fGroupFrame4, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame4->MoveResize(0, 200, 240, 50);



    
    TGGroupFrame* fGroupFrame5 = new TGGroupFrame(hf, "WaveForms");
    fGroupFrame5->SetLayoutBroken(kTRUE);
    posx = 10;
    dx = 20;
    posy = 20;

    hf->AddFrame(fGroupFrame5, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame5->MoveResize(0, 210, 240, 50);

    TGGroupFrame* fGroupFrame6 = new TGGroupFrame(hf, "Saving");
    fGroupFrame6->SetLayoutBroken(kTRUE);
    posx = 10;
    dx = 20;
    posy = 20;

    eGCBAnimation = new TGCheckButton(fGroupFrame6, "Animation");
    eGCBAnimation->MoveResize(10, posy, 80, 20);
  
    fb = new TGTextButton(fGroupFrame6, "Save animation");
    fb->MoveResize(125, posy, 100, 18);

    posy += 21;
    fb = new TGTextButton(fGroupFrame6, "Snapshot");
    fb->MoveResize(10, posy, 100, 18);
    fb = new TGTextButton(fGroupFrame6, "Save projections");
    fb->MoveResize(125, posy, 100, 18);
    // Change Background Color
    posy += 21;
    fb = new TGTextButton(fGroupFrame6, "ChangeBackgroundColour");
    fb->MoveResize(25, posy, 150, 20);
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "BackgroundColor()");
    // Exit Button
    posy += 21;
    fb = new TGTextButton(fGroupFrame6, "Exit");
    fb->MoveResize(55, posy, 100, 20);
    fb->Connect("Clicked()", "display::FaserCalDisplay", this, "DoExit()");

    hf->AddFrame(fGroupFrame6, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));
    fGroupFrame6->MoveResize(0, 260, 240, 240);

    frmMain->AddFrame(hf);

    // status bar
    Int_t parts[] = {45,20,20,15};
    fStatusBar = gEve->GetBrowser()->GetStatusBar();
    fStatusBar->SetParts(parts,4);
    fStatusBar->SetText("Welcome to the FASERCal EventDisplay",0);
    
    frmMain->MapSubwindows();
    frmMain->Resize();
    frmMain->MapWindow();

    slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    fgHtml = new TGHtml(0, 100, 100);
    TEveWindowFrame* wf = slot->MakeFrame(fgHtml);
    fgHtml->MapSubwindows();
    wf->SetElementName("Summary");

    browser->StopEmbedding();
    browser->SetTabTitle("Main", 0);
    gEve->GetBrowser()->GetTabRight()->SetTab(1);

    gEve->Redraw3D(kTRUE);

    GetDetector();
    ShowAxis();
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::ShowAxis()
  {
    TString axisName = "Coordinate system";
    Float_t axisArrowLength = 10.;
    Float_t axisArrowTubeR = 0.01;
    Color_t axisColor = kMagenta;
    Int_t axisLabelFontSize = 15;
    Int_t axisLabelOffset = 5;

    TEveElementList* axis = new TEveElementList(axisName, axisName);
    TEveArrow* xAxis = new TEveArrow(axisArrowLength, 0., 0.);
    xAxis->SetMainColor(axisColor);
    xAxis->SetTubeR(axisArrowTubeR);
    xAxis->SetConeR(axisArrowTubeR);
    xAxis->SetElementNameTitle("X Axis", "X");
    axis->AddElement(xAxis);
    TEveText* xAxisLabel = new TEveText("X");
    xAxisLabel->SetMainColor(axisColor);
    xAxisLabel->SetFontSize(axisLabelFontSize);
    xAxisLabel->SetLighting(kTRUE);
    Double_t trans[3];
    trans[0] = axisArrowLength + axisLabelOffset;
    trans[1] = 0;
    trans[2] = 0;
    xAxisLabel->PtrMainTrans()->SetPos(trans);
    axis->AddElement(xAxisLabel);

    TEveArrow* yAxis = new TEveArrow(0., axisArrowLength, 0.);
    yAxis->SetMainColor(axisColor);
    yAxis->SetTubeR(axisArrowTubeR);
    yAxis->SetConeR(axisArrowTubeR);
    yAxis->SetElementNameTitle("Y Axis", "Y");
    axis->AddElement(yAxis);
    TEveText* yAxisLabel = new TEveText("Y");
    yAxisLabel->SetMainColor(axisColor);
    yAxisLabel->SetFontSize(axisLabelFontSize);
    yAxisLabel->SetLighting(kTRUE);
    trans[0] = 0;
    trans[1] = axisArrowLength + axisLabelOffset;
    trans[2] = 0;
    yAxisLabel->PtrMainTrans()->SetPos(trans);
    axis->AddElement(yAxisLabel);

    TEveArrow* zAxis = new TEveArrow(0., 0., axisArrowLength);
    zAxis->SetMainColor(axisColor);
    zAxis->SetTubeR(axisArrowTubeR);
    zAxis->SetConeR(axisArrowTubeR);
    zAxis->SetElementNameTitle("Z Axis", "Z");
    axis->AddElement(zAxis);
    TEveText* zAxisLabel = new TEveText("Z");
    zAxisLabel->SetMainColor(axisColor);
    zAxisLabel->SetFontSize(axisLabelFontSize);
    zAxisLabel->SetLighting(kTRUE);
    trans[0] = 0;
    trans[1] = 0;
    trans[2] = axisArrowLength + axisLabelOffset;
    zAxisLabel->PtrMainTrans()->SetPos(trans);
    axis->AddElement(zAxisLabel);

    gEve->AddGlobalElement(axis);
    gEve->Redraw3D(kTRUE);

 
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::BackgroundColor()
  {
    gStyle->SetPalette(-1);
    gEve->GetDefaultGLViewer()->SetClearColor(0);
    gEve->FullRedraw3D(kTRUE);
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::DoExit()
  {
    std::cout << "Exit application..." << std::endl;
    gROOT->Reset();
    gApplication->Terminate(0);
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::SetEventNumber()
  {
     if (!fNumberEntry) {
        std::cerr << "fNumberEntry is not initialized!" << std::endl;
        return;
    }
    fEventNumber = fNumberEntry->GetNumber();
    if (!fNumberEntryRun) {
        std::cerr << "fNumberEntryRun is not initialized!" << std::endl;
        return;
    }
    fRunNumber = fNumberEntryRun->GetNumber();
    fMaskNumber = fEventTypeComboBox->GetSelected();

    if (!fStatusBar) {
        std::cerr << "fStatusBar is not initialized!" << std::endl;
        return;
    }

    fStatusBar->GetBarPart(0)->SetBackgroundColor(0xffffff);
    fStatusBar->SetText(Form("   Showing Run#: %i Event#: %i type: %i  ", fRunNumber, fEventNumber, fMaskNumber), 0);
    std::cout << "RunNumber is "<< fRunNumber <<  " Event Number is " << fEventNumber << " " << fMaskNumber << std::endl;
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::Isolate()
  {
    if (!(fIsolate->IsOn())) {
      ApplyIsolation = kFALSE;
      //gEve->GetGlobalScene()->SetRnrSelf(kTRUE);
      fDetectorElements->SetRnrState(kTRUE);
      fHitElements->SetRnrSelf(kTRUE);  // Ensure hits are shown
      //fRearECALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown                               
      //fRearHCALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown
      //fRearMuCALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown 
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
      std::cout << "Detector geometry will be shown" << std::endl;
    } else {
      ApplyIsolation = kTRUE;
      //gEve->GetGlobalScene()->SetRnrSelf(kFALSE);
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrSelf(kTRUE);  // Ensure hits are shown
      //fRearECALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown
      //fRearHCALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown 
      //fRearMuCALElements->SetRnrSelf(kTRUE);  // Ensure hits are shown 
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
      std::cout << "Detector geometry will be isolated" << std::endl;
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowRearCalos()
  { 
    if (fRear->IsOn())
    {
      std::cout << "Rear" << std::endl;
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kFALSE);
      fSecondaryShowerElements->SetRnrState(kFALSE);
      fSecondaryHadShowerElements->SetRnrState(kFALSE);
      fPixelHitElements->SetRnrState(kFALSE);  
      fShortLivedParticleHitElements->SetRnrState(kFALSE);  
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fPrimaryElements->SetRnrState(kFALSE);  
      fPionHitElements->SetRnrState(kFALSE);
      fKaonHitElements->SetRnrState(kFALSE);
      fProtonHitElements->SetRnrState(kFALSE);
      fMuonHitElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);
      fMuSpectHitElements->SetRnrState(kFALSE);
      if(!fIsolate->IsOn())
	      fDetectorElements->SetRnrState(kTRUE);
      else
	      fDetectorElements->SetRnrState(kFALSE);
      
      gEve->AddGlobalElement(fRearECALElements);
      gEve->AddGlobalElement(fRearHCALElements);
      gEve->AddGlobalElement(fRearMuCALElements);
      fRearECALElements->SetRnrState(kTRUE);
      fRearHCALElements->SetRnrState(kTRUE);
      fRearMuCALElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    } else {
      fHitElements->SetRnrState(kTRUE);
      fRearECALElements->SetRnrState(kFALSE);
      fRearHCALElements->SetRnrState(kFALSE);
      fRearMuCALElements->SetRnrState(kFALSE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    

    }
  }

  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowPrimary()
  { 
    if (fPrimary->IsOn())
    {
      std::cout << "Primary" << std::endl;
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kFALSE);
      fSecondaryShowerElements->SetRnrState(kFALSE);
      fSecondaryHadShowerElements->SetRnrState(kFALSE);
      fPixelHitElements->SetRnrState(kFALSE);  
      fShortLivedParticleHitElements->SetRnrState(kFALSE);  
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fRearECALElements->SetRnrState(kFALSE);
      fRearHCALElements->SetRnrState(kFALSE);
      fRearMuCALElements->SetRnrState(kFALSE);
      fPionHitElements->SetRnrState(kFALSE);
      fKaonHitElements->SetRnrState(kFALSE);
      fProtonHitElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);
      fMuSpectHitElements->SetRnrState(kFALSE);

      if(!fIsolate->IsOn())
	      fDetectorElements->SetRnrState(kTRUE);
      else
	      fDetectorElements->SetRnrState(kFALSE);
      
      gEve->AddGlobalElement(fPrimaryElements);
      fPrimaryElements->SetRnrState(kTRUE);  
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowMuonParticles()
  { 
    if (fMuons->IsOn())
    {
	    std::cout << "Muons" << std::endl;
	    fDetectorElements->SetRnrState(kFALSE);
	    fHitElements->SetRnrState(kFALSE);
      fSecondaryShowerElements->SetRnrState(kFALSE);
	    fSecondaryHadShowerElements->SetRnrState(kFALSE);
	    fPixelHitElements->SetRnrState(kFALSE);  
	    fShortLivedParticleHitElements->SetRnrState(kFALSE);  
      fPrimaryElements->SetRnrState(kFALSE);  
      fPixelRecoTrackElements->SetRnrState(kFALSE);
	    fRearECALElements->SetRnrState(kFALSE);
	    fRearHCALElements->SetRnrState(kFALSE);
	    fRearMuCALElements->SetRnrState(kFALSE);
	    fPionHitElements->SetRnrState(kFALSE);
      fKaonHitElements->SetRnrState(kFALSE);
      fProtonHitElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);

      if(!fIsolate->IsOn())
	      fDetectorElements->SetRnrState(kTRUE);
	    else
	      fDetectorElements->SetRnrState(kFALSE);
	
	    gEve->AddGlobalElement(fMuonHitElements);
	    fMuonHitElements->SetRnrState(kTRUE);  
      gEve->AddGlobalElement(fMuSpectHitElements);
      fMuSpectHitElements->SetRnrState(kTRUE);
      gEve->AddGlobalElement(fMuonSpectFitElements);
      fMuonSpectFitElements->SetRnrState(kTRUE);
	    gStyle->SetPalette(-1);
	    gEve->FullRedraw3D(kTRUE);
    } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }

//////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowPionParticles()
  { 
    if (fPions->IsOn())
      {
	std::cout << "Pions" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);  
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fPrimaryElements->SetRnrState(kFALSE);  
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
  fMuonHitElements->SetRnrState(kFALSE);
  fProtonHitElements->SetRnrState(kFALSE);
  fKaonHitElements->SetRnrState(kFALSE);
  fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);
	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	gEve->AddGlobalElement(fPionHitElements);
	fPionHitElements->SetRnrState(kTRUE);  
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
      } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  
//////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowKaonParticles()
  { 
    if (fKaons->IsOn())
      {
	std::cout << "Kaons" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);  
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fPrimaryElements->SetRnrState(kFALSE);  
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
  fMuonHitElements->SetRnrState(kFALSE);
  fProtonHitElements->SetRnrState(kFALSE);
  fPionHitElements->SetRnrState(kFALSE);
  fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);
	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	gEve->AddGlobalElement(fKaonHitElements);
	fKaonHitElements->SetRnrState(kTRUE);  
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
      } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  
//////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowProtonParticles()
  { 
    if (fProtons->IsOn())
      {
	std::cout << "Protons" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);  
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fPrimaryElements->SetRnrState(kFALSE);  
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
  fMuonHitElements->SetRnrState(kFALSE);
  fPionHitElements->SetRnrState(kFALSE);
  fKaonHitElements->SetRnrState(kFALSE);
  fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);
	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	gEve->AddGlobalElement(fProtonHitElements);
	fProtonHitElements->SetRnrState(kTRUE);  
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
      } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  

  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowShortLivedParticle()
  { 
    if (fShortLivedParticle->IsOn())
    {
      std::cout << "Short Lived Particle" << std::endl;
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kFALSE);
      fSecondaryShowerElements->SetRnrState(kFALSE);
      fSecondaryHadShowerElements->SetRnrState(kFALSE);
      fPixelHitElements->SetRnrState(kFALSE);  
      fPrimaryElements->SetRnrState(kFALSE);  
      fMuonHitElements->SetRnrState(kFALSE);  
      fPionHitElements->SetRnrState(kFALSE);
      fKaonHitElements->SetRnrState(kFALSE);
      fProtonHitElements->SetRnrState(kFALSE);
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fRearECALElements->SetRnrState(kFALSE);
      fRearHCALElements->SetRnrState(kFALSE);
      fRearMuCALElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);
      fMuSpectHitElements->SetRnrState(kFALSE);

      if(!fIsolate->IsOn())
	fDetectorElements->SetRnrState(kTRUE);
      else
	fDetectorElements->SetRnrState(kFALSE);
      
      gEve->AddGlobalElement(fShortLivedParticleHitElements);
      fShortLivedParticleHitElements->SetRnrState(kTRUE);  
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowSecondaryShowers()
  {
    if (fEMShowers->IsOn())
    {
      std::cout << "EM shower" << std::endl;
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kFALSE);
      fPrimaryElements->SetRnrState(kFALSE);  
      fSecondaryHadShowerElements->SetRnrState(kFALSE);
      fPixelHitElements->SetRnrState(kFALSE);  
      fShortLivedParticleHitElements->SetRnrState(kFALSE);  
      fMuonHitElements->SetRnrState(kFALSE);  
      fPionHitElements->SetRnrState(kFALSE);
      fKaonHitElements->SetRnrState(kFALSE);
      fProtonHitElements->SetRnrState(kFALSE);
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fRearECALElements->SetRnrState(kFALSE);
      fRearHCALElements->SetRnrState(kFALSE);
      fRearMuCALElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);
      fMuSpectHitElements->SetRnrState(kFALSE);

      if(!fIsolate->IsOn())
	fDetectorElements->SetRnrState(kTRUE);
      else
	fDetectorElements->SetRnrState(kFALSE);

      gEve->AddGlobalElement(fSecondaryShowerElements);
      fSecondaryShowerElements->SetRnrState(kTRUE);  
      gStyle->SetPalette(-1);
      gEve->Redraw3D(kTRUE);
    } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowSecondaryHadShowers()
  {    
    if (fHadronShowers->IsOn())
    {
      std::cout << "Had shower" << std::endl;
      fDetectorElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kFALSE);
      fSecondaryShowerElements->SetRnrState(kFALSE);
      fPrimaryElements->SetRnrState(kFALSE);
      fPixelHitElements->SetRnrState(kFALSE);  
      fShortLivedParticleHitElements->SetRnrState(kFALSE);  
      fMuonHitElements->SetRnrState(kFALSE);  
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fRearECALElements->SetRnrState(kFALSE);
      fRearHCALElements->SetRnrState(kFALSE);
      fRearMuCALElements->SetRnrState(kFALSE);
      fMuTagHitElements->SetRnrState(kFALSE);
      fMuSpectHitElements->SetRnrState(kFALSE);
      
      if(!fIsolate->IsOn())
	      fDetectorElements->SetRnrState(kTRUE);
      else
	      fDetectorElements->SetRnrState(kFALSE);
 	
      gEve->AddGlobalElement(fSecondaryHadShowerElements);
      fSecondaryHadShowerElements->SetRnrState(kTRUE);  
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);      
    } else {
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowPixelHits()
  {    
    if (fPixelTracker->IsOn())
      {
	std::cout << "Hits on PixelTracker" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fPrimaryElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);  
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fMuonHitElements->SetRnrState(kFALSE);  
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
	fMuTagHitElements->SetRnrState(kFALSE);
	fMuSpectHitElements->SetRnrState(kFALSE);

	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);
 	
	gEve->AddGlobalElement(fPixelHitElements);
	fPixelHitElements->SetRnrState(kTRUE);  
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);      
      } else {
      fPixelHitElements->SetRnrState(kFALSE);  
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowPixelRecoTrack()
  {
    TEveLine *trackLine;
    if(fPixelRecoTrack->IsOn())
      {
	std::cout << "Show Reconstructed tracks" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fPrimaryElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fMuonHitElements->SetRnrState(kFALSE);  
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
  fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);
	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	gEve->AddGlobalElement(fPixelRecoTrackElements);
	fPixelRecoTrackElements->SetRnrState(kTRUE);
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
	EnablePicking();
      } else {
      fPixelRecoTrackElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowRecoVoxelHits()
  {
    TEveLine *trackLine;
    if(fRecoVoxHits->IsOn())
      {
	std::cout << "Show Reconstructed tracks" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fPrimaryElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fMuonHitElements->SetRnrState(kFALSE);  
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);

	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	fVoxHitElements->SetRnrState(kTRUE);
	gEve->AddGlobalElement(fVoxHitElements);
	fVoxHitElements->SetRnrState(kTRUE);
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
	EnablePicking();
      } else {
      fVoxHitElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    }
  }  
   //////////////////////////////////////////////////////////
  void FaserCalDisplay::ShowOnlyGhostHits()
  {
    TEveLine *trackLine;
    if(fGhostHits->IsOn())
      {
	std::cout << "Show only ghost" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fPrimaryElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);
	fPixelHitElements->SetRnrState(kFALSE);
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fMuonHitElements->SetRnrState(kFALSE);  
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fVoxHitElements->SetRnrState(kFALSE);
	fMuTagHitElements->SetRnrState(kFALSE);
  fMuSpectHitElements->SetRnrState(kFALSE);

	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);

	fVoxGhostElements->SetRnrState(kTRUE);
	gEve->AddGlobalElement(fVoxGhostElements);
	fVoxGhostElements->SetRnrState(kTRUE);
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);
	EnablePicking();
      } else {
      fVoxGhostElements->SetRnrState(kFALSE);
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);
    }
  }   


  //////////////////////////////////////////////////////////
  /*
  void FaserCalDisplay::ShowClusterHits()
  {    
    if (fClusterHits->IsOn())
      {
	std::cout << "Cluster Hits" << std::endl;
	fDetectorElements->SetRnrState(kFALSE);
	fHitElements->SetRnrState(kFALSE);
	fSecondaryShowerElements->SetRnrState(kFALSE);
	fPrimaryElements->SetRnrState(kFALSE);
	fSecondaryHadShowerElements->SetRnrState(kFALSE);  
	fShortLivedParticleHitElements->SetRnrState(kFALSE);  
	fMuonHitElements->SetRnrState(kFALSE);  
	fPixelRecoTrackElements->SetRnrState(kFALSE);
	fRearECALElements->SetRnrState(kFALSE);
	fRearHCALElements->SetRnrState(kFALSE);
	fRearMuCALElements->SetRnrState(kFALSE);

	if(!fIsolate->IsOn())
	  fDetectorElements->SetRnrState(kTRUE);
	else
	  fDetectorElements->SetRnrState(kFALSE);
 	
	gEve->AddGlobalElement(fClusterHitElements);
	fClusterHitElements->SetRnrState(kTRUE);  
	gStyle->SetPalette(-1);
	gEve->FullRedraw3D(kTRUE);      
      } else {
      fClusterHitElements->SetRnrState(kFALSE);  
      fHitElements->SetRnrState(kTRUE);
      gStyle->SetPalette(-1);
      gEve->FullRedraw3D(kTRUE);    
    }
  }
  */
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::EnablePicking()
  {
    gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_Element);
    get_selected(1);
    
  }
  ////////////////
  TObjArray *FaserCalDisplay::get_selected(int printSel)
  {
    selected->Clear(); // Clear any previously selected items
    
    if (gEve->GetCurrentEvent() != NULL)
      {
        TEveSelection* sel = gEve->GetSelection();
        for (TEveElement::List_i it = sel->BeginChildren(); it != sel->EndChildren(); ++it)
	  {
            TEveElement* e = *it;
            if (!e) continue; // Skip null elements
	    
            // Print the element's name
            if (printSel)
	      printf("Selected Element: %s\n", e->GetElementName());
	    
            // Access the user data if available
            TObject* obj = static_cast<TObject*>(e->GetUserData());
            if (obj)
	      {
		std::cout << "User Data: " << obj->ClassName() << std::endl;
                selected->Add(obj); // Add to the selected list
	        
		// Check if the object is a TGeoVolume and look up the track information
                TGeoVolume* volume = dynamic_cast<TGeoVolume*>(obj);
                if (volume)
		  {
                    auto it = volumeTrackMap.find(volume);
                    if (it != volumeTrackMap.end())
		      {
                        DigitizedTrack* track = it->second;
                        std::cout << "Track ID: " << track->fhitIDs[0] << std::endl; // Example of accessing track information
		      }
		  }
		
	      }
	  }
      }
    
    return selected; // Return the list of selected objects 
  }
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::OnPick(TGLPhysicalShape *shape, TObject *obj, Int_t ) 
  {
    if (!obj) return;
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::OnShapeSelected(TEveElement* selectedElement) 
  {
    // Cast the selected element to TEveGeoShape
    TEveGeoShape* selectedShape = dynamic_cast<TEveGeoShape*>(selectedElement);
    if (!selectedShape) {
      std::cerr << "No valid shape selected!" << std::endl;
      return;
    }
    std::cout << "Selected shape pointer: " << selectedShape << std::endl;
    
    std::cout << "Here i am" << std::endl;
    // Highlight the selected shape
    selectedShape->SetMainColor(kYellow); // Change color to yellow
    gEve->Redraw3D(kTRUE); // Redraw the display for immediate feedback
    
    // Retrieve associated metadata from the map
    auto it = shapeTrackMap.find(selectedShape);
    if (it != shapeTrackMap.end()) {
      DigitizedTrack* curtrack = it->second;
      // Print track information
      std::cout << "Selected Track Information:" << std::endl;
      std::cout << "  PDG Code: " << curtrack->fPDG << std::endl;
      std::cout << "  Energy Deposits: " << curtrack->fEnergyDeposits[0] << std::endl;
      std::cout << "  First Hit Position: " << curtrack->fhitIDs[0] << std::endl;
      std::cout << "  Number of Hits: " << curtrack->fhitIDs.size() << std::endl;
    } else {
      std::cerr << "No track associated with the selected shape!" << std::endl;
    }
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::DoSlider(Int_t position)
  {
    ROOT::Math::XYZVector vertex(POevent->prim_vx.x(), POevent->prim_vx.y(), POevent->prim_vx.z());
    
    double zoomRadius = static_cast<double>(position) * fTcalEvent->geom_detector.fSandwichLength;
    
    // Function to process TEveElementList recursively
    auto processElement = [&](TEveElement* element, auto& processElementRef) -> void {
      TEveGeoShape* shape = dynamic_cast<TEveGeoShape*>(element);
      if (shape) {
	//std::cout << "Shape Name: " << shape->GetName() << ", Type: TEveGeoShape" << std::endl;
	TEveTrans& trans = shape->RefMainTrans();  // Get a reference to the transformation matrix
	TVector3 pos = trans.GetPos();  // Get the position as a TVector3
	double x = pos.X();  // Access the X component
	double y = pos.Y();  // Access the Y component
	double z = pos.Z();  // Access the Z component
	// Calculate the distance from the vertex
	double distance = std::sqrt(std::pow(x - vertex.X(), 2) +
				    std::pow(y - vertex.Y(), 2) +
				    std::pow(z - vertex.Z(), 2));
	
	// Show or hide the element based on the zoom radius
	if (distance <= zoomRadius) {
	  shape->SetRnrSelf(kTRUE);  // Show the shape
	} else {
	  shape->SetRnrSelf(kFALSE);  // Hide the shape
	}
      } else {
	TEveElementList* list = dynamic_cast<TEveElementList*>(element);
	if (list) {
	  //std::cout << "Processing TEveElementList: " << list->GetName() << std::endl;
	  for (TEveElement::List_i child = list->BeginChildren(); child != list->EndChildren(); ++child) {
	    processElementRef(*child, processElementRef);  // Dereference the iterator to pass TEveElement*
	  }
	} else {
	  std::cout << "Unknown element type: " << typeid(*element).name() << std::endl;
	}
      }
    };
    // Iterate over the elements and process them
    for (TEveElement::List_i element = fHitElements->BeginChildren(); element != fHitElements->EndChildren(); ++element)
      {
        processElement(*element, processElement); 
      }
    
    gEve->FullRedraw3D(kTRUE);
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::CleanViewer()
  {
    gEve->GetViewers()->DeleteAnnotations();
    gEve->GetCurrentEvent()->DestroyElements();
    TEveElement* top = gEve->GetCurrentEvent();
    gEve->GetGlobalScene()->SetRnrSelf(kFALSE);
    TEveRGBAPalette* pal = new TEveRGBAPalette(0, 50);
    gGeoManager->GetTopVolume()->SetVisRaytrace(true);
    gEve->Redraw3D(kTRUE);
  }
  //////////////////////////////////////////////////////////  
  TCanvas* FaserCalDisplay::CreateCanvas(const char* plot_name, int tb)
  {
    TCanvas* myCan = (TCanvas*)gROOT->FindObject(plot_name);
    if (myCan) {
      myCan->Clear();
    } else {
      gEve->GetBrowser()->StartEmbedding(tb);
      gROOT->ProcessLineFast("new TCanvas");
      myCan = (TCanvas*)gPad;
      myCan->SetName(plot_name);
      gEve->GetBrowser()->StopEmbedding(plot_name);
    }
    return myCan;
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::CleanCanvas()
  {
    TCanvas* myCan;
    myCan = (TCanvas*)gROOT->FindObject("EventInfo");
    if (myCan) myCan->Clear();
    myCan = (TCanvas*)gROOT->FindObject("AllEvents");
    if (myCan) myCan->Clear();
    myCan = (TCanvas*)gROOT->FindObject("HitsInCubes");
    if (myCan) myCan->Clear();
    myCan = (TCanvas*)gROOT->FindObject("Clusters");
    if (myCan) myCan->Clear();
  }
  //////////////////////////////////////////////////////////  
  TCanvas* FaserCalDisplay::CreateTabs(const char* name)
  {
    TCanvas* cx = (TCanvas*)gROOT->FindObject(name);
    if (cx) return cx;
    
    gEve->GetBrowser()->StartEmbedding(1, -1);
    TGMainFrame* fMainTabFrame = new TGMainFrame(gClient->GetRoot(), 10, 10, kMainFrame | kVerticalFrame);
    fMainTabFrame->SetName(name);
    TGTab* fMainTab = new TGTab(fMainTabFrame, 300, 300);
    fMainTab->SetTab(1);
    fMainTab->Resize();
    fMainTabFrame->AddFrame(fMainTab, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    
    fMainTabFrame->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputModeless);
    fMainTabFrame->MapSubwindows();
    fMainTabFrame->Resize();
    fMainTabFrame->MapWindow();
    gEve->GetBrowser()->StopEmbedding(name);
    
    TGCompositeFrame* subTabCompositeFrame = fMainTab->AddTab(name);
    subTabCompositeFrame->SetLayoutManager(new TGVerticalLayout(subTabCompositeFrame));
    TRootEmbeddedCanvas* rootEmbeddedCanvas = new TRootEmbeddedCanvas(0, subTabCompositeFrame, 300, 300);
    Int_t wfRootEmbeddedCanvas736a = rootEmbeddedCanvas->GetCanvasWindowId();
    cx = new TCanvas(name, 10, 10, wfRootEmbeddedCanvas736a);
    
    rootEmbeddedCanvas->AdoptCanvas(cx);
    subTabCompositeFrame->AddFrame(rootEmbeddedCanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    
    return cx;
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::SetMyStyle()
  {
    Int_t myFont = 42;
    Double_t myWidth = 2;
    Double_t myTSize = 0.05;
    
    //gROOT->SetStyle("Plain");
    TStyle* myStyle = new TStyle("myStyle", "FaserCalDisplay plots style");
    myStyle->SetPalette(1,0); 
    myStyle->SetOptStat(1110);
    myStyle->SetOptTitle(1);
    myStyle->SetOptDate(0);
    myStyle->SetStatColor(10);
    myStyle->SetStatFontSize(0.05);
    myStyle->SetStatH(0.26);
    myStyle->SetStatW(0.26);
    myStyle->SetTitleFont(42,"xyz"); // font option 
    myStyle->SetLabelFont(42,"xyz");
    myStyle->SetLabelSize(0.045,"xyz"); // size of axis value font
    myStyle->SetTitleSize(0.053,"xz"); // size of axis title font
    myStyle->SetTitleSize(0.053,"y"); // size of axis title font
    myStyle->SetTitleOffset(0.9,"x");
    myStyle->SetTitleOffset(1.3,"y");
    myStyle->SetTitleOffset(1.0,"z");
    myStyle->SetNdivisions(10, "x");
    myStyle->SetNdivisions(10, "y");
    myStyle->SetPadBottomMargin(0.13); //margins...
    myStyle->SetPadTopMargin(0.12);
    myStyle->SetPadLeftMargin(0.16);
    myStyle->SetPadRightMargin(0.16);
    myStyle->SetTitleFillColor(10);
    myStyle->SetLineWidth(2);
    // default canvas options
    myStyle->SetCanvasDefW(700);
    myStyle->SetCanvasDefH(600);
    //myStyle->SetCanvasColor(10);
    myStyle->SetCanvasColor(0);// canvas...
    myStyle->SetCanvasBorderMode(0);
    //myStyle->SetCanvasBorderMode(-1);
    myStyle->SetCanvasBorderSize(0);
    //myStyle->SetCanvasBorderSize(1);
    myStyle->SetPadColor(0);
    myStyle->SetPadBorderSize(1);
    myStyle->SetPadBorderMode(-1);
    myStyle->SetPadGridX(0); // grids, tickmarks
    myStyle->SetPadGridY(0);
    //myStyle->SetPadTickX(1);
    //myStyle->SetPadTickY(1);
    myStyle->SetFrameBorderSize(1);
    myStyle->SetFrameBorderMode(-1);
    //myStyle->SetFrameBorderMode(0);
    myStyle->SetFrameFillColor(0);
    //myStyle->SetFrameFillColor(10);
    myStyle->SetFrameLineWidth(2);
    myStyle->SetHistLineWidth(2.0);
    myStyle->SetPaperSize(20,24); // US letter size
    gROOT->SetStyle("myStyle");
  }
  ////////////////////////////////////////////////////////// 
  void FaserCalDisplayCamera::SetProjection()
  {
    TGButton* btn = (TGButton*)gTQSender;
    int id = btn->WidgetId();
    FaserCalDisplayCamera::SetCamera(id);
    std::cout << "ummmmmmmmmmmmm " << id << std::endl;
  }
  ////////////////////////////////////////////////////////// 
  void FaserCalDisplayCamera::SetCamera(int projection)
  {
    TGTabElement* tab = gEve->GetBrowser()->GetTabRight()->GetCurrentTab();
    TEveViewerList* viewers = gEve->GetViewers();
    TGLViewer* eview = nullptr;
    for (TEveElement::List_i it = viewers->BeginChildren(); it != viewers->EndChildren(); it++) {
      if (strcmp(((TEveViewer*)(*it))->GetName(), tab->GetName()) == 0) {
	eview = ((TEveViewer*)(*it))->GetGLViewer();
	break;
      }
    }
    if (eview == nullptr) {
      eview = gEve->GetDefaultGLViewer();
    }
    switch (projection) {
    case X: {
      std::cout << "XPROJECTION" << std::endl;
      gEve->GetDefaultGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
      TGLOrthoCamera& cam = (TGLOrthoCamera&)(eview->CurrentCamera());
      cam.Reset();
      cam.RotateRad(-3.14159 * 0.5, 0.0);
      cam.SetEnableRotate(1);
    } break;
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::AddCustomNucleusParticles()
  {
    struct CustomParticle {
      int pdgCode;
      const char* name;
      const char* title;
      double mass;
    };
    
    CustomParticle particles[] = {
        {1000050110, "B11", "Boron-11", 11.0093054},
        {1000060120, "C12", "Carbon-12", 12.0},
        {1000060130, "C13", "Carbon-13", 13.0033548378},
        {1000080160, "O16", "Oxygen-16", 15.99491461956},
        {1000080180, "O18", "Oxygen-18", 17.9991610},
        {1000110230, "Na23", "Sodium-23", 22.98976928},
        {1000120260, "Mg26", "Magnesium-26", 25.98259297},
        {1000130270, "Al27", "Aluminum-27", 26.9815385},
        {1000140280, "Si28", "Silicon-28", 27.97692653465},
        {1000140290, "Si29", "Silicon-29", 28.9764946653},
        {1000140300, "Si30", "Silicon-30", 29.973770136},
        {1000190390, "K39", "Potassium-39", 38.9637064864},
        {1000200400, "Ca40", "Calcium-40", 39.96259098},
        {1000200420, "Ca42", "Calcium-42", 41.95861801},
        {1000200430, "Ca43", "Calcium-43", 42.9587666},
        {1000200440, "Ca44", "Calcium-44", 43.9554818},
        {1000260540, "Fe54", "Iron-54", 53.9396105},
        {1000260560, "Fe56", "Iron-56", 55.9349375},
        {1000260570, "Fe57", "Iron-57", 56.9353940},
        {1000350810, "Br81", "Bromine-81", 80.9162906},
        {1000471070, "Ag107", "Silver-107", 106.9050916},
        {1000601420, "Nd142", "Neodymium-142", 141.9077233},
        {1000601430, "Nd143", "Neodymium-143", 142.9098143},
        {1000601440, "Nd144", "Neodymium-144", 143.9100873},
        {1000601450, "Nd145", "Neodymium-145", 144.9125736},
        {1000601460, "Nd146", "Neodymium-146", 145.9131169},
        {1000601480, "Nd148", "Neodymium-148", 147.9168930},
        {1000601500, "Nd150", "Neodymium-150", 149.920891},
        {1000741820, "W182", "Tungsten-182", 181.9482042},
        {1000741830, "W183", "Tungsten-183", 182.9502230},
        {1000741840, "W184", "Tungsten-184", 183.9509312},
        {1000741860, "W186", "Tungsten-186", 185.9543641},
        {1000822040, "Pb204", "Lead-204", 203.9730440},
        {1000822060, "Pb206", "Lead-206", 205.9744653},
        {1000822070, "Pb207", "Lead-207", 206.9758969},
        {1000822080, "Pb208", "Lead-208", 207.9766521}
    };

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();

    for (const auto& particle : particles) {
      if (pdgDB->GetParticle(particle.pdgCode) == nullptr) {	
        pdgDB->AddParticle(particle.name, particle.title, particle.mass, kFALSE, 0.0, 0, "Nucleus", particle.pdgCode);
      }
    }
    std::cout << "Custom nucleus particles added to TDatabasePDG." << std::endl;
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::DumpEvent()
  {
    CleanCanvas();
    std::cout << "Event to be dumped: " << fEventNumber <<  " from run number " << fRunNumber << " and eventtype " << fMaskNumber << std::endl;
    //FaserCalData faserCalData;
    //faserCalData.LoadEvent(fEventNumber);

    fTcalEvent = new TcalEvent();
    POevent = new TPOEvent();
    fTcalEvent->Load_event("input/", fRunNumber, fEventNumber, fMaskNumber, POevent);
    //LoadEvent(fRunNumber,fEventNumber, fMaskNumber);
    TCanvas *myCan = CreateCanvas("EventInfo",1);
    gPad->Update();
    myCan->cd();
    TPaveText* infoText = new TPaveText(0.1, 0.1, 0.9, 0.9);
    infoText->AddText(Form("Run: %d Event: %d", POevent->run_number,fEventNumber));
    std::ostringstream eventtype;
    int pdgin = POevent->in_neutrino.m_pdg_id;
    switch(pdgin) {
      case -12:
      case 12:
        eventtype << "nu_e";
        break;
      case -14:
      case  14:
        eventtype << "nu_mu";
        break;
      case -16:
      case  16:
        eventtype << "nu_tau";
        break;
      default:
        eventtype << " ?? ";
    }
    if(POevent->isCC) {
      eventtype << " CC ";
    } else {
      eventtype << " NC ";
    }

    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    AddCustomNucleusParticles();
    int chargedMultiplicity = 0;
    int neutralMultiplicity = 0;
    int totalMultiplicity = 0;
    std::string target;
    for (size_t i=0; i<POevent->n_particles(); i++) {
      struct PO& aPO = POevent->POs[i];
      TParticlePDG* particle = pdgDB->GetParticle(aPO.m_pdg_id);
      int charge = particle ? particle->Charge() : 0;
      if(aPO.m_status==1)
	{
	  totalMultiplicity++;
	  if (charge != 0) {
            chargedMultiplicity++;
	  } else {
            neutralMultiplicity++;
	  }
	}
        if(aPO.m_status==4)
	{
	  if(!POevent->is_neutrino(aPO.m_pdg_id))
	    {
	      target = particle ? particle->GetName() : "Unknown";
	      std::cout << "check target: "<< target << std::endl;
	    }
	}
    }
    infoText->AddText(Form("Interaction: %s + %s ", eventtype.str().c_str(),target.c_str()));
    infoText->AddText(Form("Neutrino Energy: %.3f GeV",POevent->in_neutrino.m_energy));
    infoText->AddText(Form("Primary Vtx: %.3f %.3f %.3f",POevent->prim_vx.x(),POevent->prim_vx.y(),POevent->prim_vx.z()));
    infoText->AddText(Form("Multiplicity: %d [Charged: %d Neutral: %d]",totalMultiplicity,chargedMultiplicity,neutralMultiplicity));
    if(POevent->isCC)
      {
	struct PO& outlepton = POevent->out_lepton;
	TParticlePDG* particle = pdgDB->GetParticle(outlepton.m_pdg_id);
	double plepton = TMath::Sqrt(outlepton.m_px*outlepton.m_px+outlepton.m_py*outlepton.m_py+outlepton.m_pz*outlepton.m_pz);
	infoText->AddText(Form("Lepton: %s P = %.1f GeV/c2 E= %.1f GeV",particle->GetName(),plepton,outlepton.m_energy));
      }
    infoText->AddText(Form("Number of segments: %lu", fTcalEvent->getfTracks().size()));
    infoText->AddText(Form("----------------------------"));
    infoText->AddText(Form("Jet: %.3f %.3f %.3f", POevent->jetpx, POevent->jetpy, POevent->jetpz));
    infoText->AddText(Form("Sum final state: %.3f %.3f %.3f", POevent->spx, POevent->spy, POevent->spz));
    infoText->AddText(Form("Sum final state (Vis): %.3f %.3f %.3f", POevent->vis_spx, POevent->vis_spy, POevent->vis_spz));
    infoText->AddText(Form("Ptmiss: %.3f  Evis: %.3f ", POevent->ptmiss, POevent->Evis));
    infoText->Draw();
    gPad->Update();
    myCan->Update();
    //DumpMCTruth();
  }
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
bool FaserCalDisplay::IsShortLivedParticle(int pdg_id, int& particleType)
    {
        particleType = 0; // 0 = not SLP, 1 = charm, 2 = tau
        // Charm particles
        if (abs(pdg_id) == 411 || abs(pdg_id) == 421 || abs(pdg_id) == 431 ||
            abs(pdg_id) == 4122 || abs(pdg_id) == 4132 || abs(pdg_id) == 4232 ||
            abs(pdg_id) == 4332) {
            particleType = 1; // Charm
            std::cout << "Charm found for this event" << std::endl;
            return true;
        }
        // Tau lepton
        if (abs(pdg_id) == 15) {
            particleType = 2; // Tau
            std::cout << "Tau found for this event" << std::endl;

            return true;
        }
        return false;
    }
    //////////////////////////////////////////////////////////  
bool FaserCalDisplay::ShortLivedParticleEvent()
{
    fCharmParentID = -1;
    fTauParentID = -1;
    fSLPParentIDs.clear();
    fSLPNames.clear();
    fSLPTypes.clear();
    
    // Reset all variables
    fcharmname = "";
    ftauname = "";
    fCharmCharge = 0;
    fTauCharge = 0;
    fCharmEnergy = 0.0;
    fTauEnergy = 0.0;
    fnumChargedDaughters = 0;
    fnumNeutralDaughters = 0;
    fnumTauChargedDaughters = 0;
    fnumTauNeutralDaughters = 0;
    fdecayMode = "";
    ftauDecayMode = "";
    
    int PrimaryTrackID = 0;
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    bool foundSLP = false;
    
    for (size_t i = 0; i < POevent->n_particles(); i++) 
    {
        struct PO& aPO = POevent->POs[i];
        if (aPO.m_status == 1)
        {
            PrimaryTrackID++;
            
            if (abs(aPO.m_pdg_id) == 13) {
                fPryMuonID = PrimaryTrackID;  
                GetPryMuon();
            }
            
            int particleType;
            if (IsShortLivedParticle(aPO.m_pdg_id, particleType)) 
            {
              std::cout << "Check for shortlived" << std::endl;
                foundSLP = true;
                fSLPParentIDs.push_back(PrimaryTrackID);
                fSLPTypes.push_back(particleType);
                
                TParticlePDG* particle = pdgDB->GetParticle(aPO.m_pdg_id);
                std::string particleName = particle ? particle->GetName() : "Unknown";
                fSLPNames.push_back(particleName);
                
                if (particleType == 1) { // Charm
                    fCharmParentID = PrimaryTrackID;
                    fCharmCharge = particle ? particle->Charge() : 0;
                    fcharmname = particleName;
                    fCharmEnergy = aPO.m_energy;
                    std::cout << "Charm found: trackID " << fCharmParentID << std::endl;

                    CharmDecayMode();
                }
                else if (particleType == 2) { // Tau
                    fTauParentID = PrimaryTrackID;
                    fTauCharge = particle ? particle->Charge() : 0;
                    ftauname = particleName;
                    fTauEnergy = aPO.m_energy;
                    std::cout << "Tau found: trackID " << fTauParentID << std::endl;
                    TauDecayMode();
                }
            }
        }
    }   
    return foundSLP;
}
/////////////////////////////////////////////////////////////

  bool FaserCalDisplay::CharmedEvent()
  {
    //fTcalEvent = new TcalEvent();
    //POevent = new TPOEvent();
    fCharmParentID = -1;
    fCharmDaughterID = -1;
    fPryMuonID = -1;
    fcharmname = "";
    fCharmCharge = 0;
    fCharmEnergy = 0.0;
    fnumChargedDaughters = 0;
    fnumNeutralDaughters = 0;
    fdecayMode = "";
    int PrimaryTrackID = 0;
     TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    for (size_t i=0; i<POevent->n_particles(); i++) 
      {
        struct PO& aPO = POevent->POs[i];
        if(aPO.m_status==1)
        {
          PrimaryTrackID++;
	  
	  if(abs(aPO.m_pdg_id)==13)
	    {
	      std::cout << "check Pry Muon" << std::endl;
	      fPryMuonID = PrimaryTrackID;  
	      GetPryMuon();
	    }
          if(IsCharmed(aPO.m_pdg_id)) 
          {
            fCharmParentID = PrimaryTrackID;
            fCharmCharge = pdgDB->GetParticle(aPO.m_pdg_id)->Charge();
            fcharmname = pdgDB->GetParticle(aPO.m_pdg_id)->GetName();
	          fCharmEnergy = aPO.m_energy;
            CharmDecayMode();
            return true;
          }
        }
      }   
    return false;
  }
  //////////////////////////////////////////////////////////  
  int FaserCalDisplay::CharmParentID()
  {
    //fTcalEvent = new TcalEvent();
    //POevent = new TPOEvent();
    int CharmTrackID = -1;
    int PrimaryTrackID = 0;
    for (size_t i=0; i<POevent->n_particles(); i++) 
      {
	struct PO& aPO = POevent->POs[i];
	if(aPO.m_status==1)
	  PrimaryTrackID++;
	if(IsCharmed(aPO.m_pdg_id)) 
	  {
	    std::cout << "found a charm " << std::endl;
	    CharmTrackID =  PrimaryTrackID;
	  }
      }   
    return CharmTrackID;
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::GetPryMuon()
  {
    for (const auto& track : fTcalEvent->getfTracks())
      {
	fmuPry_id = -1;
	//////
	// here i want to extract trackId of the muon
	// trackId might keep changing at every scattering !
	// not sure if i will get the id at Mutag 
	if(fPryMuonID==track->fprimaryID && abs(track->fPDG)==13)
	  {
	    //fmuPry_id = track->ftrackID; // it is not the one 
	    std::cout << "Muon ---- " << track->fparentID << " "
			  << track->ftrackID << " "
			  << track->fprimaryID << " "
			  << track->fPDG << " "
		      << std::endl;
	  }
	break;
      }
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::CharmDecayMode()
  {
    std::vector<std::pair<int, int>> daughterTracks; // Stores (trackID, PDG)
    std::cout << " inside charm decay mode......" << std::endl;
    TVector3 primaryVertex(POevent->prim_vx.x(), POevent->prim_vx.y(), POevent->prim_vx.z());
    TVector3 decayVertex(0.0, 0.0, 0.0);

    for (const auto& track : fTcalEvent->getfTracks())
      {
        if(fCharmParentID==track->fparentID)
        for ( size_t i = 0; i < track->fhitIDs.size(); i++)
          {
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            long module = fTcalEvent->getChannelModulefromID(track->fhitIDs[i]);

            std::cout << "CharmDaughter " << track->fparentID << " "
                << track->ftrackID << " "
                << track->fprimaryID << " "
                << track->fPDG << " "
                << track->fEnergyDeposits[i] << " "
                << position.x() << " "
                << position.y() << " "
                << position.z() << " "
                << module << " "
                << std::endl;
                daughterTracks.emplace_back(track->ftrackID, track->fPDG);
           //check if muon daughter exist and reaching to MuTag
           ROOT::Math::XYZVector positionRearHCal = fTcalEvent->getChannelXYZRearHCal(10);
           if(abs(track->fPDG)==13)
             {
	       fCharmDaughterID = track->ftrackID;
               std::cout << track->fparentID << " "
                         << track->ftrackID << " "
                         << track->fprimaryID << " "
                         << track->fPDG << " "
                         << track->fEnergyDeposits[i] << " "
                         << position.x() << " "
                         << position.y() << " "
                         << position.z() << " "
                         << module << " "
                         << std::endl;
               std::cout << "Dimuon ........" <<  std::endl;
               std::cout << positionRearHCal.X() << " " << positionRearHCal.Y() << " " << positionRearHCal.Z() << std::endl; 
             }
          /////
           if(i==0)
           {
             //std::cout << "Charm decay positions " 
             //        << position.x() << " "
             //        << position.y() << " "
             //        << position.z() << " "
             //        << module << " "
             //        << std::endl;
             fdecay_vx = position.x();
             fdecay_vy = position.y();
             fdecay_vz = position.z();
             fdecay_module = module;
             //
             decayVertex.SetXYZ(position.x(), position.y(), position.z());
             fdecayFlightLength = (decayVertex - primaryVertex).Mag();
             std::cout << "Charm Flight Length: " << fdecayFlightLength << " mm" << std::endl;
           }   
	      }   
      }
    IdentifyCharmDecayMode(daughterTracks);
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::IdentifyCharmDecayMode(const std::vector<std::pair<int, int>> &daughterTracks)
  {
    std::set<int> uniqueChargedTracks;
    std::set<int> uniqueNeutralTracks;
    int numPions = 0, numKaons = 0, numLeptons = 0, numMuons = 0, numElectrons = 0, numProtons = 0;
    bool hasPhoton = false, hasNeutrino = false;

    // Classify particles based on unique track IDs
    for (const auto &[trackID, pdg] : daughterTracks)
    {
        if (abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 2212)
        {
            uniqueChargedTracks.insert(trackID); // Store unique charged track IDs

            if (abs(pdg) == 211) numPions++;
            if (abs(pdg) == 321) numKaons++;
            if (abs(pdg) == 11) numElectrons++, numLeptons++;
            if (abs(pdg) == 13) numMuons++, numLeptons++;
            if (abs(pdg) == 2212) numProtons++;
        }
        else if (abs(pdg) == 111 || abs(pdg) == 130 || abs(pdg) == 310 || abs(pdg) == 22 || abs(pdg) == 2112 || abs(pdg) == 12 || abs(pdg) == 14)
        {
            uniqueNeutralTracks.insert(trackID); // Store unique neutral track IDs
            if (abs(pdg) == 22) hasPhoton = true;
            if (abs(pdg) == 12 || abs(pdg) == 14) hasNeutrino = true;
        }
    }
    fCharm = -1;
    // Count unique tracks
    fnumChargedDaughters = uniqueChargedTracks.size();
    fnumNeutralDaughters = uniqueNeutralTracks.size();
    if (numLeptons > 0) {
      if (numMuons > 0){
            fdecayMode = "Leptonic: Muonic";
	    fCharm = 1;
      }else if (numElectrons > 0){
            fdecayMode = "Leptonic: Electronic";
	    fCharm = 2;
      }else{
            fdecayMode = "Leptonic";
	    fCharm = 3;
      }
    }
    else if (numKaons > 0 && numPions > 0) {
        fdecayMode = "Hadronic: KaonPion";
	fCharm = 4;
    }
    else if (numPions > 1) {
        fdecayMode = "Hadronic: MultiPion";
	fCharm = 5;
    }
    else if (numKaons > 1) {
        fdecayMode = "Hadronic: Kaonic";
	fCharm = 6;
    }
    else {
        fdecayMode = "Unknown Decay Mode";
	fCharm = 0;
    }
    //std::cout << "Charm decay mode: " << fdecayMode << std::endl;
    //std::cout << "Number of charged daughters: " << fnumChargedDaughters << std::endl;
    // Add topology classification based on unique charged tracks
  }
//////////////////////////////////////////////////////////////
  void FaserCalDisplay::TauDecayMode()
{
    std::vector<std::pair<int, int>> daughterTracks; // Stores (trackID, PDG)
    std::cout << "Inside tau decay mode analysis..." << std::endl;
    
    TVector3 primaryVertex(POevent->prim_vx.x(), POevent->prim_vx.y(), POevent->prim_vx.z());
    TVector3 decayVertex(0.0, 0.0, 0.0);

    for (const auto& track : fTcalEvent->getfTracks())
    {
        if (fTauParentID == track->fparentID)
        {
            for (size_t i = 0; i < track->fhitIDs.size(); i++)
            {
                ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
                long module = fTcalEvent->getChannelModulefromID(track->fhitIDs[i]);

                std::cout << "TauDaughter " << track->fparentID << " "
                         << track->ftrackID << " "
                         << track->fprimaryID << " "
                         << track->fPDG << " "
                         << track->fEnergyDeposits[i] << " "
                         << position.x() << " "
                         << position.y() << " "
                         << position.z() << " "
                         << module << " "
                         << std::endl;
                
                daughterTracks.emplace_back(track->ftrackID, track->fPDG);
                
                // Check for muon daughter reaching MuTag
                if (abs(track->fPDG) == 13) {
                    fTauDaughterID = track->ftrackID;
                    std::cout << "Tau->Muon decay found, reaching detector..." << std::endl;
                }
                
                if (i == 0) {
                    ftau_vx = position.x();
                    ftau_vy = position.y();
                    ftau_vz = position.z();
                    ftau_decay_module = module;
                    
                    decayVertex.SetXYZ(position.x(), position.y(), position.z());
                    ftauDecayFlightLength = (decayVertex - primaryVertex).Mag();
                    std::cout << "Tau Flight Length: " << ftauDecayFlightLength << " mm" << std::endl;
                }
            }
        }
    }
    IdentifyTauDecayMode(daughterTracks);
}

void FaserCalDisplay::IdentifyTauDecayMode(const std::vector<std::pair<int, int>>& daughterTracks)
{
    std::set<int> uniqueChargedTracks;
    std::set<int> uniqueNeutralTracks;
    
    int numPions = 0, numKaons = 0, numLeptons = 0, numMuons = 0, numElectrons = 0;
    int numNeutrinos = 0, numPhotons = 0, numNeutralPions = 0;
    bool hasNeutralKaon = false;

    // Classify particles based on unique track IDs
    for (const auto& [trackID, pdg] : daughterTracks)
    {
        // Charged particles
        if (abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 11 || abs(pdg) == 13) {
            uniqueChargedTracks.insert(trackID);
            
            if (abs(pdg) == 211) numPions++;
            if (abs(pdg) == 321) numKaons++;
            if (abs(pdg) == 11) numElectrons++, numLeptons++;
            if (abs(pdg) == 13) numMuons++, numLeptons++;
        }
        // Neutral particles
        else if (abs(pdg) == 111 || abs(pdg) == 130 || abs(pdg) == 310 || 
                 abs(pdg) == 22 || abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) {
            uniqueNeutralTracks.insert(trackID);
            
            if (abs(pdg) == 111) numNeutralPions++;
            if (abs(pdg) == 130 || abs(pdg) == 310) hasNeutralKaon = true;
            if (abs(pdg) == 22) numPhotons++;
            if (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16) numNeutrinos++;
        }
    }
    
    fnumTauChargedDaughters = uniqueChargedTracks.size();
    fnumTauNeutralDaughters = uniqueNeutralTracks.size();
    fTau = -1;
    
    // Tau decay mode classification
    if (numMuons > 0) {
        ftauDecayMode = "Tau -> Muon";
        fTau = 1; // Muonic decay
    }
    else if (numElectrons > 0) {
        ftauDecayMode = "Tau -> Electron";
        fTau = 2; // Electronic decay
    }
    else if (numPions >1 && fnumTauChargedDaughters == 1) {
        ftauDecayMode = "Tau -> Pion (1-prong)";
        fTau = 3; // 1-prong hadronic
    }
    else if (fnumTauChargedDaughters == 1 && numNeutralPions > 0) {
        ftauDecayMode = "Tau -> Charged + Neutrals (1-prong)";
        fTau = 4; // 1-prong with neutrals
    }
    else if (fnumTauChargedDaughters == 3) {
        ftauDecayMode = "Tau -> 3-prong hadronic";
        fTau = 5; // 3-prong hadronic
    }
    else if (numKaons > 0) {
        ftauDecayMode = "Tau -> Kaonic decay";
        fTau = 6; // Kaonic decay
    }
    else {
        ftauDecayMode = "Tau -> Other/Unknown";
        fTau = 0; // Unknown
    }
    
    std::cout << "Tau decay mode: " << ftauDecayMode << std::endl;
    std::cout << "Number of charged daughters: " << fnumTauChargedDaughters << std::endl;
    std::cout << "Number of neutral daughters: " << fnumTauNeutralDaughters << std::endl;
    std::cout << "Number of Pions: " << numPions << std::endl;
    std::cout << "Number of Kaons: " << numKaons << std::endl;
    std::cout << "Number of Leptons: " << numLeptons << std::endl;
    std::cout << "Number of Neutrinos: " << numElectrons << std::endl;
    std::cout << "Number of Photons: " << numMuons << std::endl;
}

  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::DumpMCTruth()
  {
    std::cout << "Event to be dumped: " << fEventNumber <<  " from run number " << fRunNumber << " and eventtype " << fMaskNumber << std::endl;
    //LoadEvent(fRunNumber,fEventNumber, fMaskNumber);

    fTcalEvent = new TcalEvent();
    POevent = new TPOEvent();
    std::vector<struct PO> charmdecay;
    int CharmTrackID = 0;
    double Charm_px, Charm_py, Charm_pz;
    
    std::ofstream fileOut("log.txt");    
    // Redirecting cout to write to "log.txt"
    std::cout.rdbuf(fileOut.rdbuf());


    int PrimaryTrackID = 0;

    fTcalEvent->Load_event("input/", fRunNumber, fEventNumber, fMaskNumber, POevent);
    std::cout << POevent->run_number << " "
	      << fEventNumber << " " 
	      << POevent->in_neutrino.m_pdg_id << " "
	      << POevent->isCC << " "  
	      << POevent->n_particles() << " "  
	      << POevent->in_neutrino.m_energy << " "  
	      << POevent->prim_vx.x() << " " << POevent->prim_vx.y() << " " << POevent->prim_vx.z() << " "  
	      << std::endl; //<< " "  
    for (size_t i=0; i<POevent->n_particles(); i++) 
      {
	struct PO& aPO = POevent->POs[i];
	if(aPO.m_status==1)
	  PrimaryTrackID++;
	std::cout << i << " " << PrimaryTrackID << " " << aPO.m_pdg_id << " "
		  << aPO.m_track_id << " "
		  << aPO.m_px << " "
		  << aPO.m_py << " "
		  << aPO.m_pz << " "
		  << aPO.m_vx_decay << " "
		  << aPO.m_vy_decay << " "
		  << aPO.m_vz_decay << " "
		  << aPO.m_energy << " "
		  << aPO.nparent << " "
		  << aPO.m_status << std::endl; 
	////////
	if(IsCharmed(aPO.m_pdg_id)) 
	  {
	    std::cout << "found a charm " << std::endl;
	    CharmTrackID =  PrimaryTrackID;
	    Charm_px = aPO.m_px;
	    Charm_py = aPO.m_py;	    
	    Charm_pz = aPO.m_pz;
	  }
      }

    std::cout << "Number of hits" << fTcalEvent->getfTracks().size() << std::endl;
    for (const auto& track : fTcalEvent->getfTracks())
      {
        size_t nhits = track->fhitIDs.size();
        for ( size_t i = 0; i < nhits; i++)
          {
            long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
            if(hittype == 0 && track->fEnergyDeposits[i] < 0.5)continue;
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
	    long module = fTcalEvent->getChannelModulefromID(track->fhitIDs[i]);
	    if(CharmTrackID==track->fprimaryID)
	      {
		fmuPry_id = track->fparentID;
		std::cout << track->fparentID << " "
			  << track->ftrackID << " "
			  << track->fprimaryID << " "
			  << track->fPDG << " "
			  << position.x() << " "
			  << position.y() << " "
			  << position.z() << " "
			  << module << " "
			  << track->fEnergyDeposits[i] << " "
			  << std::endl; //<< " "
	      }
	  }
      }
  }
  //////////////////////////////////////////////////////////  
  void FaserCalDisplay::ShowEvent()
  {
    std::cout << "Event to be shown " << fEventNumber << std::endl;
    DrawEvent(fRunNumber, fEventNumber, fMaskNumber);
    //DumpEvent();
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::NextEvent()
  {
    fCurrentEventNumber = fEventNumber;
    fEventNumber = fCurrentEventNumber + 1;
    std::cout << "Next Event Number is " << fEventNumber << std::endl;
    DrawEvent(fRunNumber, fEventNumber, fMaskNumber);
    //DumpEvent();
    fStatusBar->GetBarPart(0)->SetBackgroundColor(0xffffff);
    fStatusBar->SetText(Form("   Showing Run#: %i Event#: %i type: %i  ", fRunNumber, fEventNumber, fMaskNumber), 0);
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::PreviousEvent()
  {
    fCurrentEventNumber = fEventNumber;
    fEventNumber = fCurrentEventNumber - 1;
    std::cout << "Previous Event Number is " << fEventNumber << std::endl;
    DrawEvent(fRunNumber, fEventNumber, fMaskNumber);
    //DumpEvent();
    fStatusBar->GetBarPart(0)->SetBackgroundColor(0xffffff);
    fStatusBar->SetText(Form("   Showing Run#: %i Event#: %i type: %i  ", fRunNumber, fEventNumber, fMaskNumber), 0);
  }
  //////////////////////////////////////////////////////////
  std::string FaserCalDisplay::HandleEventTypeSelection(Int_t id)
  {
    std::string selectedEventType;
    
    switch (id)
      {
      case 1:
	selectedEventType = "nueCC";
	break;
      case 2:
	selectedEventType = "numuCC";
	break;
      case 3:
	selectedEventType = "nutauCC";
	break;
      case 4:
	selectedEventType = "nuNC";
	break;
      case 5:
	selectedEventType = "nuES";
	break;
      default:
	selectedEventType = "NoMask";
	break;
      }
    
    std::cout << "Selected Event Type: " << selectedEventType << std::endl;
    return selectedEventType;
  }
  //////////////////////////////////////////////////////////
  //void FaserCalData::LoadEvent(int ievent)
  void FaserCalDisplay::LoadEvent(int irun, int ievent, int imask)
  {
    
    std::string base_path = "input/";

    std::string mask = HandleEventTypeSelection(imask);
    
    std::string file_path;
    if (imask!=0) 
      file_path = base_path + "FASERG4-Tcalevent_" +std::to_string(irun)+"_"+std::to_string(ievent)+"_"+mask+".root";
    else
      file_path = base_path + "FASERG4-Tcalevent_" +std::to_string(irun)+"_"+std::to_string(ievent)+".root";

    // Create an instance of TcalEvent and TPOEvent
    fTcalEvent = new TcalEvent();
    POevent = new TPOEvent();

    std::cout << file_path << std::endl;
    if (!gSystem->AccessPathName(file_path.c_str())) {  
      //
      fTcalEvent -> Load_event(base_path, irun, ievent, imask, POevent);
      std::cout << "Transverse size " << fTcalEvent->geom_detector.fScintillatorSizeX << " mm " << std::endl;
      std::cout << "Total size of one sandwich layer " << fTcalEvent->geom_detector.fSandwichLength << " mm " << std::endl;
      std::cout << "Number of layers " << fTcalEvent->geom_detector.NRep << std::endl;
      std::cout << "Voxel size " << fTcalEvent->geom_detector.fScintillatorVoxelSize << " mm " << std::endl;
      std::cout << " copied digitized tracks " << fTcalEvent->getfTracks().size() << std::endl;
      fTcalEvent -> fTPOEvent -> dump_event();
      ///
      /*if(CharmedEvent())
    	{
        std::cout << "Charm parent ID: " << fCharmParentID << std::endl;
        std::cout << "Charm name: " << fcharmname << std::endl;
        std::cout << "Charm charge: " << fCharmCharge << std::endl;
        std::cout << "Charm decay mode: " << fdecayMode << std::endl;
        std::cout << "Number of charged daughters: " << fnumChargedDaughters << std::endl;
        std::cout << "Number of neutral daughters: " << fnumNeutralDaughters << std::endl;  
	    }
        */
      if(ShortLivedParticleEvent())
      {
        std::cout << "Event " << ievent << " contains short-lived particles:" << std::endl;
        for (size_t i = 0; i < fSLPTypes.size(); i++) {
          std::string particleType = "";
          switch(fSLPTypes[i]) {
            case 1: particleType = "Charm"; break;
            case 2: particleType = "Tau"; break;
            case 3: particleType = "Other SLP"; break;
            default: particleType = "Unknown"; break;
          }
          std::cout << "  - " << particleType << ": " << fSLPNames[i] 
                   << " (Parent ID: " << fSLPParentIDs[i] << ")" << std::endl;
        }

      }

      fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
      // Modification based on Andre's changes
      //fPORecoEvent -> Reconstruct3DPS();
      //fPORecoEvent -> Dump();
      std::cout << "Start reconstruction of PORecs..." << std::endl;
      fPORecoEvent->ReconstructTruth();
       std::cout << "Start reconstruction of clusters..." << std::endl;
       fPORecoEvent->Reconstruct2DViewsPS();
       fPORecoEvent->Reconstruct3DPS_2();
       fPORecoEvent->ReconstructRearCals();
       // added for muon spectrometer
       std::cout << "Start reconstruction of muon spectrometer..." << std::endl;
       fPORecoEvent->ReconstructMuonSpectrometer();
       //fPORecoEvent->ReconstructMuonSpectrometerKasaKalman();

      fPORecoEvent->ReconstructClusters(0);   
      fPORecoEvent->ReconstructClusters(1);   


      std::cout << "DBSCAC in 3D" << std::endl;
      //fPORecoEvent->Reconstruct3DClusters();

      fPORecoEvent->Reconstruct3DPS_Eflow();
      fPORecoEvent->TrackReconstruct();
      fPORecoEvent->PSVoxelParticleFilter();
      //
      if(fPORecoEvent->GetPOFullRecoEvent()!=nullptr) 
	{
	  fVisibleEnergy = fPORecoEvent->GetPOFullRecoEvent()->TotalEvis();
	  fTotalEnergy = fPORecoEvent->GetPOFullRecoEvent()->TotalET();
	}
      fRearECalEnergy = fPORecoEvent->rearCals.rearCalDeposit;
      fRearHCalEnergy = fPORecoEvent->rearCals.rearHCalDeposit;
      fRearMuCalEnergy = fPORecoEvent->rearCals.rearMuCalDeposit;
      std::cout << fVisibleEnergy << " " 
		<< fTotalEnergy << " "
		<< fRearECalEnergy << " "
		<< fRearHCalEnergy << " "
		<< fRearMuCalEnergy << " "
		<< std::endl;
		
    // FastJet
    //std::cout << "##############################################" << std::endl;
    //std::cout << "Start FastJet reconstruction..." << std::endl;
    //fPORecoEvent->ReconstructJets();

    } else {
      std::cerr << "#########################################################" << std::endl;
      std::cerr << "File not found: " << file_path << std::endl;
      std::cerr << "#########################################################" << std::endl;
      if (fStatusBar) {
	      fStatusBar->GetBarPart(0)->SetBackgroundColor(kRed);
	      fStatusBar->SetText(Form("Warning: File not found for event #%d", ievent), 0);
      }
    }
  }
   //////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////
void FaserCalDisplay::PlotHistogramsFromRootFile() {
  //const std::string rootFilePath = "/home/hyperk/sw/FASERCAL/FASER_March2025/FASER/Display/build/AnalysisData.root";
  const std::string rootFilePath = "AnalysisData.root";
  
  // Check if the ROOT file exists
  TFile* inputFile = TFile::Open(rootFilePath.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ROOT file not found or is corrupted. Calling LoadAllEvents to create the ROOT file..." << std::endl;
    
    // Call LoadAllEvents to generate the ROOT file
    LoadAllEvents();
    
    // Try to reopen the file
    inputFile = TFile::Open(rootFilePath.c_str());
    if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error: Failed to create the ROOT file. Exiting." << std::endl;
      return;
    }
  }
  
  // Retrieve the tree from the ROOT file
  TTree* tree = (TTree*)inputFile->Get("MultiplicityTree");
  if (!tree) {
    std::cerr << "Error: Could not find tree in ROOT file. Exiting." << std::endl;
    inputFile->Close();
    return;
  }
  
  // Variables to read from the tree
  int totalMultiplicity, chargedMultiplicity, neutralMultiplicity, gammaMultiplicity, neutronMultiplicity;
  double neutrinoEnergy;
  std::string* mask_str = nullptr;
  tree->SetBranchAddress("TotalMultiplicity", &totalMultiplicity);
  tree->SetBranchAddress("ChargedMultiplicity", &chargedMultiplicity);
  tree->SetBranchAddress("NeutralMultiplicity", &neutralMultiplicity);
  tree->SetBranchAddress("GammaMultiplicity", &gammaMultiplicity);
  tree->SetBranchAddress("NeutronMultiplicity", &neutronMultiplicity);
  tree->SetBranchAddress("NeutrinoEnergy", &neutrinoEnergy);
  tree->SetBranchAddress("Mask", &mask_str);
  
  // Initialize histograms
  std::map<std::string, TH1F*> hMultiplicity;
  std::map<std::string, TH1F*> hMultiplicityC;
  std::map<std::string, TH1F*> hMultiplicity0;
  std::map<std::string, TH1F*> hMultiplicityGm;
  std::map<std::string, TH1F*> hMultiplicityN;
  std::map<std::string, TH1F*> hNeutrinoEnergy;
  
  std::vector<std::string> masks = {"nueCC", "numuCC", "nutauCC", "nuNC", "nuES"};
  for (const auto& mask : masks) {
    hMultiplicity[mask] = new TH1F(("hMultiplicity_" + mask).c_str(), ("Multiplicity " + mask + ";Multiplicity;Counts").c_str(), 50, -0.5, 49);
    hMultiplicityC[mask] = new TH1F(("hMultiplicityC_" + mask).c_str(), ("Charged Multiplicity " + mask + ";Multiplicity;Counts").c_str(), 50, -0.5, 49);
    hMultiplicity0[mask] = new TH1F(("hMultiplicity0_" + mask).c_str(), ("Neutral Multiplicity " + mask + ";Multiplicity;Counts").c_str(), 50, -0.5, 49);
    hMultiplicityGm[mask] = new TH1F(("hMultiplicityGm_" + mask).c_str(), ("Gamma Multiplicity " + mask + ";Multiplicity;Counts").c_str(), 50, -0.5, 49);
    hMultiplicityN[mask] = new TH1F(("hMultiplicityN_" + mask).c_str(), ("Neutron Multiplicity " + mask + ";Multiplicity;Counts").c_str(), 50, -0.5, 49);
    hNeutrinoEnergy[mask] = new TH1F(("hNeutrinoEnergy_" + mask).c_str(), ("Neutrino Energy " + mask + ";Energy (GeV);Counts").c_str(), 300, 0, 3000);
  }
  
  // Loop over tree entries and fill histograms
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    if (hMultiplicity.find(*mask_str) != hMultiplicity.end()) {
      hMultiplicity[*mask_str]->Fill(totalMultiplicity);
      hMultiplicityC[*mask_str]->Fill(chargedMultiplicity);
      hMultiplicity0[*mask_str]->Fill(neutralMultiplicity);
      hMultiplicityGm[*mask_str]->Fill(gammaMultiplicity);
      hMultiplicityN[*mask_str]->Fill(neutronMultiplicity);
      hNeutrinoEnergy[*mask_str]->Fill(neutrinoEnergy);
    }
  }
  
  // Plot histograms
  TCanvas* canvas = new TCanvas("canvas", "Multiplicity Histograms", 1200, 800);
  canvas->Divide(3, 2); // Create sub-pads for plotting
  int pad = 1;
  
  for (const auto& mask : masks) {
    canvas->cd(pad++);
    hMultiplicity[mask]->SetLineColor(kBlack);
    hMultiplicityC[mask]->SetLineColor(kBlue);
    hMultiplicity0[mask]->SetLineColor(kRed);
    hMultiplicityGm[mask]->SetLineColor(kOrange);
    hMultiplicityN[mask]->SetLineColor(kGreen);
    
    hMultiplicity[mask]->Draw("hist");
    hMultiplicityC[mask]->Draw("hist same");
    hMultiplicity0[mask]->Draw("hist same");
    hMultiplicityGm[mask]->Draw("hist same");
    hMultiplicityN[mask]->Draw("hist same");
    
    TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->AddEntry(hMultiplicity[mask], "Total", "l");
    legend->AddEntry(hMultiplicityC[mask], "Charged", "l");
    legend->AddEntry(hMultiplicity0[mask], "Neutral", "l");
    legend->AddEntry(hMultiplicityGm[mask], "Gamma", "l");
    legend->AddEntry(hMultiplicityN[mask], "Neutron", "l");
    legend->Draw();
  }
  
  TCanvas* energyCanvas = new TCanvas("energyCanvas", "Neutrino Energy Histograms", 1200, 800);
  energyCanvas->Divide(3, 2);
  pad = 1;
  
  for (const auto& mask : masks) {
    energyCanvas->cd(pad++);
    hNeutrinoEnergy[mask]->SetLineColor(kBlue);
    hNeutrinoEnergy[mask]->Draw("hist");
  }  
  // Clean up
  inputFile->Close();
  delete inputFile;
}
  ///////////////////////////////////////////
  ///////////////////////////////////////////
  //void FaserCalDisplay::DumpMuonInformations()
  //{
  //}

  ///////////////////////////////////////////
  //////////////////////////////////////////////////////////
void FaserCalDisplay::LoadAllEvents()
{
  std::ofstream fileOut("log.txt");    
  // Redirecting cout to write to "log.txt"
  std::cout.rdbuf(fileOut.rdbuf());
 
  const std::string& folder_path = "input/";
  std::vector<std::string> file_paths;
  
  // Get all file paths in the input directory
  for (const auto& entry : std::filesystem::directory_iterator(folder_path))
    {
      file_paths.push_back(entry.path().string());
    }
  // Limit the number of files to process
  size_t max_files = 100;
  if (file_paths.size() > max_files) {
    file_paths.resize(max_files); // Keep only the first 500 files
  }
  // Create ROOT file to save multiplicity data
  TFile* outputFile = new TFile("AnalysisData.root", "RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Error: Could not create output ROOT file." << std::endl;
    return;
  }
  TTree* tree = new TTree("AnalysisTree", "Storing analysis data");
  
  // Variables to store in the tree
  int irun, ievent, totalMultiplicity, chargedMultiplicity, neutralMultiplicity, gammaMultiplicity, neutronMultiplicity;
  int CharmType, CCNC, neutrinoType;
  double neutrinoEnergy;
  double prim_vx, prim_vy, prim_vz;

  std::string mask_str;
  tree->Branch("Run", &irun, "Run/I");
  tree->Branch("Event", &ievent, "Event/I");
  tree->Branch("Mask", &mask_str);
  tree->Branch("Vtx_pry",&prim_vx,"Vtx_pry/D");
  tree->Branch("Vty_pry",&prim_vy,"Vty_pry/D");
  tree->Branch("Vtz_pry",&prim_vz,"Vtz_pry/D");

  tree->Branch("TotalMultiplicity", &totalMultiplicity, "TotalMultiplicity/I");
  tree->Branch("ChargedMultiplicity", &chargedMultiplicity, "ChargedMultiplicity/I");
  tree->Branch("NeutralMultiplicity", &neutralMultiplicity, "NeutralMultiplicity/I");
  tree->Branch("GammaMultiplicity", &gammaMultiplicity, "GammaMultiplicity/I");
  tree->Branch("NeutronMultiplicity", &neutronMultiplicity, "NeutronMultiplicity/I");
  tree->Branch("NeutrinoEnergy", &neutrinoEnergy, "NeutrinoEnergy/D");
  tree->Branch("CharmType", &CharmType, "CharmType/I");
  tree->Branch("CCNC", &CCNC, "CCNC/I");
  tree->Branch("NeutrinoType", &neutrinoType, "NeutrinoType/I");
  tree->Branch("CharmProngs", &fnumChargedDaughters, "CharmProngs/I");
  tree->Branch("CharmDecayMode", &fdecayMode);
  tree->Branch("CharmParticle", &fcharmname);
  tree->Branch("CharmEnergy", &fCharmEnergy);
  tree->Branch("Vtx_dcy",&fdecay_vx,"Vtx_dcy/D");
  tree->Branch("Vty_dcy",&fdecay_vy,"Vty_dcy/D");
  tree->Branch("Vtz_dcy",&fdecay_vz,"Vtz_dcy/D");
  tree->Branch("FlightLength", &fdecayFlightLength, "FlightLength/D");
  tree->Branch("VisibleEnergy", &fVisibleEnergy,"VisibleEnergy/D");
  tree->Branch("TotalEnergy", &fTotalEnergy,"TotalEnergy/D");
  tree->Branch("RearECalEnergy", &fRearECalEnergy,"RearECalEnergy/D");
  tree->Branch("RearHCalEnergy", &fRearHCalEnergy,"RearHCalEnergy/D");
  tree->Branch("RearMuCalEnergy", &fRearMuCalEnergy,"RearMuCalEnergy/D");

  tree->Branch("TauType", &fTau, "TauType/I");
  tree->Branch("TauProngs", &fnumTauChargedDaughters, "TauProngs/I");
  tree->Branch("TauDecayMode", &ftauDecayMode);
  tree->Branch("TauParticle", &ftauname);
  tree->Branch("TauEnergy", &fTauEnergy, "TauEnergy/D");
  tree->Branch("TauVtx_dcy", &ftau_vx, "TauVtx_dcy/D");
  tree->Branch("TauVty_dcy", &ftau_vy, "TauVty_dcy/D");
  tree->Branch("TauVtz_dcy", &ftau_vz, "TauVtz_dcy/D");
  tree->Branch("TauFlightLength", &ftauDecayFlightLength, "TauFlightLength/D");
    
  // Generic SLP info
  tree->Branch("NumSLPs", &fSLPParentIDs);
  tree->Branch("SLPTypes", &fSLPTypes);
  tree->Branch("SLPNames", &fSLPNames);
 

  tree->Branch("MuTag_mom", &fmuTag_p);
  tree->Branch("MuTag_ene", &fmuTag_E);
  tree->Branch("MuTag_dist", &fmuTag_dist);
  tree->Branch("MuTag_muon", &fmuTag_muon);
  tree->Branch("MuTag_id", &fmuTag_id);

  TDatabasePDG* pdgDB = TDatabasePDG::Instance();
  AddCustomNucleusParticles();
  int cnt = 0;
  // Loop over each file and extract information
  std::cout << "Number of files " << file_paths.size() << std::endl;
  for (const std::string& file_path : file_paths)
    {
      std::cout << "Processing file: " << file_path << std::endl;
      cnt++;
      // Extract the base name from the file path
      std::string base_name = std::filesystem::path(file_path).stem().string();
      std::cout << "basename " << base_name << std::endl;    
      // Split the base_name using '_' as a delimiter
      std::istringstream ss(base_name);
      std::string token;
      std::vector<std::string> parts;
      while (std::getline(ss, token, '_'))
        {
	        parts.push_back(token);
        }
      if (parts.size() < 3)
        {
	        std::cerr << "Error: Invalid filename format, unable to parse: " << base_name << std::endl;
	        continue; // Skip this file if the format is incorrect
        }
      try {
	      mask_str = (parts.size() == 3) ? "NoMask" : parts[3];
        irun = std::stoi(parts[1]);
        ievent = std::stoi(parts[2]);
	std::cout << cnt << " Event Number: " << ievent << ", Run Number: " << irun << ", Mask: " << mask_str << std::endl;
	int imask = 0;
	if (mask_str == "nueCC") imask = 1;
	else if (mask_str == "numuCC") imask = 2;
	else if (mask_str == "nutauCC") imask = 3;
	else if (mask_str == "nuNC") imask = 4;
	else if (mask_str == "nuES") imask = 5;
	
	fTcalEvent = new TcalEvent();
	POevent = new TPOEvent();
	fTcalEvent->Load_event("input/", irun, ievent, imask, POevent);

	prim_vx = POevent->prim_vx.x();
	prim_vy = POevent->prim_vx.y();
	prim_vz = POevent->prim_vx.z();
	
	fdecay_vx = fdecay_vy = fdecay_vz = fdecayFlightLength = 0.0;
    ftau_vx = ftau_vy = ftau_vz = ftauDecayFlightLength = 0.0;

	totalMultiplicity = chargedMultiplicity = neutralMultiplicity = 0;
	gammaMultiplicity = neutronMultiplicity = 0;
	CharmType = neutrinoType = 0;
	//CharmedEvent();
  ShortLivedParticleEvent();
	// Count the particles based on their properties
	for (size_t i = 0; i < POevent->n_particles(); i++)
	  {
	    struct PO& aPO = POevent->POs[i];
	    TParticlePDG* particle = pdgDB->GetParticle(aPO.m_pdg_id);
	    int charge = particle ? particle->Charge() : 0;
	    
	    if (aPO.m_status == 1)  // Only consider final state particles
	      {
		if(IsCharmed(aPO.m_pdg_id))
		  CharmType = aPO.m_pdg_id;
		totalMultiplicity++;
		if (charge != 0) {
		  chargedMultiplicity++;
		} else {
		  neutralMultiplicity++;
		}
		if(aPO.m_pdg_id==22)
		  gammaMultiplicity++;
		if(aPO.m_pdg_id==2112)
		  neutronMultiplicity++;
	      }
	  }
	  // Get neutrino energy
        neutrinoEnergy = POevent->in_neutrino.m_energy;
	neutrinoType = POevent->in_neutrino.m_pdg_id;
	CCNC = POevent->isCC;

	fPORecoEvent = new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent);
	fPORecoEvent->ReconstructTruth();
	fPORecoEvent->Reconstruct2DViewsPS();
	fPORecoEvent->Reconstruct3DPS_2();
	fPORecoEvent->ReconstructRearCals();
	fPORecoEvent->ReconstructClusters(0); 
	fPORecoEvent->ReconstructClusters(1); 
	fPORecoEvent->Reconstruct3DPS_Eflow();
	fPORecoEvent->TrackReconstruct();
	fPORecoEvent->PSVoxelParticleFilter();
	//
	if(fPORecoEvent->GetPOFullRecoEvent()!=nullptr) 
	  {
	    fVisibleEnergy = fPORecoEvent->GetPOFullRecoEvent()->TotalEvis();
	    fTotalEnergy = fPORecoEvent->GetPOFullRecoEvent()->TotalET();
	  }
	fRearECalEnergy = fPORecoEvent->rearCals.rearCalDeposit;
	fRearHCalEnergy = fPORecoEvent->rearCals.rearHCalDeposit;
	fRearMuCalEnergy = fPORecoEvent->rearCals.rearMuCalDeposit;
	
	GetMuTagInfo();

	////////////////////////
	tree->Fill();
	delete POevent;
	delete fTcalEvent;
	delete fPORecoEvent;
      } catch (const std::invalid_argument& e) {
	std::cerr << "Error: Invalid number format in filename: " << base_name << " - " << e.what() << std::endl;
      } catch (const std::out_of_range& e) {
	std::cerr << "Error: Number out of range in filename: " << base_name << " - " << e.what() << std::endl;
      }
    }
  
  outputFile->cd();
  tree->Write();
  // Close output file
  outputFile->Close();
  std::cout << "Finished processing " << cnt << " files." << std::endl;
  DoExit();
}
  
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void FaserCalDisplay::DrawEvent(int irun, int ievent, int imask)
  {
    if (!gEve) 
    {
      std::cerr << "TEveManager is not initialized!" << std::endl;
      return;
    }
    gEve->GetViewers()->DeleteAnnotations();
    if(fHitElements)
      {
       	fHitElements->DestroyElements();
	      fPrimaryElements->DestroyElements();
	      fSecondaryShowerElements->DestroyElements();
        fSecondaryHadShowerElements->DestroyElements();
        fPixelHitElements->DestroyElements();
        fShortLivedParticleHitElements->DestroyElements();
        fMuonHitElements->DestroyElements();
        fProtonHitElements->DestroyElements();
        fPionHitElements->DestroyElements();
        fKaonHitElements->DestroyElements();
        fRearECALElements->DestroyElements();
        fRearHCALElements->DestroyElements();
        fRearMuCALElements->DestroyElements();
        fMuTagHitElements->DestroyElements();
        fMuSpectHitElements->DestroyElements();
        fMuonSpectFitElements->DestroyElements();
      }    
    TEveEventManager* currEvent = gEve->GetCurrentEvent();
    if( currEvent ) currEvent->DestroyElements();
    TEveElement* top = gEve->GetCurrentEvent();
    LoadEvent(irun, ievent, imask);
    //
    double wdx, wdy, wdz;
    if (TGeoBBox *box = dynamic_cast<TGeoBBox*>(gGeoManager->GetTopVolume()->GetShape()))
      {
	      wdx = box->GetDX();
	      wdy = box->GetDY();
	      wdz = box->GetDZ()/2.0;
	      printf("Top volume dimensions: DX = %f, DY = %f, DZ = %f\n", wdx, wdy, wdz);
      } else {
	      std::cerr << "Failed to get top volume dimensions." << std::endl;
	      exit(1);
      }
    TGeoShape *bigbox = new TGeoBBox("bigbox", wdx, wdy, wdz);
    TGeoMedium *air = gGeoManager->GetMedium("AIR");
    primary = new TGeoVolume("primary", bigbox, air);
    secondary_em = new TGeoVolume("secondary_em", bigbox, air);
    secondary_had = new TGeoVolume("secondary_had", bigbox, air);
    si_tracker = new TGeoVolume("si_tracker", bigbox, air);
    ShortLivedP = new TGeoVolume("ShortLivedP", bigbox, air);
    muonhit = new TGeoVolume("MuonHit", bigbox, air);
    pionhit = new TGeoVolume("PionHit", bigbox, air);
    protonhit = new TGeoVolume("ProtonHit", bigbox, air);
    kaonhit = new TGeoVolume("KaonHit", bigbox, air);
    // Adding rearcal, hcal and mucal
    rearecal = new TGeoVolume("RearECAL", bigbox, air);
    rearhcal = new TGeoVolume("RearHCAL", bigbox, air);
    rearmucal = new TGeoVolume("RearMuCAL", bigbox, air);

    //
    double voxelsize = fTcalEvent->geom_detector.fScintillatorVoxelSize/10.0;
    TGeoShape *box = new TGeoBBox("box", voxelsize/2.0,voxelsize/2.0,voxelsize/2.0);
    TGeoShape *trackerhitbox = new TGeoBBox("box", 0.1/2.0,0.1/2.0,0.1/2.0);
    //
    TGeoMaterial *matAluminum = new TGeoMaterial("Aluminum", 26.98, 13, 2.7);
    TGeoMedium *aluminum = new TGeoMedium("Aluminum", 2, matAluminum);
    //
    std::cout << " using copied digitized tracks " << fTcalEvent->getfTracks().size() << std::endl;
    //
    TEveElementList* hitList = new TEveElementList("Hits");
    TEveElementList* zoomhitList = new TEveElementList("ZoomHits");
    TEveElementList* muonhitList = new TEveElementList("MuonHits");
    TEveElementList* protonhitList = new TEveElementList("ProtonHits");
    TEveElementList* pionhitList = new TEveElementList("PionHits");
    TEveElementList* kaonhitList = new TEveElementList("KaonHits");

    TEveElementList* primaryList = new TEveElementList("PrimaryHits");
    TEveElementList* secondaryShowerList = new TEveElementList("SecondaryShowerHits");
    TEveElementList* secondaryHadShowerList = new TEveElementList("SecondaryHadShowerHits");
    TEveElementList* pixelhitList = new TEveElementList("PixelHits");
    TEveElementList* pixelRecoTrackList = new TEveElementList("PixelRecoTrackList");
    TEveElementList* ShortLivedParticleList = new TEveElementList("ShortLivedParticleHits");

    TEveElementList* rearecalList = new TEveElementList("RearECALHits");
    TEveElementList* rearhcalList = new TEveElementList("RearHCALHits");
    TEveElementList* rearmucalList = new TEveElementList("RearMuCALHits");

    /////////////////      
    double zoomRadius = 100.0;
    /////////////////
    for (const auto& track : fTcalEvent->getfTracks())
      {
        size_t nhits = track->fhitIDs.size();
        for ( size_t i = 0; i < nhits; i++)
          {
            long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
            
            // apply energy cut on scintillator voxel
            if(hittype == 0 && track->fEnergyDeposits[i] < 0.5)continue;
            // apply cut on pixel hit
            if(hittype == 1 && track->fEnergyDeposits[i] < 1e-3)continue;
            //
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            //
            // Create a translation matrix for the hit position
            TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, position.Y() / 10.0, position.Z() / 10.0);
            //
            auto currentTrack = track;
            //
            TGeoVolume* hitVolume = nullptr;
            TGeoVolume *hitVolumeTracker = nullptr;
            if(hittype == 0)
              {
                hitVolume = new TGeoVolume("HitVolume", box, air);
		            //hitVolume->SetFillColor(30); 
                //hitVolume->SetLineColor(kRed); 
                hitVolume->SetLineColorAlpha(kRed, 0.5);
                // Store the association between hitVolume and track
                volumeTrackMap[hitVolume] = track;
                //
                if(fabs(track->fPDG) == 11)
                  hitVolume->SetLineColorAlpha(kBlue, 0.5);
                //hitVolume->SetLineColor(kBlue); // electromagnetic is blue
                else if(fabs(track->fPDG) == 13){
                  //hitVolume->SetLineColor(kGreen); // muons
                  hitVolume->SetLineColorAlpha(kGreen,0.5); // muons
                  muonhit->AddNode(hitVolume,i,trans);
                  std::cout << "muon voxel " << track->fPDG << " " << track->fparentID << " " << track->fprimaryID << " " << track->ftrackID << ""<< fCharmParentID << " " << fTauParentID << std::endl;
                }
                else if(abs(track->fPDG) == 211){
                  hitVolume->SetLineColorAlpha(kOrange,0.5); // pions
                  pionhit->AddNode(hitVolume,i,trans);
                }
                else if(fabs(track->fPDG) == 2212){
                  hitVolume->SetLineColorAlpha(kYellow,0.5); // protons
                  protonhit->AddNode(hitVolume,i,trans);
                }
                else if(fabs(track->fPDG) == 321){
                  hitVolume->SetLineColorAlpha(kMagenta,0.5); // kaons
                  kaonhit->AddNode(hitVolume,i,trans);
                }
                // get primary tracks
                if(track->fparentID == 0)
                {
                  primary->AddNode(hitVolume, i, trans);
                  if(IsCharmed(track->fPDG)||abs(track->fPDG) == 15)
                    {
                      ShortLivedP->AddNode(hitVolume, i, trans);
                      std::cout << "short lived particle voxel " << track->fPDG << std::endl;
                    }
                }
              else 
               {
                  if(fabs(track->fPDG) == 11)
                    secondary_em->AddNode(hitVolume, i, trans);
                  else
                    secondary_had->AddNode(hitVolume, i, trans); 
                  
                    if(track->fprimaryID == fCharmParentID || track->fparentID == fTauParentID)
                  {
                    ShortLivedP->AddNode(hitVolume, i, trans);
                    std::cout << "short lived particle voxel 2 --- " << track->fPDG << std::endl;
                  }
                }
              }
      	    else if (hittype == 1)
	          {
      		    hitVolumeTracker = new TGeoVolume("TrackerHitVolume", trackerhitbox, air);
		          hitVolumeTracker->SetLineColor(kMagenta); 
		          si_tracker->AddNode(hitVolumeTracker, i, trans);
	          }
	          else
	          {
      		    std::cout << " Unknown type of hit " << std::endl;
		          delete trans;
	          }
            if (hitVolume)
	          {
              TEveGeoShape* eveShape = new TEveGeoShape(hitVolume->GetName());
              eveShape->SetShape(hitVolume->GetShape());
              eveShape->SetMainColor(hitVolume->GetLineColor());
              eveShape->SetTransMatrix(*trans);
              eveShape->SetPickable(kTRUE); // Enable selection
              eveShape->SetUserData((void*)&track);
              //std::cout << "Created shape: " << eveShape << ", pickable: " << eveShape->IsPickable() << std::endl;
              hitList->AddElement(eveShape);
              // Store the association between shape and track
              //shapeTrackMap[eveShape] = currentTrack;
	          }
            if (hitVolumeTracker)
	          {
              TEveGeoShape* eveShapeP = new TEveGeoShape(hitVolumeTracker->GetName());
              eveShapeP->SetShape(hitVolumeTracker->GetShape());
              eveShapeP->SetMainColor(hitVolumeTracker->GetLineColor());
              eveShapeP->SetTransMatrix(*trans);
              eveShapeP->SetPickable(kTRUE); // Enable selection
              pixelhitList->AddElement(eveShapeP);
              // Store the association between shape and track
              //shapeTrackMap[eveShapeP] = currentTrack;
	          }
	        }
      }
      //////////////
      if (primary) 
      {
        std::cout << "...... " << primary->GetNdaughters() << std::endl;
        for (int i = 0; i < primary->GetNdaughters(); ++i) 
        {
	        TGeoNode* node = primary->GetNode(i);
	        if (!node) continue;  
	        TGeoVolume* vol = node->GetVolume();
	        TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	        if (vol && trans) 
          {
    	      TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	          eveShape->SetShape(vol->GetShape());
	          eveShape->SetMainColor(vol->GetLineColor());
	          eveShape->SetTransMatrix(*trans);
	          primaryList->AddElement(eveShape);
    	    }
        }
        fPrimaryElements->AddElement(primaryList);
      }
      ///////////////
      if (muonhit) 
      {
        std::cout << "...... " << muonhit->GetNdaughters() << std::endl;
        for (int i = 0; i < muonhit->GetNdaughters(); ++i) 
        {
	        TGeoNode* node = muonhit->GetNode(i);
	        if (!node) continue;  
	        TGeoVolume* vol = node->GetVolume();
	        TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	        if (vol && trans) 
          {
      	    TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	          eveShape->SetShape(vol->GetShape());
	          eveShape->SetMainColor(vol->GetLineColor());
	          eveShape->SetTransMatrix(*trans);
	          muonhitList->AddElement(eveShape);
	        }
        }
        fMuonHitElements->AddElement(muonhitList);
      }
      ///////////////
      if (pionhit) 
      {
        std::cout << "...... " << pionhit->GetNdaughters() << std::endl;
        for (int i = 0; i < pionhit->GetNdaughters(); ++i) 
        {
          TGeoNode* node = pionhit->GetNode(i); 
          if (!node) continue;
          TGeoVolume* vol = node->GetVolume();
          TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
          if (vol && trans) 
          {
            TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
            eveShape->SetShape(vol->GetShape());
            eveShape->SetMainColor(vol->GetLineColor());
            eveShape->SetTransMatrix(*trans);
            pionhitList->AddElement(eveShape);
          }
        } 
        fPionHitElements->AddElement(pionhitList);
      }
      ///////////////
      if (kaonhit) 
      {
        std::cout << "...... " << kaonhit->GetNdaughters() << std::endl;
        for (int i = 0; i < kaonhit->GetNdaughters(); ++i) 
        {
          TGeoNode* node = kaonhit->GetNode(i);
          if (!node) continue;
          TGeoVolume* vol = node->GetVolume();
          TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
          if (vol && trans) 
          {
            TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
            eveShape->SetShape(vol->GetShape());
            eveShape->SetMainColor(vol->GetLineColor());
            eveShape->SetTransMatrix(*trans);
            kaonhitList->AddElement(eveShape);
          }
        }
        fKaonHitElements->AddElement(kaonhitList);
      }
      ///////////////
    if (protonhit) 
    {
      std::cout << "...... " << protonhit->GetNdaughters() << std::endl;
      for (int i = 0; i < protonhit->GetNdaughters(); ++i) 
      {
        TGeoNode* node = protonhit->GetNode(i);
        if (!node) continue;
        TGeoVolume* vol = node->GetVolume();
        TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
        if (vol && trans) 
        {
          TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
          eveShape->SetShape(vol->GetShape());
          eveShape->SetMainColor(vol->GetLineColor());
          eveShape->SetTransMatrix(*trans);
          protonhitList->AddElement(eveShape);
        }
      }
      fProtonHitElements->AddElement(protonhitList);
    }
    ///////////////
    if (secondary_em) 
    {
      std::cout << "...... " << secondary_em->GetNdaughters() << std::endl;
      for (int i = 0; i < secondary_em->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = secondary_em->GetNode(i);
	      if (!node) continue;  
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	      if (vol && trans) 
        {
	        TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        eveShape->SetShape(vol->GetShape());
	        eveShape->SetMainColor(vol->GetLineColor());
	        eveShape->SetTransMatrix(*trans);
	        secondaryShowerList->AddElement(eveShape);
        } 
      }
      fSecondaryShowerElements->AddElement(secondaryShowerList);
    }
    //////////////
    if (secondary_had) 
    {
      std::cout << "...... " << secondary_had->GetNdaughters() << std::endl;
      for (int i = 0; i < secondary_had->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = secondary_had->GetNode(i);
	      if (!node) continue;  
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	      if (vol && trans) 
        {
	        TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        eveShape->SetShape(vol->GetShape());
	        eveShape->SetMainColor(vol->GetLineColor());
	        eveShape->SetTransMatrix(*trans);
	        secondaryHadShowerList->AddElement(eveShape);
	      }
      }
      fSecondaryHadShowerElements->AddElement(secondaryHadShowerList);
    }
    //////////////
    if (si_tracker) 
    {
      std::cout << "...xx... " << si_tracker->GetNdaughters() << std::endl;
      for (int i = 0; i < si_tracker->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = si_tracker->GetNode(i);
	      if (!node) continue;
	      //std::cout<< "here i am!! " << std::endl;
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	      if (vol && trans) 
        {
	        TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        eveShape->SetShape(vol->GetShape());
	        eveShape->SetMainColor(vol->GetLineColor());
	        eveShape->SetTransMatrix(*trans);
	        pixelhitList->AddElement(eveShape);
	      }
      }
      fPixelHitElements->AddElement(pixelhitList);
    }
    //////////////
    if (ShortLivedP) 
    {
      std::cout << "...xx... " << ShortLivedP->GetNdaughters() << std::endl;
      for (int i = 0; i < ShortLivedP->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = ShortLivedP->GetNode(i);
	      if (!node) continue;
	      //std::cout<< "here i am!! " << std::endl;
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	      if (vol && trans) 
        {
	        TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        eveShape->SetShape(vol->GetShape());
	        eveShape->SetMainColor(vol->GetLineColor());
	        eveShape->SetTransMatrix(*trans);
	        ShortLivedParticleList->AddElement(eveShape);
	      }
      }
      fShortLivedParticleHitElements->AddElement(ShortLivedParticleList);
    }   
    //////////////
    fHitElements->AddElement(hitList);
    // Draw reconstructed tracks
    GetRecoTracks();
    //////
    // rearECAL
    for (const auto &it : fTcalEvent->rearCalDeposit)
    {
	    std::cout << "ïnside RearECAL" << std::endl;
	    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearCal(it.moduleID);
      double zBox = it.energyDeposit / 1e2; // 1cm is 1 GeV
	    std::cout << position << std::endl;
	    std::cout << it.energyDeposit << " " << zBox << std::endl;
	    std::cout << it.moduleID << std::endl;
   
      TGeoShape *box = new TGeoBBox("RearCALBox", fTcalEvent->geom_detector.rearCalSizeX / 20.0,
                                    fTcalEvent->geom_detector.rearCalSizeY / 20.0, zBox / 20.0);
      TGeoVolume *hitVolume = new TGeoVolume("RearCalVolume", box, air);
      hitVolume->SetLineColor(kBlue);
      TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                  position.Y() / 10.0, (position.Z() + zBox / 2.0) / 10.0);
      rearecal->AddNode(hitVolume, it.moduleID, trans);
    }
    // rearHCAL
    for (const auto &it : fTcalEvent->rearHCalDeposit)
    {
	    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearHCal(it.moduleID);
      double zBox = it.energyDeposit*10; // 1mm is 100 MeV                 
	    TGeoShape *box = new TGeoBBox("rearHCALbox", fTcalEvent->geom_detector.rearHCalSizeX / 20.0,
                                      zBox / 20.0, 
                                      fTcalEvent->geom_detector.rearHCalSizeZ / 20.0);
      TGeoVolume *hitVolume = new TGeoVolume("RearHCalVolume", box, air);
      hitVolume->SetLineColor(kRed);
      TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
                                                  (position.Y() + zBox / 2.0) / 10.0, position.Z() / 10.0);
      rearhcal->AddNode(hitVolume, it.moduleID, trans);
    }
    // rearMuCAL
   /*
    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZRearHCal(10);
    double zBox = fTcalEvent->rearMuCalDeposit*10.0; // 1mm is 100 MeV
    TGeoShape *mubox = new TGeoBBox("rearmucalbox", fTcalEvent->geom_detector.rearHCalSizeX / 20.0,
				  zBox / 20.0, 
				  fTcalEvent->geom_detector.rearHCalSizeZ / 20.0);
    TGeoVolume *hitVolume = new TGeoVolume("RearMuCalVolume", mubox, air);
    hitVolume->SetLineColor(kMagenta);
    TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0,
						 (position.Y()+zBox/2.0) / 10.0, position.Z() / 10.0);
    rearmucal->AddNode(hitVolume, 0, trans);
   */
    for (const auto *mt : fTcalEvent->fMuTagTracks)
    {
      if (!mt) continue;
      int hitsize = static_cast<int>(mt->pos.size());
      if (hitsize == 0) continue;
      std::cout << "ïnside RearMuCAL" << std::endl;
      for (const auto &hit : mt->pos)
      {
        // ROOT::Math::XYZVector accessors
        double x = hit.X() / 10.0;
        double y = hit.Y() / 10.0;
        double z = hit.Z() / 10.0;
        std::cout << x << " " << y << " " << z << std::endl;
        TGeoShape *mubox = new TGeoBBox("rearmucalbox", 0.5, 0.5, 0.5);
        TGeoVolume *hitVolume = new TGeoVolume("RearMuCalVolume", mubox, air);
        hitVolume->SetLineColor(kMagenta);
        TGeoTranslation *trans = new TGeoTranslation(x, y, z);
        rearmucal->AddNode(hitVolume, 0, trans);
      }
    }
    ////////////////////////
    // Draw rearcal, hcal and mucal
    if (rearecal) 
    {
      std::cout << "...... " << rearecal->GetNdaughters() << std::endl;
      for (int i = 0; i < rearecal->GetNdaughters(); ++i) 
      {
      	TGeoNode* node = rearecal->GetNode(i);
	      if (!node) continue;  
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
	      if (vol && trans) 
        {
	        TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        eveShape->SetShape(vol->GetShape());
	        eveShape->SetMainColor(vol->GetLineColor());
	        eveShape->SetTransMatrix(*trans);
	        rearecalList->AddElement(eveShape);
	      }
      }
      fRearECALElements->AddElement(rearecalList);
    }
    ///////////////////////
    if (rearhcal) 
    {
      std::cout << "...... " << rearhcal->GetNdaughters() << std::endl;
      for (int i = 0; i < rearhcal->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = rearhcal->GetNode(i);
	      if (!node) continue;
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
        if (vol && trans) 
        {
          TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
          eveShape->SetShape(vol->GetShape());
          eveShape->SetMainColor(vol->GetLineColor());
          eveShape->SetTransMatrix(*trans);
          rearhcalList->AddElement(eveShape);
        }
      }
      fRearHCALElements->AddElement(rearhcalList);
    }
    ///////////////////////
    if (rearmucal) 
    {
      std::cout << "...... " << rearmucal->GetNdaughters() << std::endl;
      for (int i = 0; i < rearmucal->GetNdaughters(); ++i) 
      {
	      TGeoNode* node = rearmucal->GetNode(i);
	      if (!node) continue;
	      TGeoVolume* vol = node->GetVolume();
	      TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
        if (vol && trans) 
        {
          TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
          eveShape->SetShape(vol->GetShape());
          eveShape->SetMainColor(vol->GetLineColor());
          eveShape->SetTransMatrix(*trans);
          rearmucalList->AddElement(eveShape);
        }
      }
      fRearMuCALElements->AddElement(rearmucalList);
    }
    ////////////////////////
    // Draw MuTag hits
    TGeoShape *mutagbox = new TGeoBBox("mutagbox", 1, 1, 1); // to be defined the cell size now                                                                               
    TGeoVolume* MuTagHit = new TGeoVolume("MuTagHit", mutagbox, air);
    TEveElementList* muTagHitList = new TEveElementList("MuTagHits");
    TEveElementList* muonSpectList = new TEveElementList("MuonSpect");
    for (const auto &it : fTcalEvent->fMuTagTracks)
    {
	    for(size_t i = 0; i < it->pos.size(); ++i)
	    {
	      const auto position = it->pos[i];
	      const auto momentum = it->mom[i];

        std::cout << "MuSpectrometer hit position " 
            << ievent << " "
            << it->ftrackID << " "
            << it->fPDG << " "
            << it->layerID[i] << " "
            << momentum.x() << " "
            << momentum.y() << " "
            << momentum.z() << " "
            << position.x() << " "
            << position.y() << " "
            << position.z() << " "
            << std::endl;

	      TGeoTranslation* trans = new TGeoTranslation(position.x()/10.0,
							 position.y()/10.0,
							 position.z()/10.0);
	      TGeoVolume* vol = new TGeoVolume("MuTagHit",mutagbox,air);
	      vol->SetLineColor(kOrange);
	      vol->SetLineColorAlpha(kOrange, 0.7);
	      if (abs(it->fPDG)==13)
	      {
		      vol->SetLineColor(kGreen);
		      vol->SetLineColorAlpha(kGreen, 0.7);
          
	      }
  	    TEveGeoShape* eveShape = new TEveGeoShape(Form("MuTagHit_%d_%lu",it->ftrackID,i));
	      eveShape->SetShape(vol->GetShape());
	      eveShape->SetMainColor(vol->GetLineColor());
	      eveShape->SetTransMatrix(*trans);
	      eveShape->SetPickable(kTRUE);
	      eveShape->SetUserData((void*)it);
        if (abs(it->fPDG)==13)
          muonSpectList->AddElement(eveShape);
	      muTagHitList->AddElement(eveShape);
	    }
    }
    fMuTagHitElements->AddElement(muTagHitList);
    fMuSpectHitElements->AddElement(muonSpectList);
    ////////////////////////
    // i want to check MuTag muons::
    GetMuTagInfo();
    //
    GetClusterInfo();
    //
    Get3DClusterInfo();
    //  
    GetReconstructedVoxels(bigbox, air, box);
    //
    GetMuonSpectrometerInfo(bigbox, air, box);
    //
    gEve->AddGlobalElement(fHitElements);
    gEve->AddGlobalElement(fMuTagHitElements);
    //
    gEve->FullRedraw3D(kTRUE);
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::GetMuonSpectrometerInfo(TGeoShape* bigbox, TGeoMedium* air, TGeoShape* box)
  {
    muonSpectVol = new TGeoVolume("MuonSpectVol", bigbox, air);
    TEveElementList* muonSpectFitList = new TEveElementList("MuonSpectFitList");

    std::cout << " Inside Muon Spectrometer " << std::endl;
    for (int i = 0; i < fPORecoEvent->fMuTracks.size(); ++i)
      {
        TMuTrack* muTrack = &fPORecoEvent->fMuTracks[i];
        if (!muTrack) continue;
        std::cout << "Muon track info: "
                  << muTrack->fcharge << " "
                  << muTrack->fpx << " "
                  << muTrack->fpy << " "
                  << muTrack->fpz << " "
                  << muTrack->fp << " "
                  << muTrack->ftrackID << " "
                  << muTrack->fPDG << " "
                  << std::endl;
        for (size_t j=0; j<muTrack->fpos.size(); ++j)
        {
          const auto position = muTrack->fpos[j];
          std::cout << j << " Muon Spectrometer hit position " 
                    << muTrack->ftrackID << " "
                    << muTrack->layerID[j] << " "
                    << position.x() << " "
                    << position.y() << " "
                    << position.z() << " "
                    << std::endl;

          TGeoTranslation* trans = new TGeoTranslation(position.x()/10.0,
				                        			  position.y()/10.0,
							                          position.z()/10.0);
	        TGeoVolume* vol = new TGeoVolume("MuSpectFitHit", box, air);
	        vol->SetLineColor(kCyan);
	        vol->SetLineColorAlpha(kCyan, 0.7);
          muonSpectVol->AddNode(vol, j, trans);
        }
      }
    if (muonSpectVol)
      {
        for (int i = 0; i < muonSpectVol->GetNdaughters(); ++i) 
        {
          TGeoNode* node = muonSpectVol->GetNode(i);
          if (!node) continue;
          TGeoVolume* vol = node->GetVolume();
          TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
          if (vol && trans) 
          {
            TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	          eveShape->SetShape(vol->GetShape());
            eveShape->SetMainColor(vol->GetLineColor());
            eveShape->SetTransMatrix(*trans);
            muonSpectFitList->AddElement(eveShape);
          }
      }
      fMuonSpectFitElements->AddElement(muonSpectFitList);
    }
  }
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::GetMuTagInfo()
  {
    fmuTag_p.clear();
    fmuTag_E.clear();
    fmuTag_id.clear();
    fmuTag_px.clear();
    fmuTag_py.clear();
    fmuTag_pz.clear();
    fmuTag_alltrk.clear();
    fmuTag_muonSign.clear();

    fmuTag_muon = 0;
    fmuTag_dist.clear();
    
    std::cout << " Inside MuTag " << std::endl;
    std::vector<ROOT::Math::XYZVector> muons;
    std::vector<ROOT::Math::XYZVector> muons_e;
    for (const auto &it : fTcalEvent->fMuTagTracks)
      {
	std::cout << it->ftrackID << " " 
		  << it->fPDG << " " 
		  << it->mom[0].x() << " "
		  << it->mom[0].y() << " "
		  << it->mom[0].z() << " "
		  << it->pos[0].x() << " "
		  << it->pos[0].y() << " "
		  << it->pos[0].z() << " "
		  << std::endl;

	fmuTag_alltrk.push_back(it->fPDG);
	
	// Extract kinematics if it's a muon and has data
	if (abs(it->fPDG) == 13 && !it->mom.empty()) 
	  {
	    ROOT::Math::XYZVector mom = it->mom.front();
	    ROOT::Math::XYZVector pos = it->pos.front();
	    double px = mom.X();
	    double py = mom.Y();
	    double pz = mom.Z();
	    double p = mom.R();
	    double mass_muon = 105.658; // MeV
	    double energy = sqrt(p * p + mass_muon * mass_muon);
	    double slope_xz = (pz != 0.0) ? px / pz : 999.0;
	    double slope_yz = (pz != 0.0) ? py / pz : 999.0;
	    double theta = atan2(sqrt(px*px + py*py), pz) * 180.0 / M_PI;
	    int sign = (it->fPDG > 0) ? 1 : -1;
	    fmuTag_muonSign.push_back(sign);

	    fmuTag_p.push_back(p);
	    fmuTag_E.push_back(energy);
	    fmuTag_id.push_back(it->ftrackID);
	    fmuTag_px.push_back(px);
	    fmuTag_py.push_back(py);
	    fmuTag_pz.push_back(pz);
	    
	    std::cout << "  --> Muon kinematics:"
		      << " |p| = " << p << " MeV"
		      << ", E = " << energy << " MeV"
		      << ", slope_xz = " << slope_xz
		      << ", slope_yz = " << slope_yz
		      << ", theta = " << theta << " deg"
		      << " muon sign " << sign << " "
		      << " " << it->ftrackID << " "
		      << std::endl; 
	    muons.push_back(pos);
	    muons_e.push_back(mom);
	  }
      }
    fmuTag_muon = muons.size();
    /// Check muons arriving to MuTag
    if (muons.size()<2)
      {
	std::cout << "Less than 2 muons reached to MuTag " << muons.size() << std::endl;
      } else {
      
      std::cout << "Multi-muons " << muons.size() << std::endl;
      for (size_t i = 0; i < muons.size(); i++)
	{
	  for(size_t j = i+1; j < muons.size(); j++)
	    {
	      double dist = (muons[i] - muons[j]).R(); //magnitude in 3D
	      std::cout << "Distance between muons " << dist << " mm " << std::endl;
	      fmuTag_dist.push_back(dist);
	    }
	}
    }
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::GetRecoTracks()
  {
    // Create a new TEveElementList for reconstructed tracks
    fPixelRecoTrackElements = new TEveElementList("Reconstructed Tracks");
    
    for (auto &itrk : fPORecoEvent->fTKTracks)
      {
	int nhits = itrk.tkhit.size();
	if (nhits == 0)
	  continue;
	
	TEveLine *trackLine = new TEveLine(nhits);
	trackLine->SetName(Form("Track_%d", itrk.trackID));
	trackLine->SetMainColor(kBlack);
	trackLine->SetLineWidth(2);
	
	// Add points to the TEveLine
	for (auto &hit : itrk.tkhit)
	  {
	    double x = hit.point.x() / 10.0;
	    double y = hit.point.y() / 10.0;
	    double z = hit.point.z() / 10.0;
	    trackLine->SetNextPoint(x, y, z);
	  }
	
	// Attach metadata: track ID
	trackLine->SetElementName(Form("TrackID: %d", itrk.trackID));
	trackLine->SetUserData((void *)&itrk);  // Store a pointer to the track
	
	fPixelRecoTrackElements->AddElement(trackLine);
      }
  }

  void FaserCalDisplay::GetReconstructedVoxels(TGeoShape *bigbox,TGeoMedium *air, TGeoShape *box)
  {
    voxVol = new TGeoVolume("ps_reco_vox", bigbox, air); 
    TEveElementList *voxHitList = new TEveElementList("RecoVoxelHits");
    TEveElementList *voxGhostList = new TEveElementList("GhostVoxelHits");
    int i = 0;
    for(auto &vox : fPORecoEvent->PSvoxelmap)
    {
	    i++;
	    long ID = vox.first;
	    double ehit = vox.second.RawEnergy;
	    if (ehit < 0.5) continue;
	    ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(ID);
	    TGeoTranslation *trans = new TGeoTranslation(position.X() / 10.0, position.Y() / 10.0, position.Z() / 10.0);
	    TGeoVolume *hitVolume;// = new TGeoVolume("RecoHitVolume", box, air);
	    if(vox.second.ghost)
	    {
	      hitVolume = new TGeoVolume("GhostHitVolume", box, air);
	      hitVolume->SetLineColorAlpha(kGray,0.5); // ghost
	    } 
	    else
	    {
	      hitVolume = new TGeoVolume("RecoHitVolume", box, air);
	      hitVolume->SetLineColorAlpha(kOrange,0.5); 	    
	    }
	    voxVol->AddNode(hitVolume,i,trans);
    }
    if (voxVol) 
    {
      //std::cout << "...vv... " << voxVol->GetNdaughters() << std::endl;
      for (int i = 0; i < voxVol->GetNdaughters(); ++i) {
        TGeoNode* node = voxVol->GetNode(i);
        if (!node) continue;
        TGeoVolume* vol = node->GetVolume();
        TGeoTranslation* trans = dynamic_cast<TGeoTranslation*>(node->GetMatrix());
        if (vol && trans) 
        {
          TEveGeoShape* eveShape = new TEveGeoShape(vol->GetName());
	        //std::cout << "... Ghost:: " << vol->GetName() << std::endl;
          eveShape->SetShape(vol->GetShape());
          eveShape->SetMainColor(vol->GetLineColor());
          eveShape->SetTransMatrix(*trans);
          voxHitList->AddElement(eveShape);
	        if(std::string(vol->GetName())=="GhostHitVolume")
	        {
	          voxGhostList->AddElement(eveShape);
	        }
        }
      }
      std::cout << "voxGhostList has " << voxGhostList->NumChildren() << " children" << std::endl;

      fVoxHitElements->AddElement(voxHitList);
      fVoxGhostElements->AddElement(voxGhostList);
    }
    
  }
  



  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /*
  void FaserCalDisplay::ShowClusterHits()
  {
    TGeoMaterial *matAir = new TGeoMaterial("Air", 0, 0, 0);
    TGeoMedium *air = new TGeoMedium("Air", 1, matAir);
    TGeoShape *hitShape = new TGeoBBox("ClusterHit", 0.5, 0.5, 0.5);

    TEveElementList *clusterHitList = new TEveElementList("AllClusterHits");
    int colorOffset = 0;
    
    for (const auto &c : fPORecoEvent->GetPSClusters(0)) // XZ view
      {
	int color = TColor::GetColorPalette(colorOffset % TColor::GetNumberOfColors());
	colorOffset += 5;
	
	for (size_t i = 0; i < c.hits.size(); ++i)
	  {
	    const auto &hit = c.hits[i];
	    double x, y, z;
	    fPORecoEvent->pshit2d_position(hit.id, x, y, z);

	    TGeoTranslation *trans = new TGeoTranslation(x/10.0, y/10.0, z/10.0);
	    TGeoVolume *vol = new TGeoVolume("ClusterHit", hitShape, air);
	    vol->SetLineColor(color);
	    vol->SetLineColorAlpha(color, 0.6);

	    TEveGeoShape *eveShape = new TEveGeoShape(Form("ClusterHit_%d_%lu",c.clusterID, i));
	    eveShape->SetShape(vol->GetShape());
	    eveShape->SetMainColor(vol->GetLineColor());
	    eveShape->SetTransMatrix(*trans);
	    eveShape->SetPickable(kTRUE);
	    eveShape->SetTitle(Form("Cluster %d, Hit %lu, E=%.2f MeV", c.clusterID, i, hit.EDeposit));
	    clusterHitList->AddElement(eveShape);

	  }

      }
    fClusterHitElements->AddElement(clusterHitList);
  }
  */
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /*
    void FaserCalDisplay::GetClusterInfo()
  {
    TGeoMaterial *matAir = new TGeoMaterial("Air", 0, 0, 0);
    TGeoMedium *air = new TGeoMedium("Air", 1, matAir);
    TGeoShape *hitShape = new TGeoBBox("ClusterHit", 0.5, 0.5, 0.5);
    
    TEveElementList *allClusterHits = new TEveElementList("AllClusterHits");
    
    for (int view = 0; view <= 1; ++view) // 0 = XZ, 1 = YZ
      {
        const char* viewName = (view == 0) ? "XZ" : "YZ";
        TEveElementList *clusterHitList = new TEveElementList(Form("%s Clusters", viewName));
	
        int colorOffset = (view == 0) ? 0 : 100;
	
        for (const auto &c : fPORecoEvent->GetPSClusters(view))
	  {
            int color = TColor::GetColorPalette((colorOffset + c.clusterID * 3) % TColor::GetNumberOfColors());
	    
            for (size_t i = 0; i < c.hits.size(); ++i)
	      {
                const auto &hit = c.hits[i];
                double x, y, z;
                fPORecoEvent->pshit2d_position(hit.id, x, y, z);
		
                TGeoTranslation *trans = new TGeoTranslation(x / 10.0, y / 10.0, z / 10.0);
                TGeoVolume *vol = new TGeoVolume("ClusterHit", hitShape, air);
                vol->SetLineColor(color);
                vol->SetLineColorAlpha(color, 0.6);
		
                TEveGeoShape *eveShape = new TEveGeoShape(Form("ClusterHit_%s_%d_%lu", viewName, c.clusterID, i));
                eveShape->SetShape(vol->GetShape());
                eveShape->SetMainColor(vol->GetLineColor());
                eveShape->SetTransMatrix(*trans);
                eveShape->SetPickable(kTRUE);
                eveShape->SetTitle(Form("View %s, Cluster %d, Hit %lu, E=%.2f MeV", viewName, c.clusterID, i, hit.EDeposit));
		
                clusterHitList->AddElement(eveShape);
	      }
	  }
	
        allClusterHits->AddElement(clusterHitList);
      }
    
    fClusterHitElements->AddElement(allClusterHits);
  }
  */
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::GetClusterInfo()
  {
    CleanCanvas();
    TCanvas *myCan = CreateCanvas("Clusters",1);
    gPad->Update();
    myCan->Divide(2,1);
    const double boxSize = 10.0; // 1 cm = 10 mm
    
    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nztot = fTcalEvent->geom_detector.NRep * nzlayer;
    
    std::vector<int> colors = {kRed+1, kBlue+1, kGreen+2, kMagenta+2,kOrange, kCyan, kViolet, kTeal+1};
    
    // ========= XZ VIEW =========
    myCan->cd(1);

    TH2D *xviewPS = new TH2D("xviewPS", "Scintillator xz-view", nztot, 0, nztot, nx, 0, nx);
    TH2D *yviewPS = new TH2D("yviewPS", "Scintillator yz-view", nztot, 0, nztot, ny, 0, ny);
    gStyle->SetOptStat(0);

    const auto& clustersXZ = fPORecoEvent->GetPSClusters(0);
    int colorIndex = 0;


    for (const auto& cluster : clustersXZ)
    {
      //int color = colors[colorIndex++ % colors.size()];
        for (const auto& hit : cluster.hits)
        {
            double x, y, z;
            fPORecoEvent->pshit2d_position(hit.id, x, y, z);
	    // std::cout << "........ " << cluster.clusterID << " " << hit.id << " " << x << " " << y << " " << z << std::endl;  
	    xviewPS->Fill(z,x,cluster.clusterID);
        }
    }
    xviewPS->Draw();

    // ========= ZY VIEW =========
    myCan->cd(2);
    gStyle->SetOptStat(0);

    const auto& clustersZY = fPORecoEvent->GetPSClusters(1);
    colorIndex = 0;

    for (const auto& cluster : clustersZY)
    {
      //int color = colors[colorIndex++ % colors.size()];
        for (const auto& hit : cluster.hits)
        {
            double x, y, z;
            fPORecoEvent->pshit2d_position(hit.id, x, y, z);
	    //std::cout << "xxxxxx " << cluster.clusterID << " " << hit.id << " " << x << " " << y << " " << z << std::endl;  
	    yviewPS->Fill(z,y,cluster.clusterID);
        }
    }
    yviewPS->Draw();
    myCan->Update();
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::Get3DClusterInfo()
  {
    //CleanCanvas();
    TCanvas *myCan = CreateCanvas("Clusters3D",1);
    gPad->Update();
    myCan->cd();
    //    const double boxSize = 10.0; // 1 cm = 10 mm
    
    int nx = fTcalEvent->geom_detector.fScintillatorSizeX / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int ny = fTcalEvent->geom_detector.fScintillatorSizeY / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nzlayer = fTcalEvent->geom_detector.fSandwichLength / fTcalEvent->geom_detector.fScintillatorVoxelSize;
    int nztot = fTcalEvent->geom_detector.NRep * nzlayer;
    
    const auto &clusters3D = *fPORecoEvent->GetPSClusters3D();

    std::cout << "3D histogram " << nx << " " << ny << " " << nztot << " " << std::endl;  

    TH3D *h3D = new TH3D("h3D", "3D Cluster Hits;X (cm);Y (cm);Z (cm)",
                         nx, 0, nx,
                         ny, 0, ny,
                         nztot, 0, nztot);

    std::map<int, int> clusterColorMap;
    int clusterIndex = 1;
    for (const auto &cluster : clusters3D)
      {
	if (clusterColorMap.find(cluster.clusterID)== clusterColorMap.end())
	  clusterColorMap[cluster.clusterID] = clusterIndex++;
      }
    for (const auto &cluster : clusters3D)
      {
	int colorIndex = clusterColorMap[cluster.clusterID];
        for (const auto &hit : cluster.hits)
	  {
	    double x, y, z;
	    // std::cout << "check 3D clusterIDs: " << cluster.clusterID << std::endl;
            fPORecoEvent->pshit2d_position(hit.id, x, y, z);
	    //h3D->Fill(x, y, z, cluster.clusterID);
	    h3D->Fill(x, y, z, colorIndex);
	  }
      }
    gStyle->SetPalette(kRainBow);
    h3D->GetXaxis()->SetTitleOffset(1.2);
    h3D->GetYaxis()->SetTitleOffset(1.4);
    h3D->Draw("");

    myCan->Update();
  }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  void FaserCalDisplay::CountHitsInCube()
  {
    std::map<std::tuple<int, int, int>, int> scintillatorHitCounts; // To count the hits on scintillator cubes
    std::map<std::tuple<int, int, int>, int> electronShowerHitCounts; // To count the hits on scintillator cubes for electron showers
    std::map<std::tuple<int, int, int>, int> hadronShowerHitCounts; // To count the hits on scintillator cubes for hadron showers
    double wdx, wdy, wdz;
    if (TGeoBBox *box = dynamic_cast<TGeoBBox*>(gGeoManager->GetTopVolume()->GetShape())){
      wdx = box->GetDX();
      wdy = box->GetDY();
      wdz = box->GetDZ()/2.0;
    } else {
      std::cerr << "Failed to get top volume dimensions." << std::endl;
      exit(1);
    }
    for (const auto& track : fTcalEvent->getfTracks()) {
      size_t nhits = track->fhitIDs.size();
      for (size_t i = 0; i < nhits; i++) {
	long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
	// Apply energy cut on scintillator voxel
	if (hittype == 0 && track->fEnergyDeposits[i] < 0.5) continue;
	ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
	int x = static_cast<int>(position.X());
	int y = static_cast<int>(position.Y());
	int z = static_cast<int>(position.Z());
	if (hittype == 0) {
	  
	  // Increment the hit count for the scintillator cube
	  scintillatorHitCounts[std::make_tuple(x, y, z)]++;
	  if (fabs(track->fPDG) == 11) {
	    electronShowerHitCounts[std::make_tuple(x, y, z)]++;
	  }else if (fabs(track->fPDG) != 11 && fabs(track->fPDG) != 13) {
	    hadronShowerHitCounts[std::make_tuple(x, y, z)]++;
	  }
	}
      }
    }
    SetMyStyle();
    TCanvas *myCan = CreateCanvas("HitsInCubes",1);
    gPad->Update();
    myCan->Divide(1,2);
    if (gROOT->FindObject("hScintillatorHits")) delete gROOT->FindObject("hScintillatorHits");
    if (gROOT->FindObject("hElectronShowerHits")) delete gROOT->FindObject("hElectronShowerHits");
    if (gROOT->FindObject("hHadronShowerHits")) delete gROOT->FindObject("hHadronShowerHits");
    if (gROOT->FindObject("hScintillatorHitsVsZ")) delete gROOT->FindObject("hScintillatorHitsVsZ");

    TH1F* hScintillatorHits = new TH1F("hScintillatorHits", " Hit Counts in Cubes;Hit Counts;Entries", 20, -0.5, 19.5);
    TH1F* hElectronShowerHits = new TH1F("hElectronShowerHits", "Electron Shower Hit Counts in Cubes;Hit Counts;Entries", 20, -0.5, 19.5);
    TH1F* hHadronShowerHits = new TH1F("hHadronShowerHits", "Hadron Shower Hit Counts in Cubes;Hit Counts;Entries", 20, -0.5, 19.5);
    TH2F* hScintillatorHitsVsZ = new TH2F("hScintillatorHitsVsZ", "Hit Counts vs Z Position;Z Position;Number of Hits", 2*wdz, -wdz, wdz, 20, -0.5, 19.5);
    // hit counts on scintillator cubes
    for (const auto& item : scintillatorHitCounts) {
      int x, y, z;
      std::tie(x, y, z) = item.first;
      int count = item.second;
      if(fVerbose)
	std::cout << "Scintillator cube at (" << x << ", " << y << ", " << z << ") has " << count << " hits." << std::endl;
      hScintillatorHits->Fill(count);
      hScintillatorHitsVsZ->Fill(z, count);
    }
    for (const auto& item : electronShowerHitCounts) {
      int count = item.second;
      hElectronShowerHits->Fill(count);
    }

    for (const auto& item : hadronShowerHitCounts) {
      int count = item.second;
      hHadronShowerHits->Fill(count);
    }

    myCan->cd(1);
    hScintillatorHits->Draw();
    hScintillatorHits->SetLineWidth(2);
    hElectronShowerHits->SetLineColor(kBlue);
    hElectronShowerHits->SetLineWidth(2);
    hElectronShowerHits->Draw("same");
    hHadronShowerHits->SetLineColor(kRed);
    hHadronShowerHits->SetLineWidth(2);
    hHadronShowerHits->Draw("same");
    TLegend *myl = new TLegend(0.6,0.5,0.8,0.7);
    myl->AddEntry(hScintillatorHits,"All","L");
    myl->AddEntry(hElectronShowerHits,"EM shower","L");
    myl->AddEntry(hHadronShowerHits,"Had shower","L");
    myl->SetBorderSize(0);
    myl->SetFillColor(0);
    myl->Draw("same");   
    gPad->SetLogy();
    gPad->Update();
    myCan->cd(2);
    hScintillatorHitsVsZ->SetStats(0);
    hScintillatorHitsVsZ->Draw("COLZ");
    gPad->Update();
    myCan->Update();
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

void FaserCalDisplay::PlotParticleFractionsAndHitsInCubes() {
    // Maps to count particle types and voxel hits
    std::map<int, int> particleTypeCounts; // Count of each particle type
    std::map<std::tuple<int, int, int>, int> voxelHitCounts; // Hits in each voxel

    // Loop over tracks and count particle types and voxel hits
    for (const auto& track : fTcalEvent->getfTracks()) {
        size_t nhits = track->fhitIDs.size();
        for (size_t i = 0; i < nhits; i++) {
            long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
            if (hittype == 0 && track->fEnergyDeposits[i] < 0.5) continue; // Energy cut for scintillator voxels

            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
            int x = static_cast<int>(position.X());
            int y = static_cast<int>(position.Y());
            int z = static_cast<int>(position.Z());

            // Increment voxel hit count
            voxelHitCounts[std::make_tuple(x, y, z)]++;

            // Increment particle type count
            particleTypeCounts[track->fPDG]++;
        }
    }

    // Create histograms for particle fractions
    TH1F* hParticleFractions = new TH1F("hParticleFractions", "Particle Type Fractions;Particle Type (PDG);Fraction", particleTypeCounts.size(), 0, particleTypeCounts.size());
    TH1F* hVoxelHits = new TH1F("hVoxelHits", "Voxel Hit Counts;Number of Hits;Entries", 20, -0.5, 19.5);

    // Fill particle fraction histogram
    int totalParticles = 0;
    for (const auto& [pdg, count] : particleTypeCounts) {
        totalParticles += count;
    }
    int bin = 1;
    for (const auto& [pdg, count] : particleTypeCounts) {
        double fraction = static_cast<double>(count) / totalParticles;
        hParticleFractions->GetXaxis()->SetBinLabel(bin, Form("%d", pdg));
        hParticleFractions->SetBinContent(bin, fraction);
        bin++;
    }

    // Fill voxel hit histogram
    for (const auto& [voxel, count] : voxelHitCounts) {
        hVoxelHits->Fill(count);
    }

    // Create canvas and draw histograms
    TCanvas* canvas = new TCanvas("canvas", "Particle Fractions and Voxel Hits", 1200, 600);
    canvas->Divide(2, 1);

    // Draw particle fractions
    canvas->cd(1);
    hParticleFractions->SetFillColor(kBlue - 9);
    hParticleFractions->SetBarWidth(0.8);
    hParticleFractions->SetBarOffset(0.1);
    hParticleFractions->Draw("bar");

    // Draw voxel hits
    canvas->cd(2);
    hVoxelHits->SetLineColor(kRed);
    hVoxelHits->SetLineWidth(2);
    hVoxelHits->Draw();

    // Update canvas
    canvas->Update();
}



  /*
  // Function to compute momentum from voxel position and total energy deposited
  fastjet::PseudoJet FaserCalDisplay::computeMomentumFromVoxel(ROOT::Math::XYZVector position, double totalEnergy) 
  {
    double norm = sqrt(position.x() * position.x() + position.y() * position.y() + position.z() * position.z());
    if (norm == 0) return fastjet::PseudoJet(0, 0, 0, 0); // Avoid division by zero
    double px = totalEnergy * position.x() / norm;
    double py = totalEnergy * position.y() / norm;
    double pz = totalEnergy * position.z() / norm;
    return fastjet::PseudoJet(px, py, pz, totalEnergy);
  }
  */
  /*
  void FaserCalDisplay::JetReconstructions()
  {  
    std::map<std::tuple<int, int, int>, double> voxelEnergyDeposits; // Store total energy per voxel
    
    // Aggregate energy per voxel
    for (const auto& track : fTcalEvent->getfTracks()) {
      for (size_t i = 0; i < track->fhitIDs.size(); i++) {
	long hittype = fTcalEvent->getChannelTypefromID(track->fhitIDs[i]);
        
	if (hittype == 0) { // Only consider scintillator hits
	  ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
	  int x = static_cast<int>(position.X());
	  int y = static_cast<int>(position.Y());
	  int z = static_cast<int>(position.Z());
	  
	  double energyDeposit = track->fEnergyDeposits[i];
	  
	  std::tuple<int, int, int> voxelKey = std::make_tuple(x, y, z);
	  voxelEnergyDeposits[voxelKey] += energyDeposit; // Sum energy deposits in the same voxel
	}
      }
    }
    // Convert voxel energy deposits into momentum vectors
    std::vector<fastjet::PseudoJet> input_particles;
    for (const auto& voxel : voxelEnergyDeposits) {
      int x, y, z;
      std::tie(x, y, z) = voxel.first;
      double totalEnergy = voxel.second;
      
      // Convert voxel (x, y, z) back to real space position
      ROOT::Math::XYZVector position(x, y, z);
      
      // Compute momentum using total voxel energy
      fastjet::PseudoJet p = computeMomentumFromVoxel(position, totalEnergy);
      input_particles.push_back(p);
    }
    // Jet clustering parameters
    double R = 0.4;
    double ptMin = 0.1;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    fastjet::ClusterSequence cs(input_particles, jet_def);
    
    // Store the jets above ptMin threshold
    fJets = cs.inclusive_jets(ptMin);
    
    // Debug output: Print jet properties
    for (const auto& jet : fJets) {
      std::cout << "Jet pt: " << jet.pt() 
		<< ", eta: " << jet.eta() 
		<< ", phi: " << jet.phi() 
		<< ", mass: " << jet.m() 
		<< ", n_constituents: " << jet.constituents().size()
		<< std::endl;
    }
    // Save jet data for further analysis
    std::ofstream jetFile("jets_with_voxel_energy.txt");
    if (!jetFile) {
      std::cerr << "Error opening file!" << std::endl;
      return;
    }
    
    std::map<int, TPolyMarker3D*> jetMarkers;
    std::vector<int> colors = {2, 4, 6, 8, 9, 11, 28, 46, 49}; // ROOT colors
    int colorIndex = 0;
    
    int jet_id = 0;
    for (const auto& jet : fJets) {
      jetFile << "JET " << jet_id << " " << jet.pt() << " " << jet.eta() << " " << jet.phi() << "\n";
      
      jetMarkers[jet_id] = new TPolyMarker3D();
      jetMarkers[jet_id]->SetMarkerStyle(20);
      jetMarkers[jet_id]->SetMarkerColor(colors[colorIndex % colors.size()]);
      colorIndex++;
      
      for (const auto& constituent : jet.constituents()) {
        jetFile << jet_id << " " << constituent.pt() << " " << constituent.eta() << " " << constituent.phi() << "\n";
        jetMarkers[jet_id]->SetNextPoint(constituent.eta(), constituent.phi(), constituent.pt());
      }
      jetFile << "END\n";
      jet_id++;
    }
    jetFile.close();
    
    // Create a 3D canvas for visualization
    TCanvas* c1 = new TCanvas("c1", "3D Jet Visualization", 800, 600);
    c1->cd();
    
    // Define 3D histogram frame
    TH3F* frame = new TH3F("frame", "Jets in (eta, phi, pT);#eta;#phi;pT",
			   10, -5, 5, 10, -3.14, 3.14, 10, 0, 15000);
    frame->Draw();
    
    // Draw jet constituents
    TLegend* legend = new TLegend(0.8, 0.7, 0.9, 0.9);
    for (const auto& [jet, marker] : jetMarkers) {
      marker->Draw("same");
      legend->AddEntry(marker, Form("Jet %d", jet), "p");
    }
    
    legend->Draw();
    c1->Update();
    
    ShowJetHits();
  }
  */
    //////////////////////////////////////////////////////////
    /*
    void FaserCalDisplay::JetReconstructions()
    {  
     std::vector<fastjet::PseudoJet> input_particles;
 
     for (const auto& track : fTcalEvent->getfTracks())
       {
         for ( size_t i = 0; i < track->fhitIDs.size(); i++)
           {
             ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(track->fhitIDs[i]);
             long module = fTcalEvent->getChannelModulefromID(track->fhitIDs[i]);
             double norm = TMath::Sqrt(position.x()*position.x() + position.y()*position.y() + position.z()*position.z());
             double p_x = position.x()/norm;
             double p_y = position.y()/norm;
             double p_z = position.z()/norm;
             double energy = track->fEnergyDeposits[i];
             // Create a PseudoJet from the track
             fastjet::PseudoJet p(p_x, p_y, p_z, energy);
             //p.set_user_index(ID);
             input_particles.push_back(p);
           }
         }
     // Define the jet algorithm and cluster sequence
     double R, ptMin;
     R = 0.4;
     ptMin = 0.1;
     fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
     fastjet::ClusterSequence cs(input_particles, jet_def);
 
     // Store the jets with pt above ptMin threshold
     fJets = cs.inclusive_jets(ptMin);
 
     // Debug output
     for (const auto& jet : fJets) {
         std::cout << "Jet pt: " << jet.pt() 
                   << ", eta: " << jet.eta() 
                   << ", phi: " << jet.phi() 
                   << ", mass: " << jet.m() 
                   << ", rap[y]: " << jet.rap() 
 
       << std::endl;
     }
     // Extract jets above ptMin
     std::vector<fastjet::PseudoJet> jets = cs.inclusive_jets(ptMin);
     
 
     std::ofstream jetFile("jets_with_hits.txt"); // Open file
     if (!jetFile) {
         std::cerr << "Error opening file!" << std::endl;
         return;
     }
 
     // Store jet data
     std::map<int, TPolyMarker3D*> jetMarkers;
     std::vector<int> colors = {2, 4, 6, 8, 9, 11, 28, 46, 49}; // ROOT predefined colors
     int colorIndex = 0;
 
     int jet_id = 0;
     for (const auto& jet : jets) { // Loop over jets
       jetFile << "JET " << jet_id << " " << jet.pt() << " " << jet.eta() << " " << jet.phi() << "\n";
       jetMarkers[jet_id] = new TPolyMarker3D();
       jetMarkers[jet_id]->SetMarkerStyle(20);
       jetMarkers[jet_id]->SetMarkerColor(colors[colorIndex % colors.size()]);
       colorIndex++;
       for (const auto& constituent : jet.constituents()) { // Loop over hits inside the jet
         jetFile << jet_id << " " << constituent.pt() << " " << constituent.eta() << " " << constituent.phi() << "\n";
     jetMarkers[jet_id]->SetNextPoint(constituent.eta(), constituent.phi(), constituent.pt());
     }
       jetFile << "END\n"; // Mark end of jet
       jet_id++;
     }
     jetFile.close();
 
     // Create canvas
     TCanvas* c1 = new TCanvas("c1", "3D Jet Visualization", 800, 600);
     c1->cd();
 
     // Create a 3D frame
     TH3F* frame = new TH3F("frame", "Jets in (eta, phi, pT);#eta;#phi;pT", 10, -5, 5, 10, -3.14, 3.14, 10, 0, 15000);
     frame->Draw();
 
     // Draw all jet points
     TLegend* legend = new TLegend(0.8, 0.7, 0.9, 0.9);
     for (const auto& [jet, marker] : jetMarkers) {
         marker->Draw("same");
         legend->AddEntry(marker, Form("Jet %d", jet), "p");
     }
 
     legend->Draw();
     c1->Update();
 
    }*/

  /*
    void FaserCalDisplay::ShowJetHits()
    {
        if (!gEve) gEve = TEveManager::Create();
    
        // Select the highest pT jet for visualization
        if (fJets.empty()) {
            std::cerr << "No jets found!" << std::endl;
            return;
        }
    
        // Sort jets by pT
        std::sort(fJets.begin(), fJets.end(), [](const fastjet::PseudoJet& a, const fastjet::PseudoJet& b) {
            return a.pt() > b.pt();
        });
    
        fastjet::PseudoJet selectedJet = fJets[0]; // Select the highest pT jet
    
        TEveLine *jetLine = new TEveLine();
        jetLine->SetName("Jet Trajectory");
        jetLine->SetMainColor(kRed);
        jetLine->SetLineWidth(3);
    
        TEvePointSet *jetHits = new TEvePointSet();
        jetHits->SetName("Jet Hits");
        jetHits->SetMarkerColor(kBlue);
        jetHits->SetMarkerStyle(20);
        jetHits->SetMarkerSize(1.0);
        std::cout << "here i am " << std::endl;
        // Iterate over the jet's constituents (hit voxels)
        for (const auto& constituent : selectedJet.constituents()) {
            long hitID = constituent.user_index(); // Assuming user_index stores hit ID
    
            // Retrieve hit position in detector
            ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(hitID);
            double x = position.X();
            double y = position.Y();
            double z = position.Z();
    
            // Add hit points
            jetHits->SetNextPoint(x, y, z);
    
            // Add to the jet trajectory
            jetLine->SetNextPoint(x, y, z);
        }
    
        // Add elements to Eve
        gEve->AddElement(jetLine);
        gEve->AddElement(jetHits);
    
        // Update the visualization
        gEve->FullRedraw3D(kTRUE);    

        //gEve->Redraw3D(kTRUE);
    }
  */

}  // namespace display
