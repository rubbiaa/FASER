#include "TGeoManager.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TGeoMatrix.h"

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TGeoManager.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <string>


void Geometry()
{

  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libEve");
  // --- Load the geometry ---
  TGeoManager::Import("geometry_tilted_5degree.gdml");
  
  // --- Create TEve Manager ---
  TEveManager::Create();
  
  // --- Import geometry into TEve ---
  TEveGeoTopNode* topNode =
    new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  topNode->SetVisLevel(4);  // choose detail level (2–5 works well)
  
  gEve->AddGlobalElement(topNode);
  
  // --- Show TEve GUI ---
  gEve->Redraw3D(kTRUE);
}


void Geometry_flux()
{

  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libEve");
  // --- Load the geometry ---
  TGeoManager::Import("geometry_tilted_5degree.gdml");
  
  // --- Create TEve Manager ---
  TEveManager::Create();
  
  // --- Import geometry into TEve ---
  TEveGeoTopNode* topNode =
    new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  topNode->SetVisLevel(4);  // choose detail level (2–5 works well)
  
  gEve->AddGlobalElement(topNode);

  
  // -------------------------------------------------------
  //  FLUX: read neutrino vertices + directions
  // -------------------------------------------------------
  TFile *fflux = TFile::Open("events_light_4x4.root");
  if (!fflux || fflux->IsZombie()) {
    Error("DrawGDML_TEve_WithFlux", "Cannot open flux.root");
    gEve->Redraw3D(kTRUE);
    return;
  }

  TTree *t = (TTree*) fflux->Get("flux");
  if (!t) {
    Error("DrawGDML_TEve_WithFlux", "Tree 'flux' not found in flux.root");
    gEve->Redraw3D(kTRUE);
    return;
  }
  
  Double_t vtxx, vtxy, vtxz;   // starting point
  Double_t px,  py,  pz;       // direction / momentum
  
  t->SetBranchAddress("vtxx", &vtxx);   // <-- change if needed
  t->SetBranchAddress("vtxy", &vtxy);   // <-- change if needed
  t->SetBranchAddress("vtxz", &vtxz);   // <-- change if needed
  t->SetBranchAddress("px",   &px);     // <-- change if needed
  t->SetBranchAddress("py",   &py);     // <-- change if needed
  t->SetBranchAddress("pz",   &pz);     // <-- change if needed
  
  Int_t nEntries  = t->GetEntries();
  Int_t maxToDraw = 10000;      // limit to avoid a crazy number of arrows
  Double_t L      = 10.0;     // arrow length in geometry units (e.g. cm)
  
  for (Int_t i = 0; i < nEntries && i < maxToDraw; ++i) {
    t->GetEntry(i);
    vtxx = 100*vtxx;
    vtxy = 100*vtxy;
    vtxz = 100*vtxz;
    
    // Normalize direction so arrows all have length L
    Double_t p = TMath::Sqrt(px*px + py*py + pz*pz);
    if (p == 0) continue;
    
    Double_t dx = (px / p) * L;
    Double_t dy = (py / p) * L;
    Double_t dz = (pz / p) * L;
    
    // Line from vertex to vertex + L * direction
    TEveLine *line = new TEveLine(2);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetNextPoint(vtxx,        vtxy,        vtxz       );
    line->SetNextPoint(vtxx + dx,   vtxy + dy,   vtxz + dz );
    
    gEve->AddElement(line);
  }
  
  
  
  // --- Show TEve GUI ---
  gEve->Redraw3D(kTRUE);
}


/////////////

void Geometry_flux_fasernu()
{

  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libEve");
  // --- Load the geometry ---
  TGeoManager::Import("FaserNu3.gdml");
  
  // --- Create TEve Manager ---
  TEveManager::Create();
  
  // --- Import geometry into TEve ---
  TEveGeoTopNode* topNode =
    new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
  topNode->SetVisLevel(4);  // choose detail level (2–5 works well)
  
  gEve->AddGlobalElement(topNode);  
  // -------------------------------------------------------
  //  FLUX: read neutrino vertices + directions
  // -------------------------------------------------------
  TFile *fflux = TFile::Open("events_light_4x4.root");
  if (!fflux || fflux->IsZombie()) {
    Error("DrawGDML_TEve_WithFlux", "Cannot open flux.root");
    gEve->Redraw3D(kTRUE);
    return;
  }

  TTree *t = (TTree*) fflux->Get("flux");
  if (!t) {
    Error("DrawGDML_TEve_WithFlux", "Tree 'flux' not found in flux.root");
    gEve->Redraw3D(kTRUE);
    return;
  }
  
  Double_t vtxx, vtxy, vtxz;   // starting point
  Double_t px,  py,  pz;       // direction / momentum
  
  t->SetBranchAddress("vtxx", &vtxx);   // <-- change if needed
  t->SetBranchAddress("vtxy", &vtxy);   // <-- change if needed
  t->SetBranchAddress("vtxz", &vtxz);   // <-- change if needed
  t->SetBranchAddress("px",   &px);     // <-- change if needed
  t->SetBranchAddress("py",   &py);     // <-- change if needed
  t->SetBranchAddress("pz",   &pz);     // <-- change if needed
  
  Int_t nEntries  = t->GetEntries();
  Int_t maxToDraw = 10000;      // limit to avoid a crazy number of arrows
  Double_t L      = 10.0;     // arrow length in geometry units (e.g. cm)
  
  for (Int_t i = 0; i < nEntries && i < maxToDraw; ++i) {
    t->GetEntry(i);
    vtxx = 100*vtxx;
    vtxy = 100*vtxy;
    vtxz = 100*vtxz;
    
    // Normalize direction so arrows all have length L
    Double_t p = TMath::Sqrt(px*px + py*py + pz*pz);
    if (p == 0) continue;
    
    Double_t dx = (px / p) * L;
    Double_t dy = (py / p) * L;
    Double_t dz = (pz / p) * L;
    
    // Line from vertex to vertex + L * direction
    TEveLine *line = new TEveLine(2);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetNextPoint(vtxx,        vtxy,        vtxz       );
    line->SetNextPoint(vtxx + dx,   vtxy + dy,   vtxz + dz );
    
    gEve->AddElement(line);
  }
  
  // --- Show TEve GUI ---
  gEve->Redraw3D(kTRUE);
}

/////////////////

// Helper function to recursively set transparency
void SetVolumeTransparencyRecursive(TGeoVolume* vol, Int_t transparency, Int_t maxDepth = 10, Int_t currentDepth = 0) {
  if (!vol || currentDepth > maxDepth) return;
  
  vol->SetTransparency(transparency);
  
  Int_t nDaughters = vol->GetNdaughters();
  for (Int_t i = 0; i < nDaughters; ++i) {
    TGeoNode *node = vol->GetNode(i);
    if (node) {
      TGeoVolume *daughterVol = node->GetVolume();
      SetVolumeTransparencyRecursive(daughterVol, transparency, maxDepth, currentDepth + 1);
    }
  }
}

void Geometry_flux_fasernu_fasercal()
{
// ----------------------------------------------------
  // 1. Load libraries
  // ----------------------------------------------------
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libEve");
  gSystem->Load("libGdml");   // for TGDMLParse

  // ----------------------------------------------------
  // 2. Import main geometry (FaserNu3) into ONE TGeoManager
  // ----------------------------------------------------
  TGeoManager::Import("FaserNu3.gdml");
  TGeoManager *geom = gGeoManager;
  if (!geom) {
    Error("Geometry_flux_fasernu", "Cannot load FaserNu3.gdml");
    return;
  }

  // ----------------------------------------------------
  // 3. Import second geometry as a volume via TGDMLParse
  //    and attach it under the existing top volume
  // ----------------------------------------------------
  TGDMLParse parser;
 //TGeoVolume *vol2 = parser.GDMLReadFile("/Users/ukose/sw/kits/NewFASER/FASER/GeomGDML/FASERCAL_V8_W1mm_7deg_LoSX_38cm_LoSY_17cm_3DCAL_0inXm5inY.gdml");
 TGeoVolume *vol2 = parser.GDMLReadFile("/Users/ukose/sw/kits/NewFASER/FASER/FASERG4/build/FASERCAL_V9.gdml");
  if (!vol2) {
    Error("Geometry_flux_fasernu", "Cannot load /Users/ukose/sw/kits/NewFASER/FASER/FASERG4/build/FASERCAL_V9.gdml as volume");
    return;
  }

  // Optionally rename to avoid name clashes
  vol2->SetName("CalGeometry");

  // Set position / rotation of the second geometry relative to FaserNu3
  // (currently at origin; change x,y,z and rotations as needed)
  // TGeoTranslation *tr = new TGeoTranslation("tr_cal", 11.0, 17.0, 820.0); // in cm
  //TGeoTranslation *tr = new TGeoTranslation("tr_cal", 20.0, 17.0, 0.0); // in cm
  TGeoTranslation *tr = new TGeoTranslation("tr", 0.0, .0, 600.0); // in cm
  tr->RegisterYourself();  // good practice if you later use in matrices

  geom->GetTopVolume()->AddNode(vol2, 1, tr);
  //geom->CloseGeometry();


  // ----------------------------------------------------
  // 3. Import second geometry as a volume via TGDMLParse
  //    and attach it under the existing top volume
  // ----------------------------------------------------
  //TGDMLParse parser2;
  //TGeoVolume *vol3 = parser2.GDMLReadFile("/Users/ukose/sw/kits/NewFASER/FASER/GeomGDML/test.gdml");
  //if (!vol3) {
  //  Error("Geometry_flux_fasernu", "Cannot load test.gdml as volume");
  //  return;
// }

  // Optionally rename to avoid name clashes
  //vol3->SetName("CalGeometry2");

  // Set position / rotation of the second geometry relative to FaserNu3
  // (currently at origin; change x,y,z and rotations as needed)
  // TGeoTranslation *tr = new TGeoTranslation("tr_cal", 11.0, 17.0, 820.0); // in cm
  //TGeoTranslation *tr3 = new TGeoTranslation("tr_cal", 0,0, 600.0); // in cm
  TGeoTranslation *tr3 = new TGeoTranslation("tr_cal", 0.0, 0.0, 600.0); // in cm
  //tr3->RegisterYourself();  // good practice if you later use in matrices

  //geom->GetTopVolume()->AddNode(vol3, 1, tr3);
  geom->CloseGeometry();

  // ----------------------------------------------------
  // 3b. EXPORT combined geometry to GDML file
  // ----------------------------------------------------
  std::cout << "\n=== Exporting combined geometry ===" << std::endl;
  geom->Export("combined_fasernu_fasercal.gdml");
  geom->Export("combined_fasernu_fasercal.root");
  std::cout << "Exported to: combined_fasernu_fasercal.gdml" << std::endl;
  std::cout << "You can also export to ROOT format: combined_fasernu_fasercal.root" << std::endl;
  std::cout << "===================================\n" << std::endl;

  // ----------------------------------------------------
  // 4. Start TEve and show the combined geometry
  // ----------------------------------------------------
  TEveManager::Create();

  TEveGeoTopNode *topNode =
      new TEveGeoTopNode(geom, geom->GetTopNode());
  topNode->SetVisLevel(4);
  gEve->AddGlobalElement(topNode);

  // ----------------------------------------------------
  // 4b. Find 3DCAL module and draw a line through it
  // ----------------------------------------------------
  TGeoNode *calNode = nullptr;
  TGeoVolume *topVol = geom->GetTopVolume();
  Int_t nNodes = topVol->GetNdaughters();
  
  // Find CalGeometry node
  for (Int_t i = 0; i < nNodes; ++i) {
    TGeoNode *node = topVol->GetNode(i);
    if (TString(node->GetName()).Contains("CalGeometry")) {
      calNode = node;
      break;
    }
  }

  if (calNode) {
    TGeoVolume *calVol = calNode->GetVolume();
    TGeoMatrix *calMatrix = calNode->GetMatrix();
    const Double_t *calTrans = calMatrix->GetTranslation();
    
    // First level: Look for DetectorAssemblyPV
    Int_t nDaughters = calVol->GetNdaughters();
    std::cout << "\n=== FASERCAL First Level ===" << std::endl;
    std::cout << "Total nodes: " << nDaughters << std::endl;
    
    TGeoNode *assemblyNode = nullptr;
    if (nDaughters == 1) {
      assemblyNode = calVol->GetNode(0);
      std::cout << "Found: " << assemblyNode->GetName() << std::endl;
    }
    
    if (!assemblyNode) {
      std::cout << "ERROR: DetectorAssembly node not found" << std::endl;
      return;
    }
    
    // Navigate into DetectorAssemblyPV to find the actual calorimeters
    TGeoVolume *assemblyVol = assemblyNode->GetVolume();
    TGeoMatrix *assemblyMatrix = assemblyNode->GetMatrix();
    const Double_t *assemblyTrans = assemblyMatrix->GetTranslation();
    
    // Set theyarency recursively for all calorimeter modules
    Int_t nCalModules = assemblyVol->GetNdaughters();
    std::cout << "\nSetting transparency for " << nCalModules << " calorimeter modules..." << std::endl;
    for (Int_t i = 0; i < nCalModules; ++i) {
      TGeoNode *node = assemblyVol->GetNode(i);
      TGeoVolume *vol = node->GetVolume();
      SetVolumeTransparencyRecursive(vol, 30); // Recursively set transparency
      std::cout << "  Set transparency for: " << node->GetName() << std::endl;
    }
    std::cout << "Transparency set to 60% for all volumes\n" << std::endl;
    
    // Extract rotation/tilt from the assembly matrix
    const Double_t *rotMatrix = assemblyMatrix->GetRotationMatrix();
    Double_t tiltAngle = 0.0;
    if (rotMatrix) {
      // For rotation around Y-axis: atan2(-rotMatrix[2], rotMatrix[8])
      // rotMatrix is row-major: [0,1,2; 3,4,5; 6,7,8]
      // For Y-rotation: R[0,0]=cos(θ), R[0,2]=sin(θ), R[2,0]=-sin(θ), R[2,2]=cos(θ)
      tiltAngle = TMath::ATan2(rotMatrix[2], rotMatrix[0]); // in radians
      std::cout << "\nDetectorAssembly tilt angle: " << tiltAngle * TMath::RadToDeg() << " degrees" << std::endl;
    }
    
    // Calculate tilted direction vector (rotation around Y-axis)
    Double_t dir_x = TMath::Sin(tiltAngle);
    Double_t dir_y = 0.0;
    Double_t dir_z = TMath::Cos(tiltAngle);
    std::cout << "Tilted direction vector: (" << dir_x << ", " << dir_y << ", " << dir_z << ")\n" << std::endl;
    
    std::cout << "\n=== Calorimeter Modules (inside DetectorAssembly) ===" << std::endl;
    std::cout << "Total modules: " << nCalModules << "\n" << std::endl;
    
    for (Int_t i = 0; i < nCalModules; ++i) {
      TGeoNode *node = assemblyVol->GetNode(i);
      TGeoMatrix *matrix = node->GetMatrix();
      const Double_t *localTrans = matrix->GetTranslation();
      
      // Calculate global position (local + assembly offset + CalGeometry offset)
      Double_t x = localTrans[0] + assemblyTrans[0] + calTrans[0];
      Double_t y = localTrans[1] + assemblyTrans[1] + calTrans[1];
      Double_t z = localTrans[2] + assemblyTrans[2] + calTrans[2];
      Double_t dist = TMath::Sqrt(x*x + y*y);
      
      TString nodeName(node->GetName());
      std::cout << std::setw(3) << i << ": " 
                << std::setw(50) << std::left << nodeName
                << " Pos: (" << std::setw(8) << std::right << std::fixed << std::setprecision(1) << x 
                << ", " << std::setw(8) << y 
                << ", " << std::setw(8) << z << ")"
                << "  Dist: " << std::setw(8) << dist << " cm";
      
      // Flag special detectors
      if (nodeName.Contains("Container") || nodeName.Contains("3D")) std::cout << " [3DCAL?]";
      if (nodeName.Contains("rearCal") && !nodeName.Contains("Hcal")) std::cout << " [ECAL?]";
      if (nodeName.Contains("rearHcal") || nodeName.Contains("Hcal")) std::cout << " [AHCAL?]";
      if (nodeName.Contains("Muon") || nodeName.Contains("MDT") || nodeName.Contains("Spect")) std::cout << " [MUON?]";
      std::cout << std::endl;
    }
    std::cout << "\n=== End of module list ===" << std::endl;
    
    // Search for calorimeter modules and draw lines through them
    Bool_t found3DCAL = kFALSE;
    Bool_t foundECAL = kFALSE;
    Bool_t foundAHCAL = kFALSE;
    Bool_t foundMuonSpect = kFALSE;
    
    Double_t lineLength = 2000.0; // cm
    
    std::cout << "\n=== COORDINATE SYSTEM EXPLANATION ===" << std::endl;
    std::cout << "Line of Sight (LoS): Z-axis from origin (0,0,0) in direction (0,0,1)" << std::endl;
    std::cout << "  - LoS passes through point (X=0, Y=0) for all Z values" << std::endl;
    std::cout << "  - Blue line in viewer represents the LoS" << std::endl;
    std::cout << "\nCenter Position (X, Y, Z):" << std::endl;
    std::cout << "  - Global 3D coordinates of detector center" << std::endl;
    std::cout << "  - X, Y = perpendicular offsets from LoS (dX, dY)" << std::endl;
    std::cout << "  - Z = position along beam axis" << std::endl;
    std::cout << "\nDistance from LoS = sqrt(X² + Y²) = sqrt(dX² + dY²)" << std::endl;
    std::cout << "  - Perpendicular distance from detector to LoS" << std::endl;
    std::cout << "\nTo calculate dX, dY at ANY Z position along tilted detector:" << std::endl;
    std::cout << "  Given: center (X_c, Y_c, Z_c), tilt direction (dir_x, dir_y, dir_z)" << std::endl;
    std::cout << "  At arbitrary Z_target:" << std::endl;
    std::cout << "    t = (Z_target - Z_c) / dir_z" << std::endl;
    std::cout << "    X(Z_target) = X_c + t * dir_x" << std::endl;
    std::cout << "    Y(Z_target) = Y_c + t * dir_y" << std::endl;
    std::cout << "  Then: dX = X(Z_target), dY = Y(Z_target)" << std::endl;
    std::cout << "======================================\n" << std::endl;
    
    std::cout << "\n=== Searching for Calorimeter Modules ===" << std::endl;
    
    for (Int_t i = 0; i < nCalModules; ++i) {
      TGeoNode *node = assemblyVol->GetNode(i);
      TString nodeName(node->GetName());
      TGeoMatrix *matrix = node->GetMatrix();
      const Double_t *localTrans = matrix->GetTranslation();
      
      // Calculate global position (local + assembly offset + CalGeometry offset)
      Double_t x = localTrans[0] + assemblyTrans[0] + calTrans[0];
      Double_t y = localTrans[1] + assemblyTrans[1] + calTrans[1];
      Double_t z = localTrans[2] + assemblyTrans[2] + calTrans[2];
      Double_t dist = TMath::Sqrt(x*x + y*y);
      
      // Search for 3DCAL (ContainerPlacement or ContainerBox)
      if (!found3DCAL && (nodeName.Contains("ContainerPlacement") || nodeName.Contains("ContainerBox"))) {
        found3DCAL = kTRUE;
        
        // Get volume dimensions
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *bbox = (TGeoBBox*)vol->GetShape();
        Double_t dz = bbox->GetDZ(); // half-length in z
        
        // Get combined transformation for this specific detector
        TGeoHMatrix combinedMat;
        combinedMat = *assemblyMatrix;
        combinedMat.Multiply(matrix);
        
        // Apply CalGeometry offset to combined matrix
        TGeoHMatrix finalMat;
        TGeoTranslation calTransMat(calTrans[0], calTrans[1], calTrans[2]);
        finalMat = calTransMat;
        finalMat.Multiply(&combinedMat);
        
        // Get actual global position from combined transformation
        const Double_t *globalTrans = finalMat.GetTranslation();
        x = globalTrans[0];
        y = globalTrans[1];
        z = globalTrans[2];
        dist = TMath::Sqrt(x*x + y*y);
        
        const Double_t *detRotMatrix = finalMat.GetRotationMatrix();
        
        // Calculate direction for this detector
        Double_t det_dir_x = TMath::Sin(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        Double_t det_dir_y = 0.0;
        Double_t det_dir_z = TMath::Cos(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        
        // Calculate start and end positions along tilted axis
        Double_t z_start = z - dz * det_dir_z;
        Double_t z_end = z + dz * det_dir_z;
        Double_t x_start = x - dz * det_dir_x;
        Double_t x_end = x + dz * det_dir_x;
        Double_t y_start = y - dz * det_dir_y;
        Double_t y_end = y + dz * det_dir_y;
        
        Double_t dist_start = TMath::Sqrt(x_start*x_start + y_start*y_start);
        Double_t dist_end = TMath::Sqrt(x_end*x_end + y_end*y_end);
        
        std::cout << "\n═══ 3DCAL ═══" << std::endl;
        std::cout << "  Name: " << node->GetName() << std::endl;
        std::cout << "\n  CENTER of detector:" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x << ", " << y << ", " << z << ") cm" << std::endl;
        std::cout << "    dX from LoS: " << x << " cm" << std::endl;
        std::cout << "    dY from LoS: " << y << " cm" << std::endl;
        std::cout << "    Perpendicular distance from LoS: " << dist << " cm = sqrt(" << x << "² + " << y << "²)" << std::endl;
        std::cout << "    Half-length along beam: " << dz << " cm" << std::endl;
        std::cout << "\n  UPSTREAM FACE (beam entrance):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_start << ", " << y_start << ", " << z_start << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_start << ": " << x_start << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_start << ": " << y_start << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_start << " cm" << std::endl;
        std::cout << "\n  DOWNSTREAM FACE (beam exit):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_end << ", " << y_end << ", " << z_end << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_end << ": " << x_end << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_end << ": " << y_end << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_end << " cm" << std::endl;
        
        // Draw green dashed line through 3DCAL with tilt (detector axis)
        TEveLine* cal3dLine = new TEveLine("3DCALLine");
        cal3dLine->SetLineColor(kGreen);
        cal3dLine->SetLineWidth(3);
        cal3dLine->SetLineStyle(2);
        // Line extends in tilted direction specific to this detector
        cal3dLine->SetNextPoint(x - det_dir_x * lineLength/2, y - det_dir_y * lineLength/2, z - det_dir_z * lineLength/2);
        cal3dLine->SetNextPoint(x + det_dir_x * lineLength/2, y + det_dir_y * lineLength/2, z + det_dir_z * lineLength/2);
        gEve->AddElement(cal3dLine);
        
        // Draw lines from LoS to upstream and downstream faces
        TEveLine* upstreamLine = new TEveLine("3DCAL_UpstreamToLoS");
        upstreamLine->SetLineColor(kGreen+2);
        upstreamLine->SetLineWidth(2);
        upstreamLine->SetLineStyle(1);
        upstreamLine->SetNextPoint(0, 0, z_start);  // Point on LoS at same Z
        upstreamLine->SetNextPoint(x_start, y_start, z_start);  // Point on upstream face
        gEve->AddElement(upstreamLine);
        
        TEveLine* downstreamLine = new TEveLine("3DCAL_DownstreamToLoS");
        downstreamLine->SetLineColor(kGreen-3);
        downstreamLine->SetLineWidth(2);
        downstreamLine->SetLineStyle(1);
        downstreamLine->SetNextPoint(0, 0, z_end);  // Point on LoS at same Z
        downstreamLine->SetNextPoint(x_end, y_end, z_end);  // Point on downstream face
        gEve->AddElement(downstreamLine);
        
        // Add markers at face positions
        TEvePointSet* upstreamPoint = new TEvePointSet("3DCAL_UpstreamFace");
        upstreamPoint->SetNextPoint(x_start, y_start, z_start);
        upstreamPoint->SetMarkerStyle(20);
        upstreamPoint->SetMarkerSize(1.5);
        upstreamPoint->SetMarkerColor(kGreen+2);
        gEve->AddElement(upstreamPoint);
        
        TEvePointSet* downstreamPoint = new TEvePointSet("3DCAL_DownstreamFace");
        downstreamPoint->SetNextPoint(x_end, y_end, z_end);
        downstreamPoint->SetMarkerStyle(20);
        downstreamPoint->SetMarkerSize(1.5);
        downstreamPoint->SetMarkerColor(kGreen-3);
        gEve->AddElement(downstreamPoint);
      }
      
      // Search for AHCAL first (rearHCal) - must come before ECAL check!
      if (!foundAHCAL && (nodeName.Contains("rearHCal") || nodeName.Contains("ContainerHcal"))) {
        foundAHCAL = kTRUE;
        
        // Get volume dimensions
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *bbox = (TGeoBBox*)vol->GetShape();
        Double_t dz = bbox->GetDZ(); // half-length in z
        
        // Get combined transformation for this specific detector
        TGeoHMatrix combinedMat;
        combinedMat = *assemblyMatrix;
        combinedMat.Multiply(matrix);
        
        // Apply CalGeometry offset to combined matrix
        TGeoHMatrix finalMat;
        TGeoTranslation calTransMat(calTrans[0], calTrans[1], calTrans[2]);
        finalMat = calTransMat;
        finalMat.Multiply(&combinedMat);
        
        // Get actual global position from combined transformation
        const Double_t *globalTrans = finalMat.GetTranslation();
        x = globalTrans[0];
        y = globalTrans[1];
        z = globalTrans[2];
        dist = TMath::Sqrt(x*x + y*y);
        
        const Double_t *detRotMatrix = finalMat.GetRotationMatrix();
        
        // Calculate direction for this detector
        Double_t det_dir_x = TMath::Sin(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        Double_t det_dir_y = 0.0;
        Double_t det_dir_z = TMath::Cos(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        
        // Calculate start and end positions along tilted axis
        Double_t z_start = z - dz * det_dir_z;
        Double_t z_end = z + dz * det_dir_z;
        Double_t x_start = x - dz * det_dir_x;
        Double_t x_end = x + dz * det_dir_x;
        Double_t y_start = y - dz * det_dir_y;
        Double_t y_end = y + dz * det_dir_y;
        
        Double_t dist_start = TMath::Sqrt(x_start*x_start + y_start*y_start);
        Double_t dist_end = TMath::Sqrt(x_end*x_end + y_end*y_end);
        
        std::cout << "\n═══ AHCAL ═══" << std::endl;
        std::cout << "  Name: " << node->GetName() << std::endl;
        std::cout << "\n  CENTER of detector:" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x << ", " << y << ", " << z << ") cm" << std::endl;
        std::cout << "    dX from LoS: " << x << " cm" << std::endl;
        std::cout << "    dY from LoS: " << y << " cm" << std::endl;
        std::cout << "    Perpendicular distance from LoS: " << dist << " cm = sqrt(" << x << "² + " << y << "²)" << std::endl;
        std::cout << "    Half-length along beam: " << dz << " cm" << std::endl;
        std::cout << "\n  UPSTREAM FACE (beam entrance):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_start << ", " << y_start << ", " << z_start << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_start << ": " << x_start << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_start << ": " << y_start << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_start << " cm" << std::endl;
        std::cout << "\n  DOWNSTREAM FACE (beam exit):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_end << ", " << y_end << ", " << z_end << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_end << ": " << x_end << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_end << ": " << y_end << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_end << " cm" << std::endl;
        
        // Draw cyan dashed line through AHCAL with tilt (detector axis)
        TEveLine* hcalLine = new TEveLine("AHCALLine");
        hcalLine->SetLineColor(kCyan);
        hcalLine->SetLineWidth(3);
        hcalLine->SetLineStyle(2);
        // Line extends in tilted direction specific to this detector
        hcalLine->SetNextPoint(x - det_dir_x * lineLength/2, y - det_dir_y * lineLength/2, z - det_dir_z * lineLength/2);
        hcalLine->SetNextPoint(x + det_dir_x * lineLength/2, y + det_dir_y * lineLength/2, z + det_dir_z * lineLength/2);
        gEve->AddElement(hcalLine);
        
        // Draw lines from LoS to upstream and downstream faces
        TEveLine* upstreamLine = new TEveLine("AHCAL_UpstreamToLoS");
        upstreamLine->SetLineColor(kCyan+2);
        upstreamLine->SetLineWidth(2);
        upstreamLine->SetLineStyle(1);
        upstreamLine->SetNextPoint(0, 0, z_start);
        upstreamLine->SetNextPoint(x_start, y_start, z_start);
        gEve->AddElement(upstreamLine);
        
        TEveLine* downstreamLine = new TEveLine("AHCAL_DownstreamToLoS");
        downstreamLine->SetLineColor(kCyan-3);
        downstreamLine->SetLineWidth(2);
        downstreamLine->SetLineStyle(1);
        downstreamLine->SetNextPoint(0, 0, z_end);
        downstreamLine->SetNextPoint(x_end, y_end, z_end);
        gEve->AddElement(downstreamLine);
        
        // Add markers at face positions
        TEvePointSet* upstreamPoint = new TEvePointSet("AHCAL_UpstreamFace");
        upstreamPoint->SetNextPoint(x_start, y_start, z_start);
        upstreamPoint->SetMarkerStyle(20);
        upstreamPoint->SetMarkerSize(1.5);
        upstreamPoint->SetMarkerColor(kCyan+2);
        gEve->AddElement(upstreamPoint);
        
        TEvePointSet* downstreamPoint = new TEvePointSet("AHCAL_DownstreamFace");
        downstreamPoint->SetNextPoint(x_end, y_end, z_end);
        downstreamPoint->SetMarkerStyle(20);
        downstreamPoint->SetMarkerSize(1.5);
        downstreamPoint->SetMarkerColor(kCyan-3);
        gEve->AddElement(downstreamPoint);
      }
      
      // Search for ECAL (rearCal but NOT rearHCal)
      else if (!foundECAL && (nodeName.Contains("rearCal") || nodeName.Contains("ContainerEcal"))) {
        foundECAL = kTRUE;
        
        // Get volume dimensions
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *bbox = (TGeoBBox*)vol->GetShape();
        Double_t dz = bbox->GetDZ(); // half-length in z
        
        // Get combined transformation for this specific detector
        TGeoHMatrix combinedMat;
        combinedMat = *assemblyMatrix;
        combinedMat.Multiply(matrix);
        
        // Apply CalGeometry offset to combined matrix
        TGeoHMatrix finalMat;
        TGeoTranslation calTransMat(calTrans[0], calTrans[1], calTrans[2]);
        finalMat = calTransMat;
        finalMat.Multiply(&combinedMat);
        
        // Get actual global position from combined transformation
        const Double_t *globalTrans = finalMat.GetTranslation();
        x = globalTrans[0];
        y = globalTrans[1];
        z = globalTrans[2];
        dist = TMath::Sqrt(x*x + y*y);
        
        const Double_t *detRotMatrix = finalMat.GetRotationMatrix();
        
        // Calculate direction for this detector
        Double_t det_dir_x = TMath::Sin(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        Double_t det_dir_y = 0.0;
        Double_t det_dir_z = TMath::Cos(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        
        // Calculate start and end positions along tilted axis
        Double_t z_start = z - dz * det_dir_z;
        Double_t z_end = z + dz * det_dir_z;
        Double_t x_start = x - dz * det_dir_x;
        Double_t x_end = x + dz * det_dir_x;
        Double_t y_start = y - dz * det_dir_y;
        Double_t y_end = y + dz * det_dir_y;
        
        Double_t dist_start = TMath::Sqrt(x_start*x_start + y_start*y_start);
        Double_t dist_end = TMath::Sqrt(x_end*x_end + y_end*y_end);
        
        std::cout << "\n═══ ECAL ═══" << std::endl;
        std::cout << "  Name: " << node->GetName() << std::endl;
        std::cout << "\n  CENTER of detector:" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x << ", " << y << ", " << z << ") cm" << std::endl;
        std::cout << "    dX from LoS: " << x << " cm" << std::endl;
        std::cout << "    dY from LoS: " << y << " cm" << std::endl;
        std::cout << "    Perpendicular distance from LoS: " << dist << " cm = sqrt(" << x << "² + " << y << "²)" << std::endl;
        std::cout << "    Half-length along beam: " << dz << " cm" << std::endl;
        std::cout << "\n  UPSTREAM FACE (beam entrance):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_start << ", " << y_start << ", " << z_start << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_start << ": " << x_start << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_start << ": " << y_start << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_start << " cm" << std::endl;
        std::cout << "\n  DOWNSTREAM FACE (beam exit):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_end << ", " << y_end << ", " << z_end << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_end << ": " << x_end << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_end << ": " << y_end << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_end << " cm" << std::endl;
        
        // Draw magenta dashed line through ECAL with tilt (detector axis)
        TEveLine* ecalLine = new TEveLine("ECALLine");
        ecalLine->SetLineColor(kMagenta);
        ecalLine->SetLineWidth(3);
        ecalLine->SetLineStyle(2);
        // Line extends in tilted direction specific to this detector
        ecalLine->SetNextPoint(x - det_dir_x * lineLength/2, y - det_dir_y * lineLength/2, z - det_dir_z * lineLength/2);
        ecalLine->SetNextPoint(x + det_dir_x * lineLength/2, y + det_dir_y * lineLength/2, z + det_dir_z * lineLength/2);
        gEve->AddElement(ecalLine);
        
        // Draw lines from LoS to upstream and downstream faces
        TEveLine* upstreamLine = new TEveLine("ECAL_UpstreamToLoS");
        upstreamLine->SetLineColor(kMagenta+2);
        upstreamLine->SetLineWidth(2);
        upstreamLine->SetLineStyle(1);
        upstreamLine->SetNextPoint(0, 0, z_start);
        upstreamLine->SetNextPoint(x_start, y_start, z_start);
        gEve->AddElement(upstreamLine);
        
        TEveLine* downstreamLine = new TEveLine("ECAL_DownstreamToLoS");
        downstreamLine->SetLineColor(kMagenta-3);
        downstreamLine->SetLineWidth(2);
        downstreamLine->SetLineStyle(1);
        downstreamLine->SetNextPoint(0, 0, z_end);
        downstreamLine->SetNextPoint(x_end, y_end, z_end);
        gEve->AddElement(downstreamLine);
        
        // Add markers at face positions
        TEvePointSet* upstreamPoint = new TEvePointSet("ECAL_UpstreamFace");
        upstreamPoint->SetNextPoint(x_start, y_start, z_start);
        upstreamPoint->SetMarkerStyle(20);
        upstreamPoint->SetMarkerSize(1.5);
        upstreamPoint->SetMarkerColor(kMagenta+2);
        gEve->AddElement(upstreamPoint);
        
        TEvePointSet* downstreamPoint = new TEvePointSet("ECAL_DownstreamFace");
        downstreamPoint->SetNextPoint(x_end, y_end, z_end);
        downstreamPoint->SetMarkerStyle(20);
        downstreamPoint->SetMarkerSize(1.5);
        downstreamPoint->SetMarkerColor(kMagenta-3);
        gEve->AddElement(downstreamPoint);
      }
      
      // Search for MuonSpectrometer in DetectorAssembly
      if (!foundMuonSpect && (nodeName.Contains("MuonSpectrometer") || nodeName.Contains("Muon") || 
                              nodeName.Contains("muon") || nodeName.Contains("MDT"))) {
        foundMuonSpect = kTRUE;
        
        std::cout << "\n>>> Found MuonSpectrometer: " << nodeName << " <<<" << std::endl;
        
        // Get volume dimensions
        TGeoVolume *vol = node->GetVolume();
        TGeoBBox *bbox = (TGeoBBox*)vol->GetShape();
        Double_t dx = bbox->GetDX(); // half-length in x
        Double_t dy = bbox->GetDY(); // half-length in y
        Double_t dz = bbox->GetDZ(); // half-length in z
        
        // Get combined transformation for this specific detector
        TGeoHMatrix combinedMat;
        combinedMat = *assemblyMatrix;
        combinedMat.Multiply(matrix);
        
        // Apply CalGeometry offset to combined matrix
        TGeoHMatrix finalMat;
        TGeoTranslation calTransMat(calTrans[0], calTrans[1], calTrans[2]);
        finalMat = calTransMat;
        finalMat.Multiply(&combinedMat);
        
        // Get actual global position from combined transformation
        const Double_t *globalTrans = finalMat.GetTranslation();
        x = globalTrans[0];
        y = globalTrans[1];
        z = globalTrans[2];
        dist = TMath::Sqrt(x*x + y*y);
        
        const Double_t *detRotMatrix = finalMat.GetRotationMatrix();
        
        // Calculate direction for this detector
        Double_t det_dir_x = TMath::Sin(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        Double_t det_dir_y = 0.0;
        Double_t det_dir_z = TMath::Cos(TMath::ATan2(detRotMatrix[2], detRotMatrix[0]));
        
        // Calculate start and end positions along detector axis
        Double_t z_start = z - dz * det_dir_z;
        Double_t z_end = z + dz * det_dir_z;
        Double_t x_start = x - dz * det_dir_x;
        Double_t x_end = x + dz * det_dir_x;
        Double_t y_start = y - dz * det_dir_y;
        Double_t y_end = y + dz * det_dir_y;
        
        Double_t dist_start = TMath::Sqrt(x_start*x_start + y_start*y_start);
        Double_t dist_end = TMath::Sqrt(x_end*x_end + y_end*y_end);
        
        std::cout << "\n═══ MUON SPECTROMETER ═══" << std::endl;
        std::cout << "  Name: " << node->GetName() << std::endl;
        std::cout << "\n  CENTER of detector:" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x << ", " << y << ", " << z << ") cm" << std::endl;
        std::cout << "    dX from LoS: " << x << " cm" << std::endl;
        std::cout << "    dY from LoS: " << y << " cm" << std::endl;
        std::cout << "    Perpendicular distance from LoS: " << dist << " cm = sqrt(" << x << "² + " << y << "²)" << std::endl;
        std::cout << "    Dimensions (X,Y,Z): (" << 2*dx << ", " << 2*dy << ", " << 2*dz << ") cm" << std::endl;
        std::cout << "    Half-length along beam: " << dz << " cm" << std::endl;
        std::cout << "\n  UPSTREAM FACE (beam entrance):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_start << ", " << y_start << ", " << z_start << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_start << ": " << x_start << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_start << ": " << y_start << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_start << " cm" << std::endl;
        std::cout << "\n  DOWNSTREAM FACE (beam exit):" << std::endl;
        std::cout << "    Position (X,Y,Z): (" << x_end << ", " << y_end << ", " << z_end << ") cm" << std::endl;
        std::cout << "    dX from LoS at Z=" << z_end << ": " << x_end << " cm" << std::endl;
        std::cout << "    dY from LoS at Z=" << z_end << ": " << y_end << " cm" << std::endl;
        std::cout << "    Perpendicular distance: " << dist_end << " cm" << std::endl;
        
        // List all tracking layers inside MuonSpectrometer
        Int_t nMuonLayers = vol->GetNdaughters();
        std::cout << "\n  TRACKING LAYERS inside MuonSpectrometer: " << nMuonLayers << std::endl;
        std::cout << "  " << std::string(90, '-') << std::endl;
        
        Int_t nMagnets = 0;
        Int_t nEndplates = 0;
        Int_t nSciFi = 0;
        
        // Track first and last MDTEndplate
        TGeoNode *firstEndplate = nullptr;
        TGeoNode *lastEndplate = nullptr;
        TGeoHMatrix firstEndplateMat;
        TGeoHMatrix lastEndplateMat;
        
        for (Int_t j = 0; j < nMuonLayers; ++j) {
          TGeoNode *layerNode = vol->GetNode(j);
          TString layerName(layerNode->GetName());
          TGeoMatrix *layerMatrix = layerNode->GetMatrix();
          const Double_t *layerTrans = layerMatrix->GetTranslation();
          
          // Transform layer position to global coordinates
          TGeoHMatrix layerFinalMat = finalMat;
          layerFinalMat.Multiply(layerMatrix);
          const Double_t *layerGlobalTrans = layerFinalMat.GetTranslation();
          
          Double_t layer_x = layerGlobalTrans[0];
          Double_t layer_y = layerGlobalTrans[1];
          Double_t layer_z = layerGlobalTrans[2];
          Double_t layer_dist = TMath::Sqrt(layer_x*layer_x + layer_y*layer_y);
          
          // Count component types and track endplates
          if (layerName.Contains("MDTMagnet")) nMagnets++;
          else if (layerName.Contains("MDTEndplate")) {
            nEndplates++;
            if (!firstEndplate) {
              firstEndplate = layerNode;
              firstEndplateMat = layerFinalMat;
            }
            lastEndplate = layerNode;
            lastEndplateMat = layerFinalMat;
          }
          else if (layerName.Contains("SciFi")) nSciFi++;
          
          // Mark special components
          TString marker = "  ";
          if (layerName.Contains("MDTMagnet")) marker = "🧲";
          else if (layerName.Contains("MDTEndplate")) marker = "📍";
          else if (layerName.Contains("SciFi")) marker = "🔬";
          
          std::cout << "    " << marker << " " << std::setw(3) << j << ": " 
                    << std::setw(35) << std::left << layerNode->GetName()
                    << " Z=" << std::setw(10) << std::right << std::fixed << std::setprecision(2) << layer_z
                    << " cm, Dist=" << std::setw(8) << layer_dist << " cm" << std::endl;
        }
        std::cout << "  " << std::string(90, '-') << std::endl;
        std::cout << "  Summary: " << nMagnets << " MDTMagnets, " 
                  << nEndplates << " MDTEndplates, " 
                  << nSciFi << " SciFi layers" << std::endl;
        std::cout << "  (🧲=Magnet, 📍=Endplate, 🔬=SciFi)" << std::endl;
        
        // Report detailed info on first and last MDTEndplate
        if (firstEndplate) {
          std::cout << "\n  ╔═══ FIRST MDT ENDPLATE ═══╗" << std::endl;
          const Double_t *fTrans = firstEndplateMat.GetTranslation();
          Double_t fx = fTrans[0], fy = fTrans[1], fz = fTrans[2];
          Double_t fdist = TMath::Sqrt(fx*fx + fy*fy);
          
          TGeoVolume *fVol = firstEndplate->GetVolume();
          TGeoBBox *fBox = (TGeoBBox*)fVol->GetShape();
          
          std::cout << "    Name: " << firstEndplate->GetName() << std::endl;
          std::cout << "    Position (X,Y,Z): (" << fx << ", " << fy << ", " << fz << ") cm" << std::endl;
          std::cout << "    Distance from LoS: " << fdist << " cm" << std::endl;
          std::cout << "    Dimensions: " << 2*fBox->GetDX() << " x " << 2*fBox->GetDY() << " x " << 2*fBox->GetDZ() << " cm³" << std::endl;
          std::cout << "  ╚══════════════════════════╝" << std::endl;
        }
        
        if (lastEndplate && lastEndplate != firstEndplate) {
          std::cout << "\n  ╔═══ LAST MDT ENDPLATE ═══╗" << std::endl;
          const Double_t *lTrans = lastEndplateMat.GetTranslation();
          Double_t lx = lTrans[0], ly = lTrans[1], lz = lTrans[2];
          Double_t ldist = TMath::Sqrt(lx*lx + ly*ly);
          
          TGeoVolume *lVol = lastEndplate->GetVolume();
          TGeoBBox *lBox = (TGeoBBox*)lVol->GetShape();
          
          std::cout << "    Name: " << lastEndplate->GetName() << std::endl;
          std::cout << "    Position (X,Y,Z): (" << lx << ", " << ly << ", " << lz << ") cm" << std::endl;
          std::cout << "    Distance from LoS: " << ldist << " cm" << std::endl;
          std::cout << "    Dimensions: " << 2*lBox->GetDX() << " x " << 2*lBox->GetDY() << " x " << 2*lBox->GetDZ() << " cm³" << std::endl;
          std::cout << "  ╚═════════════════════════╝" << std::endl;
        }
        
        // Set transparency for muon spectrometer
        SetVolumeTransparencyRecursive(vol, 60);
        
        // Draw orange dashed line through MuonSpectrometer (detector axis)
        TEveLine* muonLine = new TEveLine("MuonSpectLine");
        muonLine->SetLineColor(kOrange);
        muonLine->SetLineWidth(3);
        muonLine->SetLineStyle(2);
        // Line extends in detector direction
        muonLine->SetNextPoint(x - det_dir_x * lineLength/2, y - det_dir_y * lineLength/2, z - det_dir_z * lineLength/2);
        muonLine->SetNextPoint(x + det_dir_x * lineLength/2, y + det_dir_y * lineLength/2, z + det_dir_z * lineLength/2);
        gEve->AddElement(muonLine);
        
        // Draw lines from LoS to upstream and downstream faces
        TEveLine* upstreamLine = new TEveLine("MuonSpect_UpstreamToLoS");
        upstreamLine->SetLineColor(kOrange+2);
        upstreamLine->SetLineWidth(2);
        upstreamLine->SetLineStyle(1);
        upstreamLine->SetNextPoint(0, 0, z_start);
        upstreamLine->SetNextPoint(x_start, y_start, z_start);
        gEve->AddElement(upstreamLine);
        
        TEveLine* downstreamLine = new TEveLine("MuonSpect_DownstreamToLoS");
        downstreamLine->SetLineColor(kOrange-3);
        downstreamLine->SetLineWidth(2);
        downstreamLine->SetLineStyle(1);
        downstreamLine->SetNextPoint(0, 0, z_end);
        downstreamLine->SetNextPoint(x_end, y_end, z_end);
        gEve->AddElement(downstreamLine);
        
        // Add markers at face positions
        TEvePointSet* upstreamPoint = new TEvePointSet("MuonSpect_UpstreamFace");
        upstreamPoint->SetNextPoint(x_start, y_start, z_start);
        upstreamPoint->SetMarkerStyle(20);
        upstreamPoint->SetMarkerSize(1.5);
        upstreamPoint->SetMarkerColor(kOrange+2);
        gEve->AddElement(upstreamPoint);
        
        TEvePointSet* downstreamPoint = new TEvePointSet("MuonSpect_DownstreamFace");
        downstreamPoint->SetNextPoint(x_end, y_end, z_end);
        downstreamPoint->SetMarkerStyle(20);
        downstreamPoint->SetMarkerSize(1.5);
        downstreamPoint->SetMarkerColor(kOrange-3);
        gEve->AddElement(downstreamPoint);
      }
    }
    
    // Report what was found
    std::cout << "\n=== Summary ===" << std::endl;
    if (!found3DCAL) std::cout << "WARNING: 3DCAL (ContainerPlacement/ContainerBox) not found" << std::endl;
    if (!foundECAL) std::cout << "WARNING: ECAL (rearCal/ContainerEcal) not found" << std::endl;
    if (!foundAHCAL) std::cout << "WARNING: AHCAL (rearHCal/ContainerHcal) not found" << std::endl;
    if (!foundMuonSpect) std::cout << "WARNING: Muon Spectrometer not found" << std::endl;
    
    std::cout << "\n=== Visualization Legend ===" << std::endl;
    std::cout << "Blue line:         Line of Sight (LoS) along Z-axis" << std::endl;
    std::cout << "Green dashed:      3DCAL detector axis (tilted)" << std::endl;
    std::cout << "Magenta dashed:    ECAL detector axis (tilted)" << std::endl;
    std::cout << "Cyan dashed:       AHCAL detector axis (tilted)" << std::endl;
    std::cout << "Orange dashed:     Muon Spectrometer detector axis" << std::endl;
    std::cout << "\nSolid lines:       Distance from LoS to detector faces" << std::endl;
    std::cout << "  Lighter colors:  Upstream face (beam entrance)" << std::endl;
    std::cout << "  Darker colors:   Downstream face (beam exit)" << std::endl;
    std::cout << "\nDot markers:       Exact face positions" << std::endl;
    std::cout << "Red arrows:        Neutrino flux directions" << std::endl;
    std::cout << "=============================\n" << std::endl;
  }

  // ----------------------------------------------------
  // 5. Open neutrino flux and draw arrows
  // ----------------------------------------------------
  TFile *fflux = TFile::Open("events_light_4x4.root");
  if (!fflux || fflux->IsZombie()) {
    Error("Geometry_flux_fasernu", "Cannot open events_light_4x4.root");
    gEve->Redraw3D(kTRUE);
    return;
  }

  TTree *t = (TTree*) fflux->Get("flux");
  if (!t) {
    Error("Geometry_flux_fasernu", "Tree 'flux' not found in events_light_4x4.root");
    gEve->Redraw3D(kTRUE);
    return;
  }


// =====================================================
//  DRAW LINE OF SIGHT (LOS)
// =====================================================

// Define LOS direction (unit vector)
Double_t los_x = 0;
Double_t los_y = 0;
Double_t los_z = 1;

// Scale LOS to visible long line
Double_t L_los = 2000.0;  // length in cm (20 m) – choose as needed

// Starting point (IP)
Double_t ipx = 0.0;
Double_t ipy = 0.0;
Double_t ipz = 0.0;

// Create LOS line
TEveLine* los = new TEveLine("LineOfSight");
los->SetLineColor(kBlue);
los->SetLineWidth(3);

// Add two points: IP → IP + L * direction
los->SetNextPoint(ipx, ipy, ipz);
los->SetNextPoint(ipx + los_x * L_los,
                  ipy + los_y * L_los,
                  ipz + los_z * L_los);

// Add to TEve scene
gEve->AddElement(los);


  
  Double_t vtxx, vtxy, vtxz;   // starting point
  Double_t px,  py,  pz;       // direction / momentum

  t->SetBranchAddress("vtxx", &vtxx);
  t->SetBranchAddress("vtxy", &vtxy);
  t->SetBranchAddress("vtxz", &vtxz);
  t->SetBranchAddress("px",   &px);
  t->SetBranchAddress("py",   &py);
  t->SetBranchAddress("pz",   &pz);

  Int_t    nEntries  = t->GetEntries();
  Int_t    maxToDraw = 10000;   // limit for performance
  Double_t L         = 10.0;    // arrow length in geometry units (e.g. cm)

  for (Int_t i = 0; i < nEntries && i < maxToDraw; ++i) {
    t->GetEntry(i);

    // if vertices are stored in meters, convert to cm:
    vtxx *= 100.0;
    vtxy *= 100.0;
    vtxz *= 100.0;

    // Normalize momentum to unit vector for direction
    Double_t p = TMath::Sqrt(px*px + py*py + pz*pz);
    if (p == 0) continue;

    Double_t dx = (px/p) * L;
    Double_t dy = (py/p) * L;
    Double_t dz = (pz/p) * L;

    TEveLine *line = new TEveLine(2);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetNextPoint(vtxx,        vtxy,        vtxz       );
    line->SetNextPoint(vtxx + dx,   vtxy + dy,   vtxz + dz );

    gEve->AddElement(line);
  }

  // ----------------------------------------------------
  // 6. Show TEve GUI
  // ----------------------------------------------------
  gEve->Redraw3D(kTRUE);
  // Get the GL viewer
TGLViewer* glv = gEve->GetDefaultGLViewer();

// Print instructions
std::cout << "\n=== Event Viewer Controls ===" << std::endl;
std::cout << "Move mouse over 3D view to see coordinates in status bar" << std::endl;
std::cout << "Right-click on objects for detailed information" << std::endl;
std::cout << "Use mouse wheel to zoom, drag to rotate" << std::endl;
  
}

/////////////////

void ExtractFaserCalDistancesFromLoS()
{
  // ----------------------------------------------------
  // Extract distances from Line of Sight to FASERCAL subdetectors
  // LoS: z-axis from origin (0,0,0) in direction (0,0,1)
  // ----------------------------------------------------
  
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  // ----------------------------------------------------
  // 1. Import main geometry (FaserNu3)
  // ----------------------------------------------------
  TGeoManager::Import("FaserNu3.gdml");
  TGeoManager *geom = gGeoManager;
  if (!geom) {
    Error("ExtractFaserCalDistancesFromLoS", "Cannot load FaserNu3.gdml");
    return;
  }

  // ----------------------------------------------------
  // 2. Import FASERCAL_V8 as a volume and add to geometry
  // ----------------------------------------------------
  TGDMLParse parser;
  TGeoVolume *vol2 = parser.GDMLReadFile("/Users/ukose/sw/kits/NewFASER/FASER/FASERG4/build/FASERCAL_V9.gdml");
  if (!vol2) {
    Error("ExtractFaserCalDistancesFromLoS", "Cannot load /Users/ukose/sw/kits/NewFASER/FASER/FASERG4/build/FASERCAL_V9.gdml");
    return;
  }

  vol2->SetName("CalGeometry");

  // Position FASERCAL at z = 600 cm (same as visualization function)
  TGeoTranslation *tr = new TGeoTranslation("tr", 0.0, 0.0, 600.0);
  tr->RegisterYourself();

  geom->GetTopVolume()->AddNode(vol2, 1, tr);
  geom->CloseGeometry();

  // ----------------------------------------------------
  // 3. Navigate to CalGeometry and extract subdetector distances
  // ----------------------------------------------------
  std::cout << "\n========================================" << std::endl;
  std::cout << "FASERCAL Subdetector Distances from LoS" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "LoS: z-axis from origin (0,0,0) in direction (0,0,1)" << std::endl;
  std::cout << "FASERCAL offset: (0, 0, 600) cm\n" << std::endl;

  // Method 1: Try to find the node by name
  TGeoNode *calNode = nullptr;
  TGeoVolume *topVol = geom->GetTopVolume();
  Int_t nNodes = topVol->GetNdaughters();
  
  for (Int_t i = 0; i < nNodes; ++i) {
    TGeoNode *node = topVol->GetNode(i);
    if (TString(node->GetName()).Contains("CalGeometry")) {
      calNode = node;
      break;
    }
  }

  if (!calNode) {
    Error("ExtractFaserCalDistancesFromLoS", "CalGeometry node not found in top volume");
    return;
  }

  TGeoVolume *calVol = calNode->GetVolume();
  std::cout << "Found CalGeometry with " << calVol->GetNdaughters() << " daughter volumes\n" << std::endl;

  // Get the transformation matrix of CalGeometry node
  TGeoMatrix *calMatrix = calNode->GetMatrix();
  const Double_t *calTrans = calMatrix->GetTranslation();

  // ----------------------------------------------------
  // 4. Iterate through all subdetectors and calculate distances
  // ----------------------------------------------------
  Int_t nDaughters = calVol->GetNdaughters();
  
  std::cout << std::setw(5) << "Index"
            << std::setw(40) << "Subdetector Name"
            << std::setw(12) << "X (cm)"
            << std::setw(12) << "Y (cm)"
            << std::setw(12) << "Z (cm)"
            << std::setw(15) << "Dist from LoS (cm)" << std::endl;
  std::cout << std::string(96, '-') << std::endl;

  for (Int_t i = 0; i < nDaughters; ++i) {
    TGeoNode *node = calVol->GetNode(i);
    TGeoMatrix *matrix = node->GetMatrix();
    
    // Get local position within CalGeometry
    const Double_t *localTrans = matrix->GetTranslation();
    
    // Calculate global position (local + CalGeometry offset)
    Double_t x = localTrans[0] + calTrans[0];
    Double_t y = localTrans[1] + calTrans[1];
    Double_t z = localTrans[2] + calTrans[2];
    
    // Calculate perpendicular distance to LoS (z-axis)
    // For LoS along z-axis from origin: distance = sqrt(x^2 + y^2)
    Double_t distToLoS = TMath::Sqrt(x*x + y*y);
    
    std::cout << std::setw(5) << i
              << std::setw(40) << node->GetName()
              << std::setw(12) << std::fixed << std::setprecision(2) << x
              << std::setw(12) << std::fixed << std::setprecision(2) << y
              << std::setw(12) << std::fixed << std::setprecision(2) << z
              << std::setw(15) << std::fixed << std::setprecision(2) << distToLoS << std::endl;
  }

  std::cout << "\n========================================" << std::endl;
  std::cout << "Total subdetectors: " << nDaughters << std::endl;
  std::cout << "========================================\n" << std::endl;
  
  // ----------------------------------------------------
  // 5. Navigate into DetectorAssembly to find MuonSpectrometer and other modules
  // ----------------------------------------------------
  std::cout << "\n=== Searching in DetectorAssembly ===" << std::endl;
  
  // Find DetectorAssembly node (should be the only daughter of CalGeometry)
  TGeoNode *assemblyNode = nullptr;
  if (nDaughters == 1) {
    assemblyNode = calVol->GetNode(0);
    std::cout << "Found: " << assemblyNode->GetName() << std::endl;
  }
  
  if (!assemblyNode) {
    std::cout << "ERROR: DetectorAssembly node not found" << std::endl;
    return;
  }
  
  TGeoVolume *assemblyVol = assemblyNode->GetVolume();
  TGeoMatrix *assemblyMatrix = assemblyNode->GetMatrix();
  const Double_t *assemblyTrans = assemblyMatrix->GetTranslation();
  
  Int_t nAssemblyModules = assemblyVol->GetNdaughters();
  std::cout << "DetectorAssembly has " << nAssemblyModules << " modules\n" << std::endl;
  
  // Search for MuonSpectrometer in DetectorAssembly
  Bool_t foundMuonSpect = kFALSE;
  
  std::cout << "\nSearching for MuonSpectrometer among " << nAssemblyModules << " modules:" << std::endl;
  for (Int_t i = 0; i < nAssemblyModules; ++i) {
    TGeoNode *node = assemblyVol->GetNode(i);
    TString nodeName(node->GetName());
    std::cout << "  Module " << i << ": " << nodeName;
    if (nodeName.Contains("Muon") || nodeName.Contains("muon") || nodeName.Contains("MDT")) std::cout << " <- MUON CANDIDATE";
    std::cout << std::endl;
    
    if (nodeName.Contains("MuonSpectrometer") || nodeName.Contains("Muon") || 
        nodeName.Contains("muon") || nodeName.Contains("MDT")) {
      foundMuonSpect = kTRUE;
      std::cout << "\n>>> Found MuonSpectrometer: " << nodeName << " <<<\n" << std::endl;
      TGeoMatrix *matrix = node->GetMatrix();
      const Double_t *localTrans = matrix->GetTranslation();
      
      // Get combined transformation
      TGeoHMatrix combinedMat;
      combinedMat = *assemblyMatrix;
      combinedMat.Multiply(matrix);
      
      // Apply CalGeometry offset
      TGeoHMatrix finalMat;
      TGeoTranslation calTransMat(calTrans[0], calTrans[1], calTrans[2]);
      finalMat = calTransMat;
      finalMat.Multiply(&combinedMat);
      
      // Get actual global position
      const Double_t *globalTrans = finalMat.GetTranslation();
      Double_t x = globalTrans[0];
      Double_t y = globalTrans[1];
      Double_t z = globalTrans[2];
      Double_t distToLoS = TMath::Sqrt(x*x + y*y);
      
      // Get dimensions
      TGeoVolume *vol = node->GetVolume();
      TGeoBBox *bbox = (TGeoBBox*)vol->GetShape();
      Double_t dx = bbox->GetDX();
      Double_t dy = bbox->GetDY();
      Double_t dz = bbox->GetDZ();
      
      std::cout << "\n╔═══ MUON SPECTROMETER ═══╗" << std::endl;
      std::cout << "  Name: " << nodeName << std::endl;
      std::cout << "  Center Position (X,Y,Z): (" << x << ", " << y << ", " << z << ") cm" << std::endl;
      std::cout << "  Perpendicular distance from LoS: " << distToLoS << " cm" << std::endl;
      std::cout << "  Dimensions (full): " << 2*dx << " x " << 2*dy << " x " << 2*dz << " cm³" << std::endl;
      std::cout << "  Volume: " << 8*dx*dy*dz/1e6 << " m³" << std::endl;
      
      // List all tracking layers
      Int_t nLayers = vol->GetNdaughters();
      std::cout << "\n  Tracking Layers/Stations: " << nLayers << std::endl;
      std::cout << "  " << std::string(95, '-') << std::endl;
      std::cout << "  " << std::setw(3) << "ID"
                << std::setw(5) << "Type"
                << std::setw(40) << "Layer Name"
                << std::setw(12) << "Z (cm)"
                << std::setw(18) << "Dist from LoS (cm)" << std::endl;
      std::cout << "  " << std::string(95, '-') << std::endl;
      
      Int_t nMagnets = 0;
      Int_t nEndplates = 0;
      Int_t nSciFi = 0;
      
      // Track first and last MDTEndplate
      TGeoNode *firstEndplate = nullptr;
      TGeoNode *lastEndplate = nullptr;
      TGeoHMatrix firstEndplateMat;
      TGeoHMatrix lastEndplateMat;
      
      for (Int_t j = 0; j < nLayers; ++j) {
        TGeoNode *layerNode = vol->GetNode(j);
        TString layerName(layerNode->GetName());
        TGeoMatrix *layerMatrix = layerNode->GetMatrix();
        
        // Transform to global coordinates
        TGeoHMatrix layerFinalMat = finalMat;
        layerFinalMat.Multiply(layerMatrix);
        const Double_t *layerGlobalTrans = layerFinalMat.GetTranslation();
        
        Double_t layer_x = layerGlobalTrans[0];
        Double_t layer_y = layerGlobalTrans[1];
        Double_t layer_z = layerGlobalTrans[2];
        Double_t layer_dist = TMath::Sqrt(layer_x*layer_x + layer_y*layer_y);
        
        // Identify component type and track endplates
        TString type = "    ";
        if (layerName.Contains("MDTMagnet")) {
          type = "MAG ";
          nMagnets++;
        } else if (layerName.Contains("MDTEndplate")) {
          type = "END ";
          nEndplates++;
          if (!firstEndplate) {
            firstEndplate = layerNode;
            firstEndplateMat = layerFinalMat;
          }
          lastEndplate = layerNode;
          lastEndplateMat = layerFinalMat;
        } else if (layerName.Contains("SciFi")) {
          type = "SCI ";
          nSciFi++;
        }
        
        std::cout << "  " << std::setw(3) << j
                  << std::setw(5) << type
                  << std::setw(40) << layerNode->GetName()
                  << std::setw(12) << std::fixed << std::setprecision(2) << layer_z
                  << std::setw(18) << layer_dist << std::endl;
      }
      std::cout << "  " << std::string(95, '-') << std::endl;
      std::cout << "  Summary: " << nMagnets << " MDTMagnets, " 
                << nEndplates << " MDTEndplates, " 
                << nSciFi << " SciFi layers" << std::endl;
      std::cout << "  (MAG=Magnet, END=Endplate, SCI=SciFi)" << std::endl;
      
      // Report detailed info on first and last MDTEndplate
      if (firstEndplate) {
        std::cout << "\n  ╔═══ FIRST MDT ENDPLATE ═══╗" << std::endl;
        const Double_t *fTrans = firstEndplateMat.GetTranslation();
        Double_t fx = fTrans[0], fy = fTrans[1], fz = fTrans[2];
        Double_t fdist = TMath::Sqrt(fx*fx + fy*fy);
        
        TGeoVolume *fVol = firstEndplate->GetVolume();
        TGeoBBox *fBox = (TGeoBBox*)fVol->GetShape();
        
        std::cout << "    Name: " << firstEndplate->GetName() << std::endl;
        std::cout << "    Position (X,Y,Z): (" << fx << ", " << fy << ", " << fz << ") cm" << std::endl;
        std::cout << "    Distance from LoS: " << fdist << " cm" << std::endl;
        std::cout << "    Dimensions: " << 2*fBox->GetDX() << " x " << 2*fBox->GetDY() << " x " << 2*fBox->GetDZ() << " cm³" << std::endl;
        std::cout << "  ╚══════════════════════════╝" << std::endl;
      }
      
      if (lastEndplate && lastEndplate != firstEndplate) {
        std::cout << "\n  ╔═══ LAST MDT ENDPLATE ═══╗" << std::endl;
        const Double_t *lTrans = lastEndplateMat.GetTranslation();
        Double_t lx = lTrans[0], ly = lTrans[1], lz = lTrans[2];
        Double_t ldist = TMath::Sqrt(lx*lx + ly*ly);
        
        TGeoVolume *lVol = lastEndplate->GetVolume();
        TGeoBBox *lBox = (TGeoBBox*)lVol->GetShape();
        
        std::cout << "    Name: " << lastEndplate->GetName() << std::endl;
        std::cout << "    Position (X,Y,Z): (" << lx << ", " << ly << ", " << lz << ") cm" << std::endl;
        std::cout << "    Distance from LoS: " << ldist << " cm" << std::endl;
        std::cout << "    Dimensions: " << 2*lBox->GetDX() << " x " << 2*lBox->GetDY() << " x " << 2*lBox->GetDZ() << " cm³" << std::endl;
        std::cout << "  ╚═════════════════════════╝" << std::endl;
      }
      
      std::cout << "╚═════════════════════════╝\n" << std::endl;
      break;
    }
  }
  
  if (!foundMuonSpect) {
    std::cout << "Note: Muon Spectrometer not found in DetectorAssembly" << std::endl;
  }
  std::cout << "======================================\n" << std::endl;
}

/////////////////

void ViewCombinedGeometry()
{
  // ----------------------------------------------------
  // Load and visualize the exported combined geometry ROOT file
  // ----------------------------------------------------
  
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libEve");
  
  // Open the ROOT file containing the geometry
  TFile *f = TFile::Open("combined_fasernu_fasercal.root");
  if (!f || f->IsZombie()) {
    Error("ViewCombinedGeometry", "Cannot open combined_fasernu_fasercal.root");
    std::cout << "\nMake sure you've run Geometry_flux_fasernu_fasercal() first to export the file." << std::endl;
    return;
  }
  
  // List what's in the file
  std::cout << "\n=== Contents of ROOT file ===" << std::endl;
  f->ls();
  std::cout << "==============================\n" << std::endl;
  
  // Get the geometry manager from the file
  // Try different possible key names
  TGeoManager *geom = nullptr;
  
  // First try to get the first TGeoManager object
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if (TString(key->GetClassName()).Contains("TGeoManager")) {
      geom = (TGeoManager*)key->ReadObj();
      std::cout << "Found TGeoManager with key: " << key->GetName() << std::endl;
      break;
    }
  }
  
  if (!geom) {
    Error("ViewCombinedGeometry", "Cannot find TGeoManager object in file");
    std::cout << "File contents listed above. Check if geometry was exported correctly." << std::endl;
    f->Close();
    return;
  }
  
  std::cout << "\n=== Loaded Combined Geometry ===" << std::endl;
  std::cout << "Top Volume: " << geom->GetTopVolume()->GetName() << std::endl;
  std::cout << "Number of volumes: " << geom->GetListOfVolumes()->GetEntries() << std::endl;
  std::cout << "Number of materials: " << geom->GetListOfMaterials()->GetEntries() << std::endl;
  std::cout << "================================\n" << std::endl;
  
  // Create TEve Manager for visualization
  TEveManager::Create();
  
  // Import geometry into TEve
  TEveGeoTopNode *topNode = new TEveGeoTopNode(geom, geom->GetTopNode());
  topNode->SetVisLevel(4);  // Adjust detail level (2-5)
  
  gEve->AddGlobalElement(topNode);
  
  // Draw Line of Sight
  TEveLine* los = new TEveLine("LineOfSight");
  los->SetLineColor(kBlue);
  los->SetLineWidth(3);
  los->SetNextPoint(0, 0, 0);
  los->SetNextPoint(0, 0, 2000);
  gEve->AddElement(los);
  
  std::cout << "\n=== Viewer Controls ===" << std::endl;
  std::cout << "Mouse: Drag to rotate, wheel to zoom" << std::endl;
  std::cout << "Keys: r=reset view, w=wireframe, e=axes" << std::endl;
  std::cout << "Right-click objects for info" << std::endl;
  std::cout << "========================\n" << std::endl;
  
  // Show TEve GUI
  gEve->Redraw3D(kTRUE);
}

/////////////////

void ViewCombinedGeometrySimple()
{
  // ----------------------------------------------------
  // Simple view without TEve (uses TGeoManager's built-in viewer)
  // ----------------------------------------------------
  
  gSystem->Load("libGeom");
  
  TFile *f = TFile::Open("combined_fasernu_fasercal.root");
  if (!f || f->IsZombie()) {
    Error("ViewCombinedGeometrySimple", "Cannot open combined_fasernu_fasercal.root");
    return;
  }
  
  // Find TGeoManager object
  TGeoManager *geom = nullptr;
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if (TString(key->GetClassName()).Contains("TGeoManager")) {
      geom = (TGeoManager*)key->ReadObj();
      break;
    }
  }
  
  if (!geom) {
    Error("ViewCombinedGeometrySimple", "Cannot find TGeoManager object in file");
    f->Close();
    return;
  }
  
  std::cout << "Loaded geometry: " << geom->GetTopVolume()->GetName() << std::endl;
  
  // Use ROOT's built-in geometry viewer
  geom->SetVisLevel(4);
  geom->GetTopVolume()->Draw("ogl");  // OpenGL viewer
  
  std::cout << "\nViewer opened. Use mouse to rotate and zoom." << std::endl;
}


