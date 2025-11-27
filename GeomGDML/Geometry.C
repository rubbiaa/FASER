#include "TGeoManager.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"

#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TSystem.h"
#include "TMath.h"


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
  TGeoVolume *vol2 = parser.GDMLReadFile("geometry_tilted_5degree.gdml");
  if (!vol2) {
    Error("Geometry_flux_fasernu", "Cannot load geometry_tilted_5degree.gdml as volume");
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
  TGDMLParse parser2;
  TGeoVolume *vol3 = parser.GDMLReadFile("geometry_tilted_5degree.gdml");
  if (!vol3) {
    Error("Geometry_flux_fasernu", "Cannot load geometry_tilted_5degree.gdml as volume");
    return;
  }

  // Optionally rename to avoid name clashes
  vol3->SetName("CalGeometry2");

  // Set position / rotation of the second geometry relative to FaserNu3
  // (currently at origin; change x,y,z and rotations as needed)
  // TGeoTranslation *tr = new TGeoTranslation("tr_cal", 11.0, 17.0, 820.0); // in cm
  //TGeoTranslation *tr3 = new TGeoTranslation("tr_cal", 20.0, 17.0, 600.0); // in cm
  TGeoTranslation *tr3 = new TGeoTranslation("tr_cal", 0.0, 0.0, 0.0); // in cm
  tr3->RegisterYourself();  // good practice if you later use in matrices

  //geom->GetTopVolume()->AddNode(vol3, 1, tr3);
  geom->CloseGeometry();

  // ----------------------------------------------------
  // 4. Start TEve and show the combined geometry
  // ----------------------------------------------------
  TEveManager::Create();

  TEveGeoTopNode *topNode =
      new TEveGeoTopNode(geom, geom->GetTopNode());
  topNode->SetVisLevel(4);
  gEve->AddGlobalElement(topNode);

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
  
}
