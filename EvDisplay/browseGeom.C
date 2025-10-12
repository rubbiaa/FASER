gSystem->Load("libGeom");
TGeoManager::Import("../GeomGDML/geometry.gdml");
gGeoManager->GetTopVolume()->Draw("ogl");
