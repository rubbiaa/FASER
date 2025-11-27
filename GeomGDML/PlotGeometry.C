{
    // Load geometry
    TGeoManager::Import("geometry_tilted_5degree.gdml");
    gGeoManager->GetTopVolume()->Draw("ogl");
}
