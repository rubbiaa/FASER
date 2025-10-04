#include <TGFrame.h>
#include <TGWindow.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <TGeoVolume.h>
#include <RQ_OBJECT.h>
#include <TText.h>
#include <TGTextEntry.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TGButton.h>
#include <TChain.h>
#include <TGListTree.h>

#include <TcalEvent.hh>
#include "TPORecoEvent.hh"

class MyMainFrame : public TGMainFrame {
    RQ_OBJECT("MyMainFrame")
public:
    MyMainFrame(int run_number, int ievent, int mask, bool pre, const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~MyMainFrame();

    void Load_event(int run_number, int ievent, int mask);
    void Load_Recoevent(int run_number, int ievent);
    void Update_ListTree();
    void OnTrackSelected(TGListTreeItem *entry, int id);
    void Draw_event();
    void Draw_event_reco_tracks();
    void Draw_event_reco_voxel(TGeoShape *bigbox, TGeoMedium *air, TGeoShape *box);
    void Next_Event(int ievent);

    void HandleButton(); // Function to handle button click
    void toggle_prim_em();
    void toggle_prim_had();
    void toggle_sec_em();
    void toggle_sec_had();
    void toggle_reco_track();
    void toggle_reco_voxels();
    void toggle_color_fakes();
    void toggle_recon_ps_tracks();
    void only_reco();
    void next_event();
    void goto_event();

    void ZoomToPosition(Double_t x, Double_t y, Double_t z);
    void SideView();
    void ZoomIn();
    void ZoomOut();
    void MoveUp();
    void MoveDown();
    void MoveLeft();
    void MoveRight();

    ClassDef(MyMainFrame,1)

private:
    TGMainFrame *fMain;
    TGTextButton *fButton;
    TRootEmbeddedCanvas *fCanvas;
    TGListTree *listTree;
    TRootEmbeddedCanvas *fCanvas_2DPSview;
    TRootEmbeddedCanvas *fCanvas_2DPSviewZ;
    TRootEmbeddedCanvas *fCanvas_2DPSview_emhad;
    TRootEmbeddedCanvas *fCanvas_eldepo;
    TGTextEntry *textNextEventEntry;

    std::map<std::string, int> trackMap; // Maps track item names to their IDs
    int selectedTrackID = -1;            // Stores the selected track ID
    std::map<std::string, int> vertexMap; 
    int selectedVertexID = -1;            


    TText *runText = nullptr;
    TText *eventypeText = nullptr;
    TText *energyText = nullptr;
    TText *rearcalenergyText = nullptr;
    TText *rearmucalenergyText = nullptr;

    std::vector<TPolyLine3D*> polylineTracks;

    std::vector<TPolyLine3D*> polylineMagnetTracks;

    std::vector<TPolyMarker3D*> polylineVertices;

    int frun_number;
    int ievent;
    int event_mask;
    bool process_reco_event;
    bool opened_reco_event;
    TChain *reco_event_tree;

    TcalEvent* fTcalEvent;
    TPOEvent *POevent;
    TPORecoEvent* fPORecoEvent;

    TGeoShape *bigbox;

    // truth voxels
    TGeoVolume *primary_em;
    TGeoVolume *primary_had;
    TGeoVolume *secondary_em;
    TGeoVolume *secondary_had;
 
    // reconstructed tracks and vertices
    TGeoVolume *si_tracker;
    TGeoVolume *si_vertices;

    // reconstructed 3D voxel
    TGeoVolume *ps_reco_voxel;
    TGeoVolume *ps_tracks;

    // rear calorimeters
    TGeoVolume *rearcal;
    TGeoVolume *rearHcal;
    TGeoVolume *rearMucal;

    bool toggle_primary_em;
    bool toggle_primary_had;
    bool toggle_secondary_em;
    bool toggle_secondary_had;
    bool toggle_reco_voxel;
    bool toggle_color_fake_voxel;
    bool toggle_reconstructed_tracks;
    bool toggle_reconstructed_ps_tracks;

    TGCheckButton *f_fullreco_CheckBox;
    void on_fullreco_toggle(Bool_t state);

    Double_t range_min[3], range_max[3];
};
