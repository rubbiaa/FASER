#include <TGFrame.h>
#include <TGWindow.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>

#include <TcalEvent.hh>

class MyMainFrame : public TGMainFrame {
    RQ_OBJECT("MyMainFrame")
public:
    MyMainFrame(TcalEvent *ft, const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~MyMainFrame();

    void HandleButton(); // Function to handle button click

    void ZoomToPosition(Double_t x, Double_t y, Double_t z);

    ClassDef(MyMainFrame,1)

private:
    TGMainFrame *fMain;
    TGTextButton *fButton;
    TRootEmbeddedCanvas *fCanvas;

    TcalEvent* fTcalEvent;

};
