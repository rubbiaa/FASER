#ifndef _TMUONSPECT_HH_
#define _TMUONSPECT_HH_ 1

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

#include <TTree.h>

#define MAXMUTRACKS 10

class TMuonSpectrometer : public TObject {

public:
struct FEATURES {
    int ntracks;        // number of tracks in the muon spectrometer
    float charge[MAXMUTRACKS]; // the charges of the tracks
    int npoints[MAXMUTRACKS]; // number of points in each track
    float px[MAXMUTRACKS];    // the fitted momentum x component at the first point
    float py[MAXMUTRACKS];    // the fitted momentum y component at the first point
    float pz[MAXMUTRACKS];    // the fitted momentum z component at the first point
    float p[MAXMUTRACKS];     // the fitted momentum magnitude at the first point
    float chi2[MAXMUTRACKS];  // chi2 of the fit
    int nDoF[MAXMUTRACKS];    // nDoF of the fit
    float pval[MAXMUTRACKS];  // p-value of the fit
};

private:

public:
    ClassDef(TMuonSpectrometer, 1)

    struct FEATURES features;

    void Create_Sel_Tree(TTree *t);
    void Fill_Sel_Tree(TTree *t) { t->Fill();};

};


#endif
