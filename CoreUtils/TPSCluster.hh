#ifndef _TPSCLUSTER_H
#define _TPSCLUSTER_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TObject.h"

#include "TcalEvent.hh"

/// @brief A 2D cluster generated from one view (XZ or YZ) in the plastic scintillator
class TPSCluster : public TObject
{
private:
    TcalEvent* fTcalEvent;
public:
    /// @brief Longitudinal profile \frac{dE}{dt} = E_0 \cdot b \cdot \frac{(b t)^{a-1} \cdot e^{-b t}}{\Gamma(a)} (see PDG)
    struct PSCLUSTERLONGPROFILE {
        double chi2;                // fitted chi2 
        double chi2_per_ndf;
        double E0;                  // fitted total energy
        double a;                   // fitted a parameter of longitudinal profile
        double b;                   // fitted b
        double tmax;                // position of max of shower (in rad length) = (a-1)/b
        double y;                   // y = E/Ec where Ec is an adhoc critical energy
        double c;                   // c = tmax - ln(y)
    };

    struct PSCLUSTERHIT
    {
        long id;
        float EDeposit;  // MeV
    };

    int clusterID;
    int view;                                    // =0 for XZ and =1 for YZ
    double rawenergy = 0; // in MeV !
    ROOT::Math::XYZVector vtx;                   // the primary vertex of the cluster (must be set manually)
    ROOT::Math::XYZVector cog;                   // the center of gravity of the cluster
    struct PSCLUSTERLONGPROFILE longenergyprofile; // fitted parameters of energy longidutinal profile

    std::vector<struct PSCLUSTERHIT> hits;

    TPSCluster() : fTcalEvent(0), view(0) {};
    TPSCluster(int v, TcalEvent* ft) { fTcalEvent = ft; view = v;};

    ClassDef(TPSCluster, 1)

    /// @brief Compute center of gravity of cluster
    void ComputeCOG();

    /// @brief Set the assumed creation point of the cluster
    /// @param x, y, z : Origin of cluster, used for example to compute longidutinal profile
    void setVtx(double x, double y, double z) { vtx.SetX(x); vtx.SetY(y); vtx.SetZ(z); };

    /// @brief Compute longitudinal profile of shower defined by the cluster
    void ComputeLongProfile(int verbose = 0);

};

#endif
