#include <iostream>
#include <cmath>

#include "TPSCluster.hh"

#include <TVector3.h>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TMinuit.h"
#include "Fit/FitResult.h"

ClassImp(TPSCluster);

// Gamma function using std::tgamma
static double gamma_function(double a) {
    return std::tgamma(a);
}

// Differential energy deposition function
static double dE_dt(double t, double E0, double a, double b) {
    if (t < 0) return 0.0; // Energy deposition is zero for negative time/depth
    double term1 = std::pow(b * t, a - 1);
    double term2 = std::exp(-b * t);
    double gamma_a = gamma_function(a);
    return E0 * b * (term1 * term2 / gamma_a);
}

// Define the fitting function
static double energyProfileFunc(double *x, double *par) {
    double t = x[0];          // Independent variable (e.g., bin center)
    double E0 = par[0];       // E0 parameter
    double a = par[1];        // a parameter (shape)
    double b = par[2];        // b parameter (scale)
    
    double gamma_a = TMath::Gamma(a);
    double term1 = pow(b * t, a - 1);
    double term2 = exp(-b * t);
    
    return E0 * b * (term1 * term2 / gamma_a);
}

void TPSCluster::ComputeCOG() {
    double cogx = 0;
    double cogy = 0;
    double cogz = 0;
    double etot = 0;

    for (auto &it : hits) {
        long ID = it.id;
        float ehit = it.EDeposit;
        ROOT::Math::XYZVector position = fTcalEvent -> getChannelXYZfromID(ID);
        cogx += position.X()*ehit;
        cogy += position.Y()*ehit;
        cogz += position.Z()*ehit;
        etot += ehit;
    }
    cog.SetX(cogx/etot);
    cog.SetY(cogy/etot);
    cog.SetZ(cogz/etot);
}

void TPSCluster::ComputeLongProfile(int verbose) {

    const int nBins = 30;
    double mint = 0.0;
    double maxt = 30.0;    // in radiation lengths
    double binWidth = (maxt - mint) / nBins;
    double conversionX0 = 0.78/52.2;  /// in X0/mm MAKE SURE THIS IS CORRECT IF CHANGE GEOMETRY !!!

    // Initialize an array to store the energy in each bin
    std::vector<double> energyProfile(nBins, 0.0);

    ROOT::Math::XYZVector directionVec = cog - vtx;
    TVector3 direction((view==0) ? directionVec.X() : 0, (view!=0) ? directionVec.Y() : 0, directionVec.Z());
    direction = direction.Unit();

    for (auto &it : hits)
    {
        long ID = it.id;
        float ehit = it.EDeposit;
        ROOT::Math::XYZVector position = fTcalEvent->getChannelXYZfromID(ID);
        ROOT::Math::XYZVector hitdirection = position - vtx;
        if( view == 0 ) { hitdirection.SetY(0); } else { hitdirection.SetX(0); };
        TVector3 hitDirVec(hitdirection.X(), hitdirection.Y(), hitdirection.Z());
        double t = hitDirVec.Dot(direction)*conversionX0;
        int binIndex = static_cast<int>((t - mint) / binWidth);        
        if (binIndex >= 0 && binIndex < nBins) {
            energyProfile[binIndex] += ehit;
        }
    }

    if(verbose > 1) {
     for (int i = 0; i < nBins; ++i) {
        std::cout << "t: " << i*binWidth << ": " << energyProfile[i] << " MeV" << std::endl;
    }
    }

    // profile histrogram
    TH1D* hist = (TH1D*)gDirectory->Get("TPSclusterenergyProfile");
    if(hist != nullptr) {
        hist->Reset();
    } else {
       hist = new TH1D("TPSclusterenergyProfile", "Cluster Longitudinal Energy Profile", nBins, 0, nBins);
    } 
    
    for (int i = 0; i < nBins; ++i) {
        hist->SetBinContent(i+1, energyProfile[i]);
    }

    // Create the TF1 object using the fitting function
    TF1* fitFunc = (TF1*)gDirectory->Get("TPSclusterenergyfitFunc");
    if(fitFunc == nullptr)
        fitFunc = new TF1("TPSclusterenergyfitFunc", energyProfileFunc, 0, nBins, 3);

    // Set initial parameter guesses
    fitFunc->SetParameters(rawenergy, 2.0, 0.6);  
    fitFunc->FixParameter(2, 0.6);

    // Perform the fit and capture the result
    TFitResultPtr fitStatus = hist->Fit("TPSclusterenergyfitFunc", "S R"); // "S" option is required to get the fit status
    if (fitStatus != 0) {
        std::cout << "Fit failed with status: " << fitStatus << std::endl;
        // Handle the failure case
    }

    // Extract the fit parameters
    double E0 = fitFunc->GetParameter(0);
    double a = fitFunc->GetParameter(1);
    double b = fitFunc->GetParameter(2);
    double tmax = (a-1)/b;
    ///////////////////////////////////////////////////////////////////////
    double critical_energy = 23.82; // ad-hoc  (W 8 MeV)
    ///////////////////////////////////////////////////////////////////////
    double y = rawenergy / critical_energy;
    double c = tmax-std::log(y);

    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double chi2_per_ndf = chi2 / ndf;

    longenergyprofile.chi2 = chi2;
    longenergyprofile.chi2_per_ndf = chi2_per_ndf;
    longenergyprofile.E0 = E0;
    longenergyprofile.a = a;
    longenergyprofile.b = b;
    longenergyprofile.tmax = tmax;
    longenergyprofile.y = y;
    longenergyprofile.c = c;

    if(verbose < 1) return;

    std::cout << "Fitted parameters:" << std::endl;
    std::cout << "E0 = " << E0 << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "tmax = " << tmax << std::endl;
    std::cout << "c = " << c << std::endl;

#if 0
    // Visualize the results
    TCanvas *c1 = new TCanvas("TPSclusterenergycanvas", "Energy Profile Fit", 800, 600);
    hist->Draw();
    fitFunc->Draw("SAME");
#endif
}
