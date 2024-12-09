
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooStats/SPlot.h"
#include "TMath.h"
#include "TTree.h"
#include "TLegend.h"
#include "THStack.h"
#include "TFile.h"
#include "TStyle.h"
#include <iostream>

TTree* event;

using namespace RooFit;

void signal_background_significance(std::string what, std::string cut = "") {
   // gSystem->Load("libRooFit");
   // gSystem->Load("libRooFitCore");

    double xmin = 0;
    double xmax = 4000.0;
    TH1F *h1 = new TH1F("h1", "nueCC", 10, xmin, xmax);
    TH1F *h2 = new TH1F("h2", "numuCC", 10, xmin, xmax);
    TH1F *h3 = new TH1F("h3", "NC", 10, xmin, xmax);
    event->Draw((what + std::string(">>h1")).c_str(), (cut + std::string("t_reaction == 1")).c_str());
    event->Draw((what + std::string(">>h2")).c_str(), (cut + std::string("t_reaction == 2")).c_str());
    event->Draw((what + std::string(">>h3")).c_str(), (cut + std::string("t_reaction > 3")).c_str());

    // Create signal and background histograms
    TH1F *h_signal = h1;
    TH1F *h_background = h2;
    h_background->Add(h3);

    // Desired area (e.g., normalize to an area of 500)
    double desired_area = 35*0.67; // 500.0;
    double current_area = h1->Integral();
    std::cout << "Current area: " << current_area << std::endl;
    double scale_factor = desired_area / current_area;
    // Rescale the histogram
    h1->Scale(scale_factor);
    h2->Scale(scale_factor);
    h3->Scale(scale_factor);

    // Observable (x-axis, e.g., energy)
    RooRealVar x("x", "Reconstructed Energy", 0, 4000);

    // Create RooDataHist from TH1 histograms
    RooDataHist data_signal("data_signal", "Data Signal", x, h_signal);
    RooDataHist data_background("data_background", "Data Background", x, h_background);

    // Create PDFs (Probability Distribution Functions) for signal and background
    RooHistPdf pdf_signal("pdf_signal", "Signal PDF", x, data_signal);
    RooHistPdf pdf_background("pdf_background", "Background PDF", x, data_background);

    // Define yield variables (number of signal and background events)
    RooRealVar n_signal("n_signal", "Signal events", 100, 0, 1000);
    RooRealVar n_background("n_background", "Background events", 500, 0, 2000);

    // Combine signal and background into a single model
    RooAddPdf model("model", "Signal + Background", RooArgList(pdf_signal, pdf_background), RooArgList(n_signal, n_background));

    // Generate data (for the purpose of this example, use combined model to generate data)
    RooDataHist data_combined("data_combined", "Combined Data", x, RooFit::Import(*h_signal));
    data_combined.add(RooDataHist("background_data", "Background Data", x, RooFit::Import(*h_background)));

    // Fit the model to the data
    RooFitResult* result = model.fitTo(data_combined);

    // Print the fit result
 //   result->Print();
    // Plot the data and fit result
    TCanvas *c1 = new TCanvas("c1", "Signal + Background Fit", 800, 600);
    RooPlot *xframe = x.frame();
    data_combined.plotOn(xframe);
    model.plotOn(xframe);
    model.plotOn(xframe, Components(pdf_background), LineStyle(kDashed), LineColor(kRed));  // Plot background
    model.plotOn(xframe, Components(pdf_signal), LineStyle(kDashed), LineColor(kBlue));     // Plot signal
    xframe->Draw();

    c1->Update();
    c1->SaveAs("plots/xframe.jpg");

    // Calculate the significance using the fitted signal and background yields
    double signal_yield = n_signal.getVal();
    double background_yield = n_background.getVal();
    
    // Using the approximation S/sqrt(B) for significance
    double significance = signal_yield / sqrt(background_yield);

    std::cout << "Signal yield: " << signal_yield << std::endl;
    std::cout << "Background yield: " << background_yield << std::endl;
    std::cout << "Significance: " << significance << std::endl;
}

TH2D* transform_migration(TH2D *h1, const char* histname, const char* title) {
    TH2D *mh1 = new TH2D(histname, title,
                         h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),  // X-axis (True Energy)
                         h1->GetNbinsY(), h1->GetYaxis()->GetXmin(), h1->GetYaxis()->GetXmax()); // Y-axis (Reconstructed Energy)

    // Loop over the bins of the original histogram h1
    for (int i = 1; i <= h1->GetNbinsX(); ++i)
    { // Loop over X-axis (True Energy)
        for (int j = 1; j <= h1->GetNbinsY(); ++j)
        { // Loop over Y-axis (Reconstructed Energy)

            // Get the content of the bin (i,j) in the original histogram
            double content = h1->GetBinContent(i, j);

            // Fill the new histogram with inverted Y-axis:
            // The true energy (X) remains the same, but we swap the reconstructed energy (Y) bins
            int inverted_j = h1->GetNbinsY() - j + 1;

            // Set the content in the new histogram mh1
            mh1->SetBinContent(i, inverted_j, content);
        }
    }
    // Change the Y-axis labels to reflect the inverted order
    for (int j = 1; j <= h1->GetNbinsY(); ++j)
    {
        // Get the center of the Y-axis bin in h1
        double originalLabel = h1->GetYaxis()->GetBinCenter(j);

        // Set the inverted label in mh1 for the corresponding bin
        int inverted_j = h1->GetNbinsY() - j + 1;
        mh1->GetYaxis()->SetBinLabel(inverted_j, Form("%.1f", originalLabel)); // Set the inverted label
    }
    // Change the X-axis labels so it looks similar
    for (int j = 1; j <= h1->GetNbinsY(); ++j)
    {
        // Get the center of the Y-axis bin in h1
        double originalLabel = h1->GetYaxis()->GetBinCenter(j);
        mh1->GetXaxis()->SetBinLabel(j, Form("%.1f", originalLabel)); // Set the inverted label
    }
    return mh1;
}

void nueselplot(int d, std::string what, std::string cut = "", std::string savefile = "plots/plot.png",
    double xmin = 0., double xmax = 2000., double ymin = 0., double ymax = 40., bool inverty = false, bool overlap = false) {
    TCanvas *c1 = new TCanvas("c1", "nueCC selection", 1024 , 800);
    gStyle->SetStatFontSize(0.05);  // Set the font size for the statistics box
    if(d==1) {
        TH1D *h1 = new TH1D("h1","nueCC", 10,xmin,xmax);
        TH1D *h2 = new TH1D("h2","numuCC", 10,xmin,xmax);
        TH1D *h3 = new TH1D("h3","NC", 10, xmin, xmax);
        event->Draw((what + std::string(">>h1")).c_str(), (cut + std::string("t_reaction == 1")).c_str());
        event->Draw((what + std::string(">>h2")).c_str(), (cut + std::string("t_reaction == 2")).c_str());
        event->Draw((what + std::string(">>h3")).c_str(), (cut + std::string("t_reaction > 3")).c_str());
        c1->Clear();
            h1->SetMarkerStyle(20);     // Set marker style to circle
            h1->SetMarkerSize(1.0);     // Set marker size
            h1->SetMarkerColor(kBlack); // Set marker color to black
            h1->GetXaxis()->SetTitle(what.c_str());
            h2->SetMarkerStyle(20);     // Set marker style to circle
            h2->SetMarkerSize(1.0);     // Set marker size
            h2->SetMarkerColor(kBlack); // Set marker color to black
            h2->GetXaxis()->SetTitle(what.c_str());
            h3->SetMarkerStyle(20);     // Set marker style to circle
            h3->SetMarkerSize(1.0);     // Set marker size
            h3->SetMarkerColor(kBlack); // Set marker color to black
            h3->GetXaxis()->SetTitle(what.c_str());
        if(!overlap){
            c1->Divide(2, 2);
            c1->cd(1);
            h1->Draw("ep");
            h1->Draw("hist,same");
            c1->cd(2);
            h2->Draw("ep");
            h2->Draw("hist,same");
            c1->cd(3);
            h3->Draw("ep");
            h3->Draw("hist,same");
            std::cout << what << " " << cut << std::endl;
            std::cout << " nueCC " << h1->Integral() << std::endl;
            std::cout << " numuCC " << h2->Integral() << std::endl;
            std::cout << " nuNC " << h3->Integral() << std::endl;
        }
        else
        {
            c1->cd();
            c1->SetLogy();
            h1->SetFillColor(kRed);
            h2->SetFillColor(kBlue);
            h3->SetFillColor(kGreen);
            THStack *hs = new THStack("hs","nue selection");
            hs->Add(h3);
            hs->Add(h2);
            hs->Add(h1);
            hs->Draw("hist");
            hs->GetXaxis()->SetTitle("Reconstructed energy");
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            legend->AddEntry(h1, "nueCC", "f");
            legend->AddEntry(h2, "numuCC", "f");
            legend->AddEntry(h3, "nuNC", "f");
            legend->Draw();
        }
    }
    else
    {
        TH2D *h1 = new TH2D("h1","nueCC", 10,xmin,xmax,10,ymin,ymax);
        TH2D *h2 = new TH2D("h2","numuCC", 10.,xmin,xmax,10,ymin,ymax);
        TH2D *h3 = new TH2D("h3","NC", 10, xmin, xmax,10,ymin,ymax);
        event->Draw((what + std::string(">>h1")).c_str(), (cut + std::string("t_reaction == 1")).c_str());
        event->Draw((what + std::string(">>h2")).c_str(), (cut + std::string("t_reaction == 2")).c_str());
        event->Draw((what + std::string(">>h3")).c_str(), (cut + std::string("t_reaction > 3")).c_str());
        if(inverty) {
            TH2D* mh1 = transform_migration(h1, "mh1", "nueCC");
            TH2D* mh2 = transform_migration(h2, "mh2", "numuCC");
            TH2D* mh3 = transform_migration(h3, "mh3", "nuNC");
            c1->Clear();
            c1->Divide(2, 2);
            c1->cd(1);
            mh1->SetStats(0);  // Disable the stats box for this histogram
            mh1->Draw("COLZ");
            mh1->Draw("TEXT,same");
            mh1->GetXaxis()->SetTitle(what.c_str());
            c1->cd(2);
            mh2->SetStats(0);  // Disable the stats box for this histogram
            mh2->Draw("COLZ");
            mh2->Draw("TEXT,same");
            mh2->GetXaxis()->SetTitle(what.c_str());
            c1->cd(3);
            mh3->SetStats(0);  // Disable the stats box for this histogram
            mh3->Draw("COLZ");
            mh3->Draw("TEXT,same");
            mh3->GetXaxis()->SetTitle(what.c_str());
        } else {
            c1->Clear();
            c1->Divide(2, 2);
            c1->cd(1);
            h1->Draw("COLZ,TEXT");
            h1->GetXaxis()->SetTitle(what.c_str());
            c1->cd(2);
            h2->Draw("COLZ,TEXT");
            h2->GetXaxis()->SetTitle(what.c_str());
            c1->cd(3);
            h3->Draw("COLZ,TEXT");
            h3->GetXaxis()->SetTitle(what.c_str());
        }
    }
    c1->Modified();
    c1->Update();
    c1->SaveAs(savefile.c_str());
}

int main() {
 //   gSystem->Load("libTPORec.so");

    TFile f("Analysis_1.root","READ");
    event = (TTree*)f.Get("Event");

    const char *enecompensation = "c_E1*1.1+rear_Cal*5.7";
    nueselplot(2,"t_Eneutrino:c_E1+rear_Cal","","plots/migrationEnu.png", 0., 2000., 0., 4000.);
    std::string p1 = std::string("(") + std::string(enecompensation) + std::string("-t_Eneutrino)/t_Eneutrino");
    nueselplot(1,p1.c_str(),"","plots/migrationEnu1.png", -2.,2.);
    std::string p2 = std::string("t_Eneutrino:") + std::string(enecompensation);
    nueselplot(2,p2.c_str(),"","plots/migrationEnu2.png", 0., 4000., 0., 4000.,true);
 
    nueselplot(2,"rear_MuCal:c_E1+rear_Cal","","plots/rearmucal.png");

    nueselplot(1,"c_E1", "", "plots/rearmucal0.png");
    nueselplot(1,"c_E1", "rear_MuCal<4 &&", "plots/rearmucal1.png");
    nueselplot(1,"c_E1+rear_Cal", "rear_MuCal<4 &&", "plots/rearmucal2.png");
    nueselplot(1,enecompensation, "rear_MuCal<4 &&", "plots/rearmucal3.png", 0., 4000.);
    nueselplot(1,enecompensation, "rear_MuCal<4 &&", "plots/rearmucal4.png", 0., 4000., 0., 0., false, true);

    nueselplot(1,"c_E2","rear_MuCal<4 &&","plots/cE2.png", 0.,10.);
    nueselplot(1,"c_chi2_1","rear_MuCal<4 &&","plots/cchi2.png", 0.,2000.);
    nueselplot(1,"c_a_1","rear_MuCal<4 &&","plots/ca.png", 0.,20.);

    nueselplot(1,enecompensation, "rear_MuCal<4 && c_a_1>12 && ", "plots/compiled1.png", 0., 4000., 0., 0., false, false);
    nueselplot(1,enecompensation, "rear_MuCal<4 && c_a_1>12 && ", "plots/compiled2.png", 0., 4000., 0., 0., false, true);

    signal_background_significance(enecompensation, "rear_MuCal<4 && c_a_1>12 &&");

}