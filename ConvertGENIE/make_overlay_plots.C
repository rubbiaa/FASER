// make_overlay_plots.C -- overlays the per-detector reference histograms
// (h_z_all, h_x_all, h_y_all, h_xy_all, h_xz_all, h_yz_all) ConvertGENIE.cc
// writes into each detector's own interaction_plots_<detector>.root, so
// e.g. all three of 3DCAL/ECAL/AHCAL's Z-vertex distributions end up on
// one canvas/PNG instead of three separate ones.
//
// Invoked by run_convertgenie.py (never by hand, normally) as:
//   root -l -b -q 'make_overlay_plots.C("<manifest>", "<outdir>")'
// where <manifest> is a small text file written by run_convertgenie.py,
// one line per detector that actually produced output:
//   <detector_name> <path to that detector's interaction_plots_X.root>
//
// Why a manifest file rather than passing the detector/file lists directly
// as macro arguments: ROOT's CLI argument parsing for a function invoked
// as `root 'macro.C(args)'` handles strings/ints fine but is finicky and
// version-dependent for anything array-like (std::vector<std::string>,
// C arrays of TString, ...) -- a one-line-per-entry manifest file sidesteps
// that entirely and is trivial to generate from Python and to inspect by
// hand if something goes wrong.
//
// A 1D histogram (h_z_all/h_x_all/h_y_all) overlay is a straightforward
// same-canvas HIST/HIST SAME with each detector in its own line color and
// a legend. A 2D histogram doesn't overlay the same way -- stacking
// several transparent COLZ color maps on one pad isn't legible -- so
// instead this divides one canvas into one pad per detector and draws
// each with COLZ, the same style ConvertGENIE.cc itself uses for its own
// single-detector PNGs; each pad's title (now "<detector> (L=... fb^-1)",
// set by ConvertGENIE.cc) identifies which detector that pad is.

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

void make_overlay_plots(const char *manifest_path, const char *outdir)
{
  std::vector<std::string> detector_names;
  std::vector<std::string> file_paths;

  std::ifstream in(manifest_path);
  if (!in.is_open())
  {
    std::cerr << "make_overlay_plots: could not open manifest " << manifest_path << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty())
      continue;
    std::istringstream iss(line);
    std::string det, path;
    iss >> det >> path;
    if (det.empty() || path.empty())
      continue;
    detector_names.push_back(det);
    file_paths.push_back(path);
  }
  in.close();

  if (detector_names.empty())
  {
    std::cerr << "make_overlay_plots: no usable entries in manifest " << manifest_path << std::endl;
    return;
  }

  std::vector<TFile *> files;
  std::vector<std::string> opened_detector_names;
  for (size_t i = 0; i < file_paths.size(); i++)
  {
    TFile *f = TFile::Open(file_paths[i].c_str(), "READ");
    if (!f || f->IsZombie())
    {
      std::cerr << "make_overlay_plots: could not open " << file_paths[i]
                << " for detector " << detector_names[i] << " -- skipping it" << std::endl;
      continue;
    }
    files.push_back(f);
    opened_detector_names.push_back(detector_names[i]);
  }

  if (files.empty())
  {
    std::cerr << "make_overlay_plots: none of the manifest's ROOT files could be opened" << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir(outdir, true);

  int colors[5] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1};

  // 1D overlays: h_z_all/z_distribution, h_x_all/x_distribution, h_y_all/y_distribution.
  const char *hist1d_names[3] = {"h_z_all", "h_x_all", "h_y_all"};
  const char *plot1d_names[3] = {"z_distribution", "x_distribution", "y_distribution"};
  const char *axis1d_labels[3] = {"Z [mm];Events", "X [mm];Events", "Y [mm];Events"};

  for (int ih = 0; ih < 3; ih++)
  {
    TCanvas c(Form("c_overlay_%s", plot1d_names[ih]), "ConvertGENIE overlay", 800, 600);
    TLegend leg(0.70, 0.70, 0.89, 0.89);
    leg.SetBorderSize(0);

    std::vector<TH1 *> hs;
    double ymax = 0;
    for (size_t i = 0; i < files.size(); i++)
    {
      TH1 *h = dynamic_cast<TH1 *>(files[i]->Get(hist1d_names[ih]));
      if (!h)
      {
        std::cerr << "make_overlay_plots: " << hist1d_names[ih] << " not found in "
                  << file_paths[i] << std::endl;
        continue;
      }
      h->SetDirectory(0); // detach from its (about-to-close) file
      h->SetLineColor(colors[i % 5]);
      h->SetLineWidth(2);
      hs.push_back(h);
      leg.AddEntry(h, opened_detector_names[i].c_str(), "l");
      if (h->GetMaximum() > ymax)
        ymax = h->GetMaximum();
    }
    if (hs.empty())
      continue;

    // hs[0]'s own title (e.g. "3DCAL (L=1000 fb^-1)") is correct for its
    // single-detector PNG but would misrepresent this multi-detector
    // overlay -- replace just the in-memory (detached, not written back)
    // copy's title with a generic one; the legend is what actually
    // identifies each detector's curve here.
    hs[0]->SetTitle((std::string(plot1d_names[ih]) + " -- detector overlay;" + axis1d_labels[ih]).c_str());
    hs[0]->SetMaximum(ymax * 1.2);
    hs[0]->Draw("HIST");
    for (size_t i = 1; i < hs.size(); i++)
      hs[i]->Draw("HIST SAME");
    leg.Draw();
    c.SaveAs((std::string(outdir) + "/" + plot1d_names[ih] + "_overlay.png").c_str());
  }

  // 2D "overlays": one pad per detector, each drawn with COLZ (same style
  // as ConvertGENIE.cc's own single-detector PNGs) -- see the header
  // comment for why this isn't a literal single-pad overlay like the 1D
  // case above.
  const char *hist2d_names[3] = {"h_xy_all", "h_xz_all", "h_yz_all"};
  const char *plot2d_names[3] = {"xy_distribution", "xz_distribution", "yz_distribution"};

  for (int ih = 0; ih < 3; ih++)
  {
    std::vector<TH2 *> hs;
    for (size_t i = 0; i < files.size(); i++)
    {
      TH2 *h = dynamic_cast<TH2 *>(files[i]->Get(hist2d_names[ih]));
      if (!h)
      {
        std::cerr << "make_overlay_plots: " << hist2d_names[ih] << " not found in "
                  << file_paths[i] << std::endl;
        continue;
      }
      h->SetDirectory(0);
      hs.push_back(h);
    }
    if (hs.empty())
      continue;

    TCanvas c(Form("c_overlay_%s", plot2d_names[ih]), "ConvertGENIE overlay",
              400 * (int)hs.size(), 600);
    c.Divide((int)hs.size(), 1);
    for (size_t i = 0; i < hs.size(); i++)
    {
      c.cd((int)i + 1);
      hs[i]->Draw("COLZ");
    }
    c.cd(0);
    c.SaveAs((std::string(outdir) + "/" + plot2d_names[ih] + "_overlay.png").c_str());
  }

  for (size_t i = 0; i < files.size(); i++)
    files[i]->Close();

  std::cout << "make_overlay_plots: wrote overlay PNGs for " << files.size()
            << " detector(s) to " << outdir << std::endl;
}
