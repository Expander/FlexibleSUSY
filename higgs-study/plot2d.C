#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>

#include "TH1.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"

#include "plot.h"

void set_plot_style()
{
    static const Int_t NRGBs = 5;
    static const Int_t NCont = 30;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void plot2d(const TString& file_name = "higgs-study/data/scan_MSSM.dat",
            const TString& title = "m_{h}^{\\text{pole}}\\text{ / GeV}")
{
   std::ifstream ifs(file_name.Data());
   std::string line;

   if (!ifs.good()) {
      cout << "Error: could not open file: " << file_name << endl;
      return;
   }

   const int nbinsx = 50, nbinsy = 50;
   const double xlow = 0., xhigh = 50.;
   const double ylow = 0., yhigh = 10.; // in TeV

   TProfile2D* h = new TProfile2D("h", file_name, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
   h->SetStats(0);

   while (getline(ifs,line)) {
      std::istringstream input(line);
      std::string word;
      input >> word;

      if (word.find("#") == std::string::npos) {
         std::istringstream kk(line);
         double tanb, m0, mh, error;
         bool stream_ok = true;

         kk >> tanb >> m0 >> mh >> error;

         stream_ok = kk.good();
         if (!stream_ok) {
            cout << "Error: invalid line: " << line << endl;
            continue;
         }

         if (error)
            continue;

         m0 /= 1000.; // convert to TeV

         h->Fill(tanb, m0, mh, 1.0);
      }
   }

   h->SetTitle(title);
   h->GetXaxis()->SetTitle("\\tan\\beta");
   h->GetYaxis()->SetTitle("$m_0$ / TeV");
   h->SetMinimum(95.);
   h->SetMaximum(135.);
   // SetZminZmax(h);

   TCanvas* canvas = new TCanvas("canvas", title, 400, 300);
   canvas->cd(1);
   // gStyle->SetPaperSize(10.,10.); // in cm
   set_plot_style();
   h->Draw("colz");

   TString tex_file(file_name);
   tex_file.ReplaceAll(".dat",".tex");
   TString pdf_file(file_name);
   pdf_file.ReplaceAll(".dat",".pdf");
   canvas->Print(tex_file);

   gSystem->Exec(TString("mv ") + tex_file + " img.tex");
   gSystem->Exec(TString("pdflatex ") + " img_container.tex");
   gSystem->Exec(TString("mv img_container.pdf ") + pdf_file);
}
