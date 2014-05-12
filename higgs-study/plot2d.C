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
#include "TROOT.h"
#include "TGraph.h"

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

bool has_empty_neighbor(const TProfile2D* h, Double_t x, Double_t y)
{
   Int_t bin = h->FindFixBin(x, y);
   Int_t binx, biny, binz;
   h->GetBinXYZ(bin, binx, biny, binz);

   if (h->GetBinContent(binx + 1, biny) == 0.)
      return true;
   if (h->GetBinContent(binx - 1, biny) == 0.)
      return true;
   if (h->GetBinContent(binx, biny + 1) == 0.)
      return true;
   if (h->GetBinContent(binx, biny - 1) == 0.)
      return true;

   return false;
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

   TCanvas* canvas = new TCanvas("canvas", title, 400, 300);
   canvas->cd(1);
   // gStyle->SetPaperSize(10.,10.); // in cm

   // define contours
   Double_t contours[1];
   contours[0] = 125.9;
   h->SetContour(1, contours);
   // Draw contours as filled regions, and Save points
   h->Draw("CONT Z LIST");
   canvas->Update(); // Needed to force the plotting and retrieve the contours in TGraphs

   // Get Contours
   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   TList* contLevel = NULL;
   TGraph* curv     = NULL;
   TGraph* gc       = NULL;
   Int_t TotalConts = 0, nGraphs = 0;
   if (conts == NULL){
      printf("*** No Contours Were Extracted!\n");
      TotalConts = 0;
      return;
   } else {
      TotalConts = conts->GetSize();
   }
   printf("TotalConts = %d\n", TotalConts);
   for(int i = 0; i < TotalConts; i++){
      contLevel = (TList*)conts->At(i);
      printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
      nGraphs += contLevel->GetSize();
   }

   // draw the histogram
   h->SetContour(20); // set back to default
   h->SetMinimum(95.);
   h->SetMaximum(135.);
   // SetZminZmax(h);
   set_plot_style();
   h->Draw("colz");

   // draw contours
   for(int i = 0; i < TotalConts; i++){
      contLevel = (TList*)conts->At(i);
      // Get first graph from list on curves on this level
      curv = (TGraph*)contLevel->First();
      for(int j = 0; j < contLevel->GetSize(); j++){
         Double_t x0, y0;
         curv->GetPoint(1, x0, y0);
         // check if neighbor bins in h are empty
         bool empty_neighbor = has_empty_neighbor(h, x0, y0);

         if (empty_neighbor) {
            curv = (TGraph*)contLevel->After(curv); // Get Next graph
            continue;
         }

         gc = (TGraph*)curv->Clone();
         gc->Draw("C same");

         curv = (TGraph*)contLevel->After(curv); // Get Next graph
      }
   }
   canvas->Update();
   printf("\n\n\tExtracted %d Contours and %d Graphs \n", TotalConts, nGraphs );

   TString tex_file(file_name);
   tex_file.ReplaceAll(".dat",".tex");
   TString pdf_file(file_name);
   pdf_file.ReplaceAll(".dat",".pdf");
   canvas->Print(tex_file);

   gSystem->Exec(TString("mv ") + tex_file + " img.tex");
   gSystem->Exec(TString("pdflatex ") + " img_container.tex");
   gSystem->Exec(TString("mv img_container.pdf ") + pdf_file);
}
