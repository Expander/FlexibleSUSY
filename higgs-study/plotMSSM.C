#include <iostream>
#include <fstream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>

#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

void SetZminZmax(TProfile2D* hist)
{
   Double_t min = DBL_MAX;
   Double_t max = DBL_MIN;

   for (Int_t x = 1; x <= hist->GetNbinsX(); ++x) {
      for (Int_t y = 1; y <= hist->GetNbinsY(); ++y) {
         const Int_t bin = hist->GetBin(x, y);
         const Double_t entries = hist->GetBinEntries(bin);
         if (entries > 0) {
            const Double_t content = hist->GetBinContent(bin);
            if (content < min)
               min = content;
            if (content > max)
               max = content;
         }
      }
   }
   min -= 0.01 * TMath::Abs(min);
   max += 0.01 * TMath::Abs(max);

   hist->SetMinimum(static_cast<Double_t>(min));
   hist->SetMaximum(static_cast<Double_t>(max));
}

void plotMSSM(const TString& file_name = "higgs-study/scanMSSM.dat")
{
   std::ifstream ifs(file_name.Data());
   std::string line;

   if (!ifs.good()) {
      cout << "Error: could not open file: " << file_name << endl;
      return;
   }

   const int nbinsx = 10, nbinsy = 10;
   const double xlow = 0., xhigh = 100.;
   const double ylow = 0., yhigh = 500.;

   TProfile2D* h = new TProfile2D("h", file_name, nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
   h->SetStats(0);

   while (getline(ifs,line)) {
      std::istringstream input(line);
      std::string word;
      input >> word;

      if (word.find("#") == std::string::npos) {
         std::istringstream kk(line);
         double tanb, a0, mh, error;
         bool stream_ok = true;

         kk >> tanb >> a0 >> mh >> error;

         stream_ok = kk.good();
         if (!stream_ok) {
            cout << "Error: invalid line: " << line << endl;
            continue;
         }

         if (error)
            continue;

         h->Fill(tanb, a0, mh, 1.0);
      }
   }

   h->SetTitle("CMSSM $m_{h}$ / GeV");
   h->GetXaxis()->SetTitle("tan(beta)");
   h->GetYaxis()->SetTitle("a0 / GeV");
   SetZminZmax(h);

   TCanvas* canvas = new TCanvas("canvas", "CMSSM Higgs mass", 800, 600);
   canvas->cd(1);
   h->Draw("colz");

   TString tex_file(file_name.ReplaceAll(".dat",".tex"));
   canvas->Print(tex_file);
}
