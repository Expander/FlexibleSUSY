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

#include "plot.h"

void plotMSSM(const TString& file_name = "higgs-study/data/scanMSSM.dat")
{
   std::ifstream ifs(file_name.Data());
   std::string line;

   if (!ifs.good()) {
      cout << "Error: could not open file: " << file_name << endl;
      return;
   }

   const int nbinsx = 30, nbinsy = 30;
   const double xlow = 0., xhigh = 50.;
   const double ylow = 0., yhigh = 10.;

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

         m0 /= 1000.; // convet to TeV

         h->Fill(tanb, m0, mh, 1.0);
      }
   }

   h->SetTitle("CMSSM $m_{h}$ / GeV");
   h->GetXaxis()->SetTitle("$\tan\beta$");
   h->GetYaxis()->SetTitle("$m_0$ / TeV");
   SetZminZmax(h);

   TCanvas* canvas = new TCanvas("canvas", "CMSSM Higgs mass", 800, 600);
   canvas->cd(1);
   h->Draw("colz");

   TString tex_file(file_name);
   tex_file.ReplaceAll(".dat",".tex");
   canvas->Print(tex_file);
}
