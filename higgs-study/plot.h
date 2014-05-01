
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
