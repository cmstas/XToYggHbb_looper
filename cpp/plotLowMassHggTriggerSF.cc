#include "../NanoCORE/Tools/lowMassHggTriggerSF.h"
using namespace lowMassHggTriggerSF;

int plotLowMassHggTriggerSF() {
  gStyle->SetOptStat(0);
  set_lowMassHggTriggerSF();

  std::vector<TString> years = { "2018", "2017", "2016", "2016EBHiR9" };
  std::vector<TString> photons = { "Lead", "Sublead" };
  std::vector<float> etas = { 0.0/*barrel*/, 2.0/*endcap*/ };

  for ( auto &year : years ){
    bool EBHiR9 = year.Contains("EBHiR9");

    for ( const auto &photon : photons ) {
      set_etabins();

      for ( const auto &eta : etas ) {
        if ( EBHiR9 ) year = "2016";

        set_r9bins(eta, photon, EBHiR9, year);

        std::vector<float> r9Vec;
        for ( unsigned int ir9=0;  ir9<nr9bins; ir9++ ) {
          if ( thresholds_r9[ir9]==999.0 ) {
            r9Vec.push_back(1); // Dummy maximum for R9
            break;
          }
          r9Vec.push_back(thresholds_r9[ir9]);
        }

        const unsigned int r9Len = r9Vec.size();
        float *r9Arr = &r9Vec[0];

        for ( unsigned int ir9=1; ir9<r9Len; ir9++ ) {
          float r9 = (r9Arr[ir9-1] + r9Arr[ir9]) / 2; // Center of R9 bin
          set_ptbins(eta, r9, photon, EBHiR9, year);

          std::vector<float> ptVec;
          for ( unsigned int ipt=0; ipt<nptbins; ipt++ ) {
            if ( thresholds_pt[ipt]==999999.0 ) {
              ptVec.push_back(100); // Dummy maximum for pt
              break;
            }
            ptVec.push_back(thresholds_pt[ipt]);
          }

          const unsigned int ptLen = ptVec.size();
          float *ptArr = &ptVec[0];

          if ( EBHiR9 ) year = "2016EBHiR9";
          TString name = "lowMassTrigSF_"+year+"_"+photon+"_"+(eta < 1.0 ? "barrel" : "endcap")+"_r9"+to_string(r9Arr[ir9-1])+"-"+to_string(r9Arr[ir9]);
          TCanvas *c = new TCanvas(name);
          TH1F *h = new TH1F(name, name, ptLen-1, ptArr);
          for ( unsigned int ipt=1; ipt<ptLen; ipt++ ) {
            float pt = (ptArr[ipt-1] + ptArr[ipt]) / 2; // Center of pt bin
            TString etabin = get_etaBin(eta);
            TString r9bin = get_r9Bin(r9);
            TString ptbin = get_ptBin(pt);
            h->SetBinContent(ipt, photon == "Lead" ? idsf_lead[year][etabin][r9bin][ptbin] : idsf_sublead[year][etabin][r9bin][ptbin]);
            h->SetBinError(ipt, photon == "Lead" ? idsfunc_lead[year][etabin][r9bin][ptbin] : idsfunc_sublead[year][etabin][r9bin][ptbin]);
          }

          h->Draw("E1");
          h->SetMinimum(0.0);
          h->SetMaximum(1.2);
          h->GetXaxis()->SetTitle("p_{T} [GeV]");
          h->GetYaxis()->SetTitle("Efficiency");
          c->SaveAs(".png");
        }
      }
    }
  }

  return 0;
}
