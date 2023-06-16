#include "TFile.h"
#include "TChain.h"

int getFakePhotonsFromGJetsFromPresel() {
  TFile *fout = new TFile("fakePhotonsFromGJetsFromPresel.root","RECREATE");
  TChain *chain = new TChain("tout");

  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/preselectedNanoAOD/fakePhotonsInPresel/output_GJets_HT-*To*_*.root");

  chain->Draw("LeadPhoton_mvaID>>h1","LeadPhoton_genPartFlav!=1","");
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  chain->Draw("SubleadPhoton_mvaID>>h2","SubleadPhoton_genPartFlav!=1","");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2");

  TH1F *hist = (TH1F*)h1->Clone("htemp");
  hist->Add(h2);

  hist->Write();
  fout->Write();
  fout->Close();

  return 0;
}
