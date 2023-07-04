#include "TFile.h"
#include "TChain.h"

int getFakePhotonsFromGJetsFromPresel(TString year="", unsigned int incl_barrel_endcap=0) {
  if (year != "") year = "_"+year;

  TString etaLabel, etaCutLead, etaCutSublead;
  if (incl_barrel_endcap==0) { // Inclusive in eta
    etaLabel = "";
    etaCutLead = "";
    etaCutSublead = "";
  }
  else if (incl_barrel_endcap==1) { // Barrel
    etaLabel = "_barrel";
    etaCutLead = "&& fabs(LeadPhoton_eta)<1.5";
    etaCutSublead = "&& fabs(SubleadPhoton_eta)<1.5";
  }
  else if (incl_barrel_endcap==2) { // Endcap
    etaLabel = "_endcap";
    etaCutLead = "&& fabs(LeadPhoton_eta)>1.5";
    etaCutSublead = "&& fabs(SubleadPhoton_eta)>1.5";
  }
  else {
    std::cout << "Invalid option for incl_barrel_endcap!\n";
    return 1;
  }

  TFile *fout = new TFile("fakePhotonsFromGJetsFromPresel"+year+etaLabel+".root","RECREATE");
  TChain *chain = new TChain("tout");

  if (year == "") year = "_*";
  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/preselectedNanoAOD/fakePhotonsInPresel/output_GJets_HT-*To*"+year+".root");

  for ( auto f : *chain->GetListOfFiles() ) std::cout<<f->GetTitle()<<"\n";
  chain->Draw("LeadPhoton_mvaID>>h1","LeadPhoton_genPartFlav!=1 "+etaCutLead,"");
  TH1F *h1 = (TH1F*)gDirectory->Get("h1");
  chain->Draw("SubleadPhoton_mvaID>>h2","SubleadPhoton_genPartFlav!=1 "+etaCutSublead,"");
  TH1F *h2 = (TH1F*)gDirectory->Get("h2");

  TH1F *hist = (TH1F*)h1->Clone("htemp");
  hist->Add(h2);

  hist->Write();
  fout->Write();
  fout->Close();

  return 0;
}
