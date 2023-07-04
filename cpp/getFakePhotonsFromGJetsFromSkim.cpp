#include "TFile.h"
#include "TChain.h"

int getFakePhotonsFromGJetsFromSkim() {
  TFile *fout = new TFile("fakePhotonsFromGJetsFromSkim.root","RECREATE");
  TChain *chain = new TChain("Events");

  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD/2018/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v0/tree_*");
  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD/2018/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v0/tree_*");
  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD/2018/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v0/tree_*");
  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD/2018/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v0/tree_*");
  chain->Add("/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD/2018/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v0/tree_*");

  chain->Draw("Photon_mvaID","Photon_genPartFlav!=1","");
  TH1F *hist = (TH1F*)gPad->GetPrimitive("htemp"); 

  hist->Write();
  fout->Write();
  fout->Close();

  return 0;
}
