#include "TROOT.h"
#include "ScanChain_Hgg.C"
#include <fstream>

unsigned char2unsigned(const char *c) {
  char p = *c;
  unsigned res = 0;
  while (p) {
      res = res*10 + (p - '0');
      c++;
      p = *c;
  }
  return res;
}
int char2int(const char *c) {
  return (*c == '-') ? -char2unsigned(c+1) : char2unsigned(c);
}

double getSumOfGenEventSumw(TChain *chaux, bool isMC)
{
  double genEventSumw, sumOfGenEventSumw=0.0;
  if (isMC) {
    chaux->SetBranchAddress("genEventSumw",&genEventSumw);
    for (unsigned int run = 0; run < chaux->GetEntriesFast(); run++) {
      chaux->GetEntry(run);
      sumOfGenEventSumw += genEventSumw;
    }
  }
  else {
    sumOfGenEventSumw=1.0;
  }
  return sumOfGenEventSumw;
}


int main(int argc, char **argv) {
  //Arguments
  const char* outdir      = ( argc > 1  ? argv[1]            : "temp_data" );
  TString yearArg         = ( argc > 2  ? argv[2]            : "all" );
  int run_data            = ( argc > 3  ? char2int(argv[3])  : 1 );
  int run_MCbkg           = ( argc > 4  ? char2int(argv[4])  : 1 );
  int run_signal          = ( argc > 5  ? char2int(argv[5])  : 1 );
  TString sampleArg       = ( argc > 6  ? argv[6]            : "all" );
  int lowMassMode         = ( argc > 7  ? char2int(argv[7])  : 1 );
  // Set to 0 to disable, set to 1 for central value, set to +/-2 to get uncertainty
  int prefireWeight       = ( argc > 8  ? char2int(argv[8])  : 1 );
  int PUWeight            = ( argc > 9  ? char2int(argv[9])  : 1 );
  int electronVetoSF      = ( argc > 10 ? char2int(argv[10]) : 1 );
  int triggerSF           = ( argc > 11 ? char2int(argv[11]) : 1 );
  int preselSF            = ( argc > 12 ? char2int(argv[12]) : 0 );
  int phoMVAIDWP90SF      = ( argc > 13 ? char2int(argv[13]) : 1 );
  int bTagSF              = ( argc > 14 ? char2int(argv[14]) : 1 ); // Set to +/-X to get uncertainty, X in [2,10]
  int fnufUnc             = ( argc > 15 ? char2int(argv[15]) : 0 ); // No central value
  int materialUnc         = ( argc > 16 ? char2int(argv[16]) : 0 ); // No central value
  int PhoScaleUnc         = ( argc > 17 ? char2int(argv[17]) : 0 ); // No central value
  int PhoSmearUnc         = ( argc > 18 ? char2int(argv[18]) : 0 ); // No central value
  int JESUnc              = ( argc > 19 ? char2int(argv[19]) : 0 ); // Set to +/-X to get uncertainty, X in [2,12]
  int JERUnc              = ( argc > 20 ? char2int(argv[20]) : 0 ); // No central value
  int HEMCheck            = ( argc > 21 ? char2int(argv[21]) : 0 ); // Set to 1 to check
  // Option to reproduce the summary.json
  int onlyCreateJSON      = ( argc > 22 ? char2int(argv[22]) : 0 );
  

  // Map definitions
  vector<TString> years = { };
  if ( yearArg=="all" ) {
    years.push_back("2018");
    years.push_back("2017");
    years.push_back("2016APV");
    years.push_back("2016nonAPV");
  }
  else if ( yearArg=="2018" || yearArg=="2017" || yearArg=="2016APV" || yearArg=="2016nonAPV" ) {
    years.push_back(yearArg);
  }
  else {
    cout << "\nInvalid option for year, exiting...\n\n";
    return 1;
  }
  vector<TString> samples = { };
  map<TString,TString> sample_names = { };
  map<TString,int> sample_procids = { };
  map<TString,map<TString,vector<TString>>> sample_prod = { };


  // Sample list: Data
  if (run_data) {
    if ( sampleArg=="Data" || sampleArg=="all" ) {
      TString sampleName = "Data";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 0});
      sample_names.insert({sampleName, sampleName});
      sample_prod.insert({sampleName, { { "2018",       { "Run2018A-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018B-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018C-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018D-UL2018_MiniAODv2_GT36-v2" } },
                                        { "2017",       { "Run2017B-UL2017_MiniAODv2-v1",
                                                          "Run2017C-UL2017_MiniAODv2-v2",
                                                          "Run2017D-UL2017_MiniAODv2-v1",
                                                          "Run2017E-UL2017_MiniAODv2-v1",
                                                          "Run2017F-UL2017_MiniAODv2-v2" } },
                                        { "2016APV",    { "Run2016B-ver1_HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016B-ver2_HIPM_UL2016_MiniAODv2-v3",
                                                          "Run2016C-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016D-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016E-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016F-HIPM_UL2016_MiniAODv2-v1" } },
                                        { "2016nonAPV", { "Run2016F-UL2016_MiniAODv2-v1",
                                                          "Run2016G-UL2016_MiniAODv2-v1",
                                                          "Run2016H-UL2016_MiniAODv2-v1" } } } });
    }
  }


  // Sample list: MC
  if (run_MCbkg) {
    if ( sampleArg=="DiPhoton" || sampleArg=="all" ) {
      TString sampleName = "DiPhoton";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 1});
      sample_names.insert({sampleName, "DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });
    }

    if ( sampleArg=="TT" || sampleArg=="all" ) {
      TString sampleName = "TTGG";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 2});
      sample_names.insert({sampleName, "TTGG_0Jets_TuneCP5_13TeV-amcatnlo-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });

      sampleName = "TTGJets";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 3});
      sample_names.insert({sampleName, "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });

      sampleName = "TTJets";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 4});
      sample_names.insert({sampleName, "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });
    }

    if ( sampleArg=="H" || sampleArg=="all" ) {
      TString sampleName = "VH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 5});
      sample_names.insert({sampleName, "VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });  

      sampleName = "VBFH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 6});
      sample_names.insert({sampleName,"VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });

      sampleName = "ttH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 7});
      sample_names.insert({sampleName, "ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v3" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });

      sampleName = "ggHToDiPhoM125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 8});
      sample_names.insert({sampleName, "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });
    }

    if ( sampleArg=="GJets" || sampleArg=="all" ) {
      vector<TString> M = { "40", "100", "200", "400", "600", "Inf" };
      map<TString,int> processId_M = { {"40",0}, {"100",1}, {"200",2}, {"400",3}, {"600",4} };
      for ( unsigned int imass=0; imass<M.size()-1; imass++ ) {
        TString sampleName = "GJets_HT-"+M[imass]+"To"+M[imass+1];
        samples.push_back(sampleName);
        sample_procids.insert({sampleName, 9+processId_M[M[imass]]});
        sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-madgraphMLM-pythia8"});
        sample_prod.insert({sampleName, { { "2018",       { ( M[imass] =="100" ?
                                                              "RunIISummer20UL18MiniAODv2-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2" :
                                                              "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" ) } },
                                          { "2017",       { ( M[imass] =="100" ?
                                                              "RunIISummer20UL17MiniAODv2-4cores5k_106X_mc2017_realistic_v9-v2" :
                                                              "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" ) } },
                                          { "2016APV",    { ( M[imass] =="100" ?
                                                              "RunIISummer20UL16MiniAODAPVv2-4cores5k_106X_mcRun2_asymptotic_preVFP_v11-v1" :
                                                              "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" ) } },
                                          { "2016nonAPV", { ( M[imass] =="100" ?
                                                              "RunIISummer20UL16MiniAODv2-4cores5k_106X_mcRun2_asymptotic_v17-v1" :
                                                              "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" ) } } } });
      }
    }

    if ( sampleArg=="DY" || sampleArg=="all" ) {
      TString sampleName = "DY";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 14});
      sample_names.insert({sampleName, "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });
    }

    if ( sampleArg=="VG" || sampleArg=="all" ) {
      vector<TString> V = { "W", "Z" };
      map<TString,int> processId_V = { {"W",0}, {"Z",1} };
      for ( unsigned int iV=0; iV<V.size(); iV++ ) {
        TString sampleName = V[iV]+"G";
        samples.push_back(sampleName);
        sample_procids.insert({sampleName, 15+processId_V[V[iV]]});
        sample_names.insert({sampleName, ( V[iV]=="W" ?
                                           sampleName+"ToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8" :
                                           sampleName+"ToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8" ) });
        sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                          { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                          { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                          { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });
      }
    }

    if ( sampleArg=="HHbbgg" || sampleArg=="all" ) {
      TString sampleName = "HHbbgg";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 17});
      sample_names.insert({sampleName, "GluGluToHHTo2B2G_node_cHHH1_TuneCP5_13TeV-powheg-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "privateUL18" } },
                                        { "2017",       { "privateUL17" } },
                                        { "2016APV",    { "privateUL16APV" } },
                                        { "2016nonAPV", { "privateUL16" } } } });
    }

    if ( sampleArg=="DiPhotonLow" || sampleArg=="all" ) {
      TString sampleName = "DiPhotonLow";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 18});
      sample_names.insert({sampleName, "DiPhotonJetsBox_MGG-40to80_13TeV-sherpa"});
      sample_prod.insert({sampleName, { { "2018",       { "privateUL18" } },
                                        { "2017",       { "privateUL17" } },
                                        { "2016APV",    { "privateUL16APV" } },
                                        { "2016nonAPV", { "privateUL16" } } } });
    }

    if ( sampleArg=="VV" || sampleArg=="all" ) {
      vector<TString> V = { "WW", "WZ", "ZZ" };
      map<TString,int> processId_V = { {"WW",0}, {"WZ",1}, {"ZZ",2} };
      for ( unsigned int iV=0; iV<V.size(); iV++ ) {
        TString sampleName = V[iV];
        samples.push_back(sampleName);
        sample_procids.insert({sampleName, 19+processId_V[V[iV]]});
        sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-pythia8"});
        sample_prod.insert({sampleName, { { "2018",       { ( V[iV] == "ZZ" ?
                                                              "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" :
                                                              "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" ) } },
                                          { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                          { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                          { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });
      }
    }

    if ( sampleArg=="QCD" || sampleArg=="all" ) {
      TString sampleName = "QCD_Pt-30to40_MGG-80toInf";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 22});
      sample_names.insert({sampleName, "QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });

      sampleName = "QCD_Pt-30toInf_MGG-40to80";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 23});
      sample_names.insert({sampleName, "QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });

      sampleName = "QCD_Pt-40toInf_MGG-80toInf";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 24});
      sample_names.insert({sampleName, "QCD_Pt-40ToInf_DoubleEMEnriched_MGG-80ToInf_TuneCP5_13TeV-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });
    }

    if ( sampleArg=="DDQCDGJets" || sampleArg=="all" ) {
      TString sampleName = "DDQCDGJets";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 25});
      sample_names.insert({sampleName, sampleName});
      sample_prod.insert({sampleName, { { "2018",       { "Run2018A-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018B-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018C-UL2018_MiniAODv2_GT36-v1",
                                                          "Run2018D-UL2018_MiniAODv2-v2" } },
                                        { "2017",       { "Run2017B-UL2017_MiniAODv2-v1",
                                                          "Run2017C-UL2017_MiniAODv2-v2",
                                                          "Run2017D-UL2017_MiniAODv2-v1",
                                                          "Run2017E-UL2017_MiniAODv2-v1",
                                                          "Run2017F-UL2017_MiniAODv2-v2" } },
                                        { "2016APV",    { "Run2016B-ver1_HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016B-ver2_HIPM_UL2016_MiniAODv2-v3",
                                                          "Run2016C-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016D-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016E-HIPM_UL2016_MiniAODv2-v1",
                                                          "Run2016F-HIPM_UL2016_MiniAODv2-v1" } },
                                        { "2016nonAPV", { "Run2016F-UL2016_MiniAODv2-v1",
                                                          "Run2016G-UL2016_MiniAODv2-v1",
                                                          "Run2016H-UL2016_MiniAODv2-v1" } } } });
    }

    if ( sampleArg=="TXX" || sampleArg=="all" ) {
      TString sampleName = "tZq";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 26});
      sample_names.insert({sampleName, "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });

      sampleName = "TTZ";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 27});
      sample_names.insert({sampleName, "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v1" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });

      sampleName = "TTW";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 28});
      sample_names.insert({sampleName, "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v1" } },
                                        { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                        { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1" } } } });
    }
  }

  // Sample list: Signal
  if (run_signal) {
    map<TString,TString> processId_MX = { {"240","1"},  {"280","2"},  {"300","3"},  {"320","4"},  {"360","5"},
                                          {"400","6"},  {"450","7"},  {"500","8"},  {"550","9"},  {"600","10"},
                                          {"650","11"}, {"700","12"}, {"750","13"}, {"800","14"}, {"850","15"},
                                          {"900","16"}, {"950","17"}, {"1000","18"} };
    map<TString,TString> processId_MY = { {"70","11"},  {"80","12"},  {"90","13"},  {"100","14"}, {"125","15"},
                                          {"150","16"}, {"170","17"}, {"190","18"}, {"250","19"}, {"300","20"},
                                          {"350","21"}, {"400","22"}, {"450","23"}, {"500","24"}, {"550","25"},
                                          {"600","26"}, {"650","27"}, {"700","28"}, {"800","29"} };

    if ( sampleArg.Contains("all") || sampleArg.Contains("nominal") ) {
      map<TString,vector<TString>> MComb = { { "240",  { "70","80","90","100" } },
                                             { "280",  { "70","80","90","100","125","150" } },
                                             { "300",  { "70","80","90","100","125","150","170" } },
                                             { "320",  { "70","80","90","100","125","150","170","190" } },
                                             { "360",  { "70","80","90","100","125","150","170","190" } },
                                             { "400",  { "70","80","90","100","125","150","170","190","250" } },
                                             { "450",  { "70","80","90","100","125","150","170","190","250","300" } },
                                             { "500",  { "70","80","90","100","125","150","170","190","250","300","350" } },
                                             { "550",  { "70","80","90","100","125","150","170","190","250","300","350","400" } },
                                             { "600",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450" } },
                                             { "650",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500" } },
                                             { "700",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550" } },
                                             { "750",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600" } },
                                             { "800",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600","650" } },
                                             { "850",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600","650","700" } },
                                             { "900",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600","650","700" } },
                                             { "950",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600","650","700","800" } },
                                             { "1000",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500","550","600","650","700","800" } } };

      map<TString,vector<TString>> M;
      if ( sampleArg.Contains("_low" ) ) // Includes MX of 1000, 240 - 400 GeV
        M.insert( MComb.begin(), next(MComb.begin(),7) );
      else if ( sampleArg.Contains("_med" ) ) // Includes MX of 450 - 750 GeV
        M.insert( next(MComb.begin(),7), next(MComb.begin(),14) );
      else if ( sampleArg.Contains("_high" ) ) // Includes MX of 800 - 950 GeV
        M.insert( next(MComb.begin(),14), MComb.end() );
      else
        M.insert( MComb.begin(),  MComb.end() );

      map<TString,vector<TString>>::iterator it = M.begin();

      while ( it != M.end() ) {
        TString MX = it->first;
        vector<TString> MYs = it->second;
        for ( auto MY : MYs ) {
          TString sampleName = "NMSSM_XToYHTo2G2B_MX-"+MX+"_MY-"+MY;
          TString v2018 = ( sampleName == "NMSSM_XToYHTo2G2B_MX-1000_MY-90"  ? "v1" : "v2" );
          TString v2017 = ( sampleName == "NMSSM_XToYHTo2G2B_MX-1000_MY-100" ? "v1" : "v2" );
          samples.push_back(sampleName);
          sample_procids.insert({sampleName, (processId_MX[MX]+processId_MY[MY]).Atoi()});
          sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-madgraph-pythia8"});
          sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-"+v2018 } },
                                            { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-"+v2017 } },
                                            { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                            { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });
        }
        ++it;
      }
    }

    if ( sampleArg.Contains("all") || sampleArg.Contains("inverted") ) {
      map<TString,vector<TString>> MComb = { { "650",  { "70","80","90","100","125" } } };

      map<TString,vector<TString>>::iterator it = MComb.begin();
      while ( it != MComb.end() ) {
        TString MX = it->first;
        vector<TString> MYs = it->second;
        for ( auto MY : MYs ) {
          TString sampleName = "NMSSM_XToYHTo2B2G_MX-"+MX+"_MY-"+MY;
          samples.push_back(sampleName);
          sample_procids.insert({sampleName, ("10"+processId_MX[MX]+processId_MY[MY]).Atoi()});
          sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-madgraph-pythia8"});
          sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                            { "2017",       { "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9-v2" } },
                                            { "2016APV",    { "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11-v2" } },
                                            { "2016nonAPV", { "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v2" } } } });
        }
        ++it;
      }
    }
  }

  if ( samples.size()==0 ) {
    std::cout << "No samples selected! Exiting...\n";
    return 0;
  }


  if (onlyCreateJSON) {
    // Create summary json
    ofstream file;
    file.open("summary.json");
    file << "{" << endl;
    file << "\t\"sample_id_map\": {" << endl;
    file << "\t\t\""<<samples[0]<<"\": "<<sample_procids[samples[0]];
    for ( int isample=1; isample<samples.size(); isample++ ) {
      file << "," <<endl;
      TString sample = samples[isample];
      if ( sample.Contains("NMSSM_XToYHTo") ) {
        sample.ReplaceAll("-","_");
      }
      file << "\t\t\""<<sample<<"\": "<<sample_procids[samples[isample]];
    }
    file << endl;
    file << "\t}" << endl;
    file << "}" << endl;
    file.close();
  }
  else {
    // Main loops
    TString baseDir = "/store/group/Hgg/XToYHToggbb/skimmedNanoAOD";
    TString version = "v0";

    for ( int iyear=0; iyear<years.size(); iyear++ ) {
      TString year = years[iyear];
      std::cout<<"\n";
      std::cout<<"Year: "<<year<<"\n";
      std::cout<<"------------\n";

      for ( int isample=0; isample<samples.size(); isample++ ) {
        TString sample = samples[isample];
        TString dataformat = "MINIAODSIM";
        bool isMC = 1;

        if ( sample == "Data" || sample == "DDQCDGJets" ) {
          sample_names[sample] = ( year=="2018" ? "EGamma" : "DoubleEG" );
          dataformat = "MINIAOD";
          isMC = 0;
        }
        TString sample_name = sample_names[sample];
        int sample_procid = sample_procids[sample];

        TChain *ch_temp = new TChain("Events");
        TChain *chaux_temp = new TChain("Runs");
        for ( unsigned int d=0; d<sample_prod[sample][year].size(); d++ ) {
          TString trees = "/ceph/cms"+baseDir+"/"+year+"/"+sample_name+"_"+sample_prod[sample][year][d]+"_"+dataformat+"_"+version+"/"+"tree_*.root"; // Local access
          //TString trees = "davs://redirector.t2.ucsd.edu:1095"+baseDir+"/"+year+"/"+sample_name+"_"+sample_prod[sample][year][d]+"_"+dataformat+"_"+version+"/"+"tree_*.root"; // Global access

          ch_temp->Add(trees);
          chaux_temp->Add(trees);
        }

        std::cout<<"Sample: "<<sample<<" --> Process ID: "<<sample_procid<<"\n\n";
        if ( ch_temp->GetEntries()==0 )
          std::cout << "##### Zero entries for the sample above => Need to check! #####\n";
        ScanChain_Hgg(ch_temp,getSumOfGenEventSumw(chaux_temp, isMC),year,sample,sample_procid,outdir,lowMassMode,prefireWeight,PUWeight,electronVetoSF,triggerSF,preselSF,phoMVAIDWP90SF,bTagSF,fnufUnc,materialUnc,PhoScaleUnc,PhoSmearUnc,JESUnc,JERUnc,HEMCheck);
      }
    }
  }

  return 0;
}
