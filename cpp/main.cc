#include "TROOT.h"
#include "ScanChain_Hgg.C"
#include <fstream>


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
  const char* outdir = ( argc > 1 ? argv[1]               : "temp_data" );
  TString yearArg    = ( argc > 2 ? argv[2]               : "all" );
  int run_data       = ( argc > 3 ? ((int)*argv[3] - 48 ) : 1 );  // '0' has the ASCII code of 48
  int run_MCbkg      = ( argc > 4 ? ((int)*argv[4] - 48 ) : 1 );  // '0' has the ASCII code of 48
  int run_signal     = ( argc > 5 ? ((int)*argv[5] - 48 ) : 1 );  // '0' has the ASCII code of 48
  TString sampleArg  = ( argc > 6 ? argv[6]               : "all" );
  int onlyCreateJSON = ( argc > 7 ? ((int)*argv[7] - 48 ) : 0 );  // '0' has the ASCII code of 48
  
  int PUWeight=1;
  int bTagSF=1;
  int JECUnc=0; // No central value, set to +/-2 to get
  

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
                                                          "Run2018D-UL2018_MiniAODv2-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
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
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }

    if ( sampleArg=="TT" || sampleArg=="all" ) {
      TString sampleName = "TTGG";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 2});
      sample_names.insert({sampleName, "TTGG_0Jets_TuneCP5_13TeV-amcatnlo-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "TTGJets";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 3});
      sample_names.insert({sampleName, "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "TTJets";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 4});
      sample_names.insert({sampleName, "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }

    if ( sampleArg=="H" || sampleArg=="all" ) {
      TString sampleName = "VH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 5});
      sample_names.insert({sampleName, "VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });  

      sampleName = "VBFH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 6});
      sample_names.insert({sampleName,"VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "ttH_M125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 7});
      sample_names.insert({sampleName, "ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "ggHToDiPhoM125";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 8});
      sample_names.insert({sampleName, "GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }

    if ( sampleArg=="GJets" || sampleArg=="all" ) {
      vector<TString> M = { "40", "100", "200", "400", "600", "Inf" };
      map<TString,int> processId_M = { {"40",0}, {"100",1}, {"200",2}, {"400",3}, {"600",4} };
      for ( unsigned int imass=0; imass<M.size()-1; imass++ ) {
        TString sampleName = "GJets_HT-"+M[imass]+"To"+M[imass+1];
        samples.push_back(sampleName);
        sample_procids.insert({sampleName, 9+processId_M[M[imass]]});
        sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-madgraphMLM-pythia8"});
        sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                          { "2017",       { "" } },
                                          { "2016APV",    { "" } },
                                          { "2016nonAPV", { "" } } } });
      }
    }

    if ( sampleArg=="DY" || sampleArg=="all" ) {
      TString sampleName = "DY";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 14});
      sample_names.insert({sampleName, "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }

    if ( sampleArg=="VG" || sampleArg=="all" ) {
      vector<TString> V = { "W", "Z" };
      map<TString,int> processId_V = { {"W",0}, {"Z",1} };
      for ( unsigned int iV=0; iV<V.size(); iV++ ) {
        TString sampleName = V[iV]+"G";
        samples.push_back(sampleName);
        sample_procids.insert({sampleName, 15+processId_V[V[iV]]});
        sample_names.insert({sampleName, sampleName+"ToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8"});
        sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                          { "2017",       { "" } },
                                          { "2016APV",    { "" } },
                                          { "2016nonAPV", { "" } } } });
      }
    }

    if ( sampleArg=="HHbbgg" || sampleArg=="all" ) {
      TString sampleName = "HHbbgg";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 17});
      sample_names.insert({sampleName, "GluGluToHHTo2B2G_node_cHHH1_TuneCP5_13TeV-powheg-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "privateUL18" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }

    if ( sampleArg=="DiPhotonLow" || sampleArg=="all" ) {
      TString sampleName = "DiPhotonLow";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 18});
      sample_names.insert({sampleName, "DiPhotonJetsBox_MGG-40to80_13TeV-sherpa"});
      sample_prod.insert({sampleName, { { "2018",       { "privateUL18Leonardo" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
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
                                          { "2017",       { "" } },
                                          { "2016APV",    { "" } },
                                          { "2016nonAPV", { "" } } } });
      }
    }

    if ( sampleArg=="TVX" || sampleArg=="all" ) {
      TString sampleName = "TTW";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 22});
      sample_names.insert({sampleName, "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "TTZ";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 23});
      sample_names.insert({sampleName, "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });

      sampleName = "tZq";
      samples.push_back(sampleName);
      sample_procids.insert({sampleName, 24});
      sample_names.insert({sampleName, "tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8"});
      sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2" } },
                                        { "2017",       { "" } },
                                        { "2016APV",    { "" } },
                                        { "2016nonAPV", { "" } } } });
    }
  }

  // Sample list: Signal
  if (run_signal) {
    if ( sampleArg.Contains("all" ) ) {
      map<TString,TString> processId_MX = { {"240","1"},  {"280","2"},  {"300","3"},  {"320","4"},  {"360","5"},
                                            {"400","6"},  {"450","7"},  {"500","8"},  {"550","9"},  {"600","10"},
                                            {"650","11"}, {"700","12"}, {"750","13"}, {"800","14"}, {"850","15"},
                                            {"900","16"}, {"950","17"}, {"1000","18"} };
      map<TString,TString> processId_MY = { {"70","11"},  {"80","12"},  {"90","13"},  {"100","14"}, {"125","15"},
                                            {"150","16"}, {"170","17"}, {"190","18"}, {"250","19"}, {"300","20"},
                                            {"350","21"}, {"400","22"}, {"450","23"}, {"500","24"}, {"550","25"},
                                            {"600","26"}, {"650","27"}, {"700","28"}, {"800","29"} };

      map<TString,vector<TString>> MComb = { { "240",  { "70","80","90","100" } },
                                             { "280",  { "70",/*"80",*/"90","100","125","150" } },
                                             { "300",  { "70",/*"80","90","100",*/"125","150","170" } },
                                             { "320",  { "70","80","90","100","125","150","170","190" } },
                                             { "360",  { "70",/*"80",*/"90","100","125","150","170","190" } },
                                             { "400",  { /*"70","80","90","100","125",*/"150","170","190","250" } },
                                             { "450",  { "70","80","90",/*"100",*/"125","150","170",/*"190",*/"250","300" } },
                                             { "500",  { "70","80","90",/*"100",*/"125",/*"150","170",*/"190","250","300","350" } },
                                             { "550",  { "70","80","90","100","125","150","170","190","250","300",/*"350",*/"400" } },
                                             { "600",  { /*"70",*/"80","90",/*"100",*/"125","150","170","190","250",/*"300","350",*/"400",
                                                         "450" } },
                                             { "650",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500" } },
                                             { "700",  { "70","80","90","100","125","150","170","190","250","300","350","400",
                                                         "450","500"/*,"550"*/ } },
                                             { "750",  { "70","80","90","100","125",/*"150",*/"170","190","250","300","350","400",
                                                         "450",/*"500",*/"550","600" } },
                                             { "800",  { "70",/*"80","90","100",*/"125","150",/*"170",*/"190","250","300","350","400",
                                                         "450","500",/*"550",*/"600","650" } },
                                             { "850",  { /*"70",*/"80","90",/*"100",*/"125","150","170","190","250","300","350","400",
                                                         "450","500",/*"550",*/"600"/*,"650","700"*/ } },
                                             { "900",  { /*"70",*/"80","90","100","125","150","170","190","250","300",/*"350",*/"400",
                                                         "450",/*"500",*/"550","600","650"/*,"700"*/ } },
                                             { "950",  { /*"70",*/"80","90","100","125","150",/*"170","190","250",*/"300","350","400",
                                                         "450","500","550","600",/*"650",*/"700","800" } },
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
          TString v = ( sampleName == "NMSSM_XToYHTo2G2B_MX-1000_MY-90" ? "v1" : "v2" );
          samples.push_back(sampleName);
          sample_procids.insert({sampleName, (processId_MX[MX]+processId_MY[MY]).Atoi()});
          sample_names.insert({sampleName, sampleName+"_TuneCP5_13TeV-madgraph-pythia8"});
          sample_prod.insert({sampleName, { { "2018",       { "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-"+v } },
                                            { "2017",       { "" } },
                                            { "2016APV",    { "" } },
                                            { "2016nonAPV", { "" } } } });
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
      if ( sample.Contains("NMSSM_XToYHTo2G2B") ) {
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
    TString baseDir = "/ceph/cms/store/group/Hgg/XToYHToggbb/skimmedNanoAOD";
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

        if (sample == "Data") {
          sample_names[sample] = ( year=="2018" ? "EGamma" : "DoubleEG" );
          dataformat = "MINIAOD";
          isMC = 0;
        }
        TString sample_name = sample_names[sample];
        int sample_procid = sample_procids[sample];

        TChain *ch_temp = new TChain("Events");
        TChain *chaux_temp = new TChain("Runs");
        for ( unsigned int d=0; d<sample_prod[sample][year].size(); d++ ) {
          TString trees = baseDir+"/"+year+"/"+sample_name+"_"+sample_prod[sample][year][d]+"_"+dataformat+"_"+version+"/"+"tree_*.root";

          ch_temp->Add(trees);
          chaux_temp->Add(trees);
        }

        std::cout<<"Sample: "<<sample<<" --> Process ID: "<<sample_procid<<"\n\n";
        ScanChain_Hgg(ch_temp,getSumOfGenEventSumw(chaux_temp, isMC),year,sample,sample_procid,outdir,PUWeight,bTagSF,JECUnc);
      }
    }
  }

  return 0;
}
