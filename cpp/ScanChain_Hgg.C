#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"
#include "TLorentzVector.h"
#include "math.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/Config.h"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/XYMETCorrection_withUL17andUL18andUL16.h"
#include "../NanoCORE/Tools/goodrun.h"
#include "../NanoCORE/Tools/dorky.h"
#include "../NanoCORE/Tools/bTagEff.h"
#include "../NanoCORE/Tools/electronVetoSF.h"
#include "../NanoCORE/Tools/lowMassHggTriggerSF.h"
#include "../NanoCORE/Tools/highMassHggTriggerSF.h"
#include "../NanoCORE/Tools/lowMassHggPreselSF.h"
#include "../NanoCORE/Tools/phoMVAIDWP90SF.h"
#include "../NanoCORE/DiPhotonSelections.h"
#include "../NanoCORE/LeptonSelections.h"
#include "../NanoCORE/DiJetSelections.h"
#include "../NanoCORE/GenPart.h"

#include <iostream>
#include <iomanip>
#include <sys/stat.h>
#include <fstream>

#define H1(name,nbins,low,high,xtitle) TH1D *h_##name = new TH1D(#name,"",nbins,low,high); h_##name->GetXaxis()->SetTitle(xtitle); h_##name->GetYaxis()->SetTitle("Events");

#define Zmass 91.1876

// General flags
bool removeSpikes = true;
bool removeDataDuplicates = false;
bool usePuppiMET = true;

using namespace std;
using namespace tas;
using namespace duplicate_removal;
int count_test=0;

ofstream txtout("evtnb.txt", ofstream::app);

void imputePhotonID(const float minPhotonID_cut, const float maxPhotonID, TF1* fakePhotonID_shape, float &minPhotonID, float &weight) {
    minPhotonID = fakePhotonID_shape->GetRandom(minPhotonID_cut, maxPhotonID);
    weight *= fakePhotonID_shape->Integral(minPhotonID_cut, maxPhotonID) / fakePhotonID_shape->Integral(-0.9, minPhotonID_cut);
    return;
}

TF1* get_fakePhotonIDShape(TString year, bool inclusive, bool isEndcap=true) {
  float lowerBound = -0.9;
  float upperBound = 1.0;

  TF1* fakePhotonMVAIDShape = new TF1("fakePhotonMVAIDShape", "pol9", lowerBound, upperBound);
  if ( year == "" ) {
    if ( inclusive ) {
      fakePhotonMVAIDShape->SetParameter(0, 1316.87);
      fakePhotonMVAIDShape->SetParameter(1,-2282.82);
      fakePhotonMVAIDShape->SetParameter(2,  1426.4);
      fakePhotonMVAIDShape->SetParameter(3, 2367.66);
      fakePhotonMVAIDShape->SetParameter(4, 23557.3);
      fakePhotonMVAIDShape->SetParameter(5,-38991.9);
      fakePhotonMVAIDShape->SetParameter(6,-58465.1);
      fakePhotonMVAIDShape->SetParameter(7, 85949.8);
      fakePhotonMVAIDShape->SetParameter(8, 59732.7);
      fakePhotonMVAIDShape->SetParameter(9,-74589.9);
    }
    else {
      if ( isEndcap ) {
        fakePhotonMVAIDShape->SetParameter(0, 232.801);
        fakePhotonMVAIDShape->SetParameter(1,-209.519);
        fakePhotonMVAIDShape->SetParameter(2, 121.689);
        fakePhotonMVAIDShape->SetParameter(3, 1152.68);
        fakePhotonMVAIDShape->SetParameter(4, 2567.51);
        fakePhotonMVAIDShape->SetParameter(5,-8273.67);
        fakePhotonMVAIDShape->SetParameter(6,-6191.54);
        fakePhotonMVAIDShape->SetParameter(7, 15893.3);
        fakePhotonMVAIDShape->SetParameter(8, 6767.25);
        fakePhotonMVAIDShape->SetParameter(9,-12053.2);
      }
      else {
        fakePhotonMVAIDShape->SetParameter(0, 1081.74);
        fakePhotonMVAIDShape->SetParameter(1,-2061.18);
        fakePhotonMVAIDShape->SetParameter(2, 1352.68);
        fakePhotonMVAIDShape->SetParameter(3, 1004.92);
        fakePhotonMVAIDShape->SetParameter(4, 20678.1);
        fakePhotonMVAIDShape->SetParameter(5,-29726.8);
        fakePhotonMVAIDShape->SetParameter(6,-51609.3);
        fakePhotonMVAIDShape->SetParameter(7, 68346.1);
        fakePhotonMVAIDShape->SetParameter(8, 52524.9);
        fakePhotonMVAIDShape->SetParameter(9,-61567.5);
      }
    }
  }
  else if ( year == "2018" ) {
    if ( inclusive ) {
      fakePhotonMVAIDShape->SetParameter(0, 417.324);
      fakePhotonMVAIDShape->SetParameter(1,-558.751);
      fakePhotonMVAIDShape->SetParameter(2, 228.889);
      fakePhotonMVAIDShape->SetParameter(3,  2616.9);
      fakePhotonMVAIDShape->SetParameter(4, 4602.03);
      fakePhotonMVAIDShape->SetParameter(5,  -16952);
      fakePhotonMVAIDShape->SetParameter(6,-10282.8);
      fakePhotonMVAIDShape->SetParameter(7, 31664.8);
      fakePhotonMVAIDShape->SetParameter(8, 10209.6);
      fakePhotonMVAIDShape->SetParameter(9,-21936.5);
    }
    else {
      if ( isEndcap ) {
        fakePhotonMVAIDShape->SetParameter(0, 62.5195);
        fakePhotonMVAIDShape->SetParameter(1,-53.1982);
        fakePhotonMVAIDShape->SetParameter(2,  85.929);
        fakePhotonMVAIDShape->SetParameter(3,  609.37);
        fakePhotonMVAIDShape->SetParameter(4, 138.294);
        fakePhotonMVAIDShape->SetParameter(5,-3274.68);
        fakePhotonMVAIDShape->SetParameter(6,-300.989);
        fakePhotonMVAIDShape->SetParameter(7, 5694.07);
        fakePhotonMVAIDShape->SetParameter(8, 514.687);
        fakePhotonMVAIDShape->SetParameter(9,-3471.98);
      }
      else {
        fakePhotonMVAIDShape->SetParameter(0,   353.3);
        fakePhotonMVAIDShape->SetParameter(1,-491.686);
        fakePhotonMVAIDShape->SetParameter(2, 169.783);
        fakePhotonMVAIDShape->SetParameter(3, 1752.05);
        fakePhotonMVAIDShape->SetParameter(4, 4265.64);
        fakePhotonMVAIDShape->SetParameter(5,-12428.9);
        fakePhotonMVAIDShape->SetParameter(6,-9515.45);
        fakePhotonMVAIDShape->SetParameter(7, 23796.1);
        fakePhotonMVAIDShape->SetParameter(8, 9360.25);
        fakePhotonMVAIDShape->SetParameter(9,-17236.2);
      }
    }
  }
  else if ( year == "2017" ) {
    if ( inclusive ) {
      fakePhotonMVAIDShape->SetParameter(0, 544.888);
      fakePhotonMVAIDShape->SetParameter(1,-1055.37);
      fakePhotonMVAIDShape->SetParameter(2, 999.277);
      fakePhotonMVAIDShape->SetParameter(3,-990.145);
      fakePhotonMVAIDShape->SetParameter(4, 10466.5);
      fakePhotonMVAIDShape->SetParameter(5,-10142.1);
      fakePhotonMVAIDShape->SetParameter(6,-27407.7);
      fakePhotonMVAIDShape->SetParameter(7, 28106.2);
      fakePhotonMVAIDShape->SetParameter(8, 29077.8);
      fakePhotonMVAIDShape->SetParameter(9,-29592.1);
    }
    else {
      if ( isEndcap ) {
        fakePhotonMVAIDShape->SetParameter(0,  85.5531);
        fakePhotonMVAIDShape->SetParameter(1,  -82.159);
        fakePhotonMVAIDShape->SetParameter(2,-0.582958);
        fakePhotonMVAIDShape->SetParameter(3,  269.774);
        fakePhotonMVAIDShape->SetParameter(4,  1266.19);
        fakePhotonMVAIDShape->SetParameter(5, -2380.84);
        fakePhotonMVAIDShape->SetParameter(6, -2755.71);
        fakePhotonMVAIDShape->SetParameter(7,  4592.77);
        fakePhotonMVAIDShape->SetParameter(8,  2861.17);
        fakePhotonMVAIDShape->SetParameter(9, -3853.95);
      }
      else {
        fakePhotonMVAIDShape->SetParameter(0, 457.937);
        fakePhotonMVAIDShape->SetParameter(1,-972.254);
        fakePhotonMVAIDShape->SetParameter(2, 1010.73);
        fakePhotonMVAIDShape->SetParameter(3,-1298.21);
        fakePhotonMVAIDShape->SetParameter(4, 9155.69);
        fakePhotonMVAIDShape->SetParameter(5,-7514.89);
        fakePhotonMVAIDShape->SetParameter(6,-24576.8);
        fakePhotonMVAIDShape->SetParameter(7, 23014.8);
        fakePhotonMVAIDShape->SetParameter(8, 26169.4);
        fakePhotonMVAIDShape->SetParameter(9,-25429.9);
      }
    }
  }
  else if ( year == "2016APV" ) {
    if ( inclusive ) {
      fakePhotonMVAIDShape->SetParameter(0, 175.978);
      fakePhotonMVAIDShape->SetParameter(1,-330.624);
      fakePhotonMVAIDShape->SetParameter(2, 230.033);
      fakePhotonMVAIDShape->SetParameter(3, 239.253);
      fakePhotonMVAIDShape->SetParameter(4, 3305.56);
      fakePhotonMVAIDShape->SetParameter(5,-5200.57);
      fakePhotonMVAIDShape->SetParameter(6,-8306.41);
      fakePhotonMVAIDShape->SetParameter(7, 11571.7);
      fakePhotonMVAIDShape->SetParameter(8, 8885.88);
      fakePhotonMVAIDShape->SetParameter(9,-10579.4);
    }
    else {
      if ( isEndcap ) {
        fakePhotonMVAIDShape->SetParameter(0,  39.861);
        fakePhotonMVAIDShape->SetParameter(1,-16.9761);
        fakePhotonMVAIDShape->SetParameter(2, 125.197);
        fakePhotonMVAIDShape->SetParameter(3,-109.174);
        fakePhotonMVAIDShape->SetParameter(4, 12.2005);
        fakePhotonMVAIDShape->SetParameter(5,-380.418);
        fakePhotonMVAIDShape->SetParameter(6,-518.982);
        fakePhotonMVAIDShape->SetParameter(7, 1344.25);
        fakePhotonMVAIDShape->SetParameter(8, 1087.55);
        fakePhotonMVAIDShape->SetParameter(9,-1574.44);
      }
      else {
        fakePhotonMVAIDShape->SetParameter(0, 135.327);
        fakePhotonMVAIDShape->SetParameter(1,-310.684);
        fakePhotonMVAIDShape->SetParameter(2, 111.795);
        fakePhotonMVAIDShape->SetParameter(3, 300.753);
        fakePhotonMVAIDShape->SetParameter(4, 3222.67);
        fakePhotonMVAIDShape->SetParameter(5,-4593.07);
        fakePhotonMVAIDShape->SetParameter(6,   -7618);
        fakePhotonMVAIDShape->SetParameter(7, 9830.08);
        fakePhotonMVAIDShape->SetParameter(8, 7682.17);
        fakePhotonMVAIDShape->SetParameter(9,-8776.85);
      }
    }
  }
  else if ( year == "2016nonAPV" ) {
    if ( inclusive ) {
      fakePhotonMVAIDShape->SetParameter(0, 171.478);
      fakePhotonMVAIDShape->SetParameter(1,-322.973);
      fakePhotonMVAIDShape->SetParameter(2, 156.348);
      fakePhotonMVAIDShape->SetParameter(3, 147.162);
      fakePhotonMVAIDShape->SetParameter(4,  3848.3);
      fakePhotonMVAIDShape->SetParameter(5,-4779.66);
      fakePhotonMVAIDShape->SetParameter(6,-9530.07);
      fakePhotonMVAIDShape->SetParameter(7, 10969.1);
      fakePhotonMVAIDShape->SetParameter(8, 9591.86);
      fakePhotonMVAIDShape->SetParameter(9,-10248.4);
    }
    else {
      if ( isEndcap ) {
        fakePhotonMVAIDShape->SetParameter(0, 42.0966);
        fakePhotonMVAIDShape->SetParameter(1,-55.5476);
        fakePhotonMVAIDShape->SetParameter(2,-71.3172);
        fakePhotonMVAIDShape->SetParameter(3, 355.575);
        fakePhotonMVAIDShape->SetParameter(4, 1011.38);
        fakePhotonMVAIDShape->SetParameter(5,-2024.01);
        fakePhotonMVAIDShape->SetParameter(6,-2281.21);
        fakePhotonMVAIDShape->SetParameter(7, 3764.46);
        fakePhotonMVAIDShape->SetParameter(8, 2064.11);
        fakePhotonMVAIDShape->SetParameter(9, -2814.4);
      }
      else {
        fakePhotonMVAIDShape->SetParameter(0, 128.576);
        fakePhotonMVAIDShape->SetParameter(1,-267.889);
        fakePhotonMVAIDShape->SetParameter(2, 217.171);
        fakePhotonMVAIDShape->SetParameter(3,-177.811);
        fakePhotonMVAIDShape->SetParameter(4, 2922.14);
        fakePhotonMVAIDShape->SetParameter(5,-2913.87);
        fakePhotonMVAIDShape->SetParameter(6,-7441.72);
        fakePhotonMVAIDShape->SetParameter(7,  7476.9);
        fakePhotonMVAIDShape->SetParameter(8, 7656.73);
        fakePhotonMVAIDShape->SetParameter(9,-7587.61);
      }
    }
  }
  else {
    cout<<"Non-valid process: No fakePhotonMVAIDShape!"<<endl;
    return nullptr;
  }

  return fakePhotonMVAIDShape;
}

int ScanChain_Hgg(TChain *ch, double genEventSumw, TString year, TString process, int process_id, const char* outdir="temp_data", int lowMassMode=1, int prefireWeight=1, int PUWeight=1, int electronVetoSF=1, int triggerSF=1, int preselSF=1, int phoMVAIDWP90SF=1, int bTagSF=1, int fnufUnc=0, int materialUnc=0, int PhoScaleUnc=0, int PhoSmearUnc=0, int JESUnc=0, int JERUnc=0, int HEMCheck=0, int enrichDY=0) {
// Event weights / scale factors:
//  0: Do not apply
//  1: Apply central value
// +X: Apply positive variation X, X>1 
// -X: Apply negative variation X, X>1 

  float factor = 1.0;
  float lumi = 1.0;
  float xsec = 1.0;
  bool isMC = true;
  if ( process == "Data" || process == "DDQCDGJets" ) {
    isMC = false;
  }
   
  // Processes and cross-sections (in fb):
  // set this in a different file
  else if ( process == "ttbar" )                      { xsec = 87310.0;                        }
  else if ( process == "DY" )                         { xsec = 5765400.0;                      }
  else if ( process == "WW" )                         { xsec = 118700.0;                       }
  else if ( process == "WZ" )                         { xsec = 47130.0;                        }
  else if ( process == "ZZ" )                         { xsec = 16523.0;                        }
  else if ( process == "tW" )                         { xsec = 19550;                          }
  else if ( process == "tbarW" )                      { xsec = 19550;                          }
  else if ( process == "tZq" )                        { xsec = 75.8;                           }
  else if ( process == "TTW" )                        { xsec = 204.3;                          }
  else if ( process == "TTZ" )                        { xsec = 252.9;                          }
  else if ( process == "TTHToNonbb" )                 { xsec = 507.5*(1-0.575);                }
  else if ( process == "TTHTobb" )                    { xsec = 507.5*0.575;                    }
  else if ( process == "TTGG" )                       { xsec = 0.01687 * 1000;                 }
  else if ( process == "TTGJets" )                    { xsec = 4.078 * 1000;                   }
  else if ( process == "TTJets" )                     { xsec = 831.76 * 1000;                  }
  else if ( process == "VBFH_M125" )                  { xsec = 0.00858514 *1000;               }
  else if ( process == "VH_M125" )                    { xsec = 0.00512 *1000;                  }
  else if ( process == "ggHToDiPhoM125" )             { xsec = 0.1118429*1000 ;                }
  else if ( process == "ttH_M125" )                   { xsec = 0.5071 * 1000 * 0.00227;        }
  else if ( process == "GJets_HT-40To100" )           { xsec = 23100*1000;                     }
  else if ( process == "GJets_HT-100To200" )          { xsec = 8631.0*1000;                    }
  else if ( process == "GJets_HT-200To400" )          { xsec = 2280.0*1000;                    }
  else if ( process == "GJets_HT-400To600" )          { xsec = 273*1000;                       }
  else if ( process == "GJets_HT-600ToInf" )          { xsec = 1*1000;                         }
  else if ( process == "DiPhoton" )                   { xsec = 84.4*1000;                      }
  else if ( process == "DiPhotonLow" )                { xsec = 303.2*1000;                     }
  else if ( process == "HHbbgg" )                     { xsec = 0.03105*1000*0.00262230;        }
  else if ( process == "WG" )                         { xsec = 191.4*1000;                     }
  else if ( process == "ZG" )                         { xsec = 55.6*1000;                      }
  else if ( process == "WW" )                         { xsec = 75.8*1000;                      }
  else if ( process == "WZ" )                         { xsec = 27.6*1000;                      }
  else if ( process == "ZZ" )                         { xsec = 12.14*1000;                     }
  else if ( process == "QCD_Pt-30to40_MGG-80toInf" )  { xsec = 24810.0*1000;                   }
  else if ( process == "QCD_Pt-30toInf_MGG-40to80" )  { xsec = 241400.0*1000;                  }
  else if ( process == "QCD_Pt-40toInf_MGG-80toInf" ) { xsec = 118100.0*1000;                  }
  else if ( process == "DDQCDGJets" )                 { xsec = 1;                              }
  else if ( process == "Data" )                       { xsec = 1;                              }
  else if ( process.Contains("NMSSM_XToYHTo2G2B") )   { xsec = 1; process.ReplaceAll("-","_"); }
  else if ( process.Contains("NMSSM_XToYHTo2B2G") )   { xsec = 1; process.ReplaceAll("-","_"); }
  else {
    cout<<"Non-valid process: Exiting!"<<endl;
    return 1;
  }

  // Configuration setup: NanoCORE/Config.{h,cc}
  gconf.nanoAOD_ver = 9;
  gconf.GetConfigs(year.Atoi());
  lumi = gconf.lumi;
  if (year == "2018") lumi = lowMassMode ? 54.5 : 59.8;
  if (year == "2017") lumi = 41.5; // Same for low mass and high mass triggers
  if (year == "2016APV") lumi = 19.5; // Same for low mass and high mass triggers
  if (year == "2016nonAPV") lumi = 16.8; // Same for low mass and high mass triggers

  // Golden JSON files
  if ( !isMC ) {
    if ( year == "2016APV" )
      set_goodrun_file_json("../utils/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    if ( year == "2016nonAPV" )
      set_goodrun_file_json("../utils/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    if ( year == "2017" )
      set_goodrun_file_json("../utils/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt");
    if ( year == "2018" )
      set_goodrun_file_json("../utils/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt");
  }

  if ( isMC )
    factor = xsec*lumi/genEventSumw;


  //all Modify the name of the output file to include arguments of ScanChain function (i.e. process, year, etc.)
  int mdir = mkdir(outdir,0755);
  TString oDir(outdir);
  TFile* fout = new TFile(oDir+"/output_"+process+"_"+year+".root", "RECREATE");
  TTree* tout = new TTree("tout","Tree with photon variables");

  // define histograms, to be put in a different file TODO
  H1(LeadPhoton_sieie, 20, 0, 0.05, "");
  H1(LeadPhoton_pfPhoIso03, 20, 0, 10, "");
  H1(LeadPhoton_chargedHadronIso, 20, 0, 10, "");
  H1(LeadPhoton_trkSumPtHollowConeDR03, 20, 0, 10, "");
  H1(SubleadPhoton_sieie, 20, 0, 0.05, "");
  H1(SubleadPhoton_pfPhoIso03, 20, 0, 10, "");
  H1(SubleadPhoton_chargedHadronIso, 20, 0, 10, "");
  H1(SubleadPhoton_trkSumPtHollowConeDR03, 20, 0, 10, "");


  // Variables for output branches
  float xcand_pt, xcand_eta, xcand_phi, xcand_mass;

  float LeadPhoton_pt, LeadPhoton_eta, LeadPhoton_phi, LeadPhoton_mvaID;
  float SubleadPhoton_pt, SubleadPhoton_eta, SubleadPhoton_phi, SubleadPhoton_mvaID;
  float Diphoton_pt, Diphoton_eta, Diphoton_phi, Diphoton_mass, Diphoton_pt_mgg, Diphoton_dR;
  float LeadPhoton_sieie, LeadPhoton_pfPhoIso03, LeadPhoton_chargedHadronIso, LeadPhoton_r9, LeadPhoton_trkSumPtHollowConeDR03;
  float SubleadPhoton_sieie, SubleadPhoton_pfPhoIso03, SubleadPhoton_chargedHadronIso, SubleadPhoton_r9, SubleadPhoton_trkSumPtHollowConeDR03;
  bool LeadPhoton_pixelSeed, SubleadPhoton_pixelSeed;

  int n_jets;
  float dijet_lead_pt, dijet_lead_eta, dijet_lead_phi, dijet_lead_mass, dijet_lead_btagDeepFlavB;
  float dijet_sublead_pt, dijet_sublead_eta, dijet_sublead_phi, dijet_sublead_mass, dijet_sublead_btagDeepFlavB;
  float dijet_pt, dijet_eta, dijet_phi, dijet_mass, dijet_dR;
  float pfmet_pt, puppimet_pt;
  int year_out, eventNum;
  float weight_central, weight_central_initial, weight_central_no_lumi, weight_beforeBTagSF, weight_afterBTagSF;

  unsigned int LeadPhoton_genPartFlav, SubleadPhoton_genPartFlav;
  int n_gen_matched_jets, n_gen_matched_in_dijet;
  bool dijet_lead_gen_match, dijet_sublead_gen_match;
  float GenHiggs_pt, GenHiggs_eta, GenHiggs_phi, GenHiggs_mass, GenHiggs_dR;
  float GenY_pt, GenY_eta, GenY_phi, GenY_mass, GenY_dR;
  float GenX_pt, GenX_eta, GenX_phi, GenX_mass, GenX_dR;
  float GenBFromHiggs_1_pt, GenBFromHiggs_1_eta, GenBFromHiggs_1_phi, GenBFromHiggs_1_mass;
  float GenBFromHiggs_2_pt, GenBFromHiggs_2_eta, GenBFromHiggs_2_phi, GenBFromHiggs_2_mass;


  // Branch booking
  tout->Branch("xcand_pt", &xcand_pt, "xcand_pt/F");
  tout->Branch("xcand_eta", &xcand_eta, "xcand_eta/F");
  tout->Branch("xcand_phi", &xcand_phi, "xcand_phi/F");
  tout->Branch("xcand_mass", &xcand_mass, "xcand_mass/F");

  tout->Branch("LeadPhoton_pt",&LeadPhoton_pt,"LeadPhoton_pt/F");
  tout->Branch("LeadPhoton_eta",&LeadPhoton_eta,"LeadPhoton_eta/F");
  tout->Branch("LeadPhoton_phi",&LeadPhoton_phi,"LeadPhoton_phi/F");
  tout->Branch("LeadPhoton_pixelSeed",&LeadPhoton_pixelSeed,"LeadPhoton_pixelSeed/B");
  tout->Branch("LeadPhoton_r9",&LeadPhoton_r9,"LeadPhoton_r9/F");
  tout->Branch("LeadPhoton_sieie",&LeadPhoton_sieie,"LeadPhoton_sieie/F");
  tout->Branch("LeadPhoton_pfPhoIso03",&LeadPhoton_pfPhoIso03,"LeadPhoton_pfPhoIso03/F");
  tout->Branch("LeadPhoton_trkSumPtHollowConeDR03",&LeadPhoton_trkSumPtHollowConeDR03,"LeadPhoton_trkSumPtHollowConeDR03/F");
  tout->Branch("LeadPhoton_chargedHadronIso",&LeadPhoton_chargedHadronIso,"LeadPhoton_chargedHadronIso/F");
  tout->Branch("LeadPhoton_mvaID",&LeadPhoton_mvaID,"LeadPhoton_mvaID/F");

  tout->Branch("SubleadPhoton_pt",&SubleadPhoton_pt,"SubleadPhoton_pt/F");
  tout->Branch("SubleadPhoton_eta",&SubleadPhoton_eta,"SubleadPhoton_eta/F");
  tout->Branch("SubleadPhoton_phi",&SubleadPhoton_phi,"SubleadPhoton_phi/F");
  tout->Branch("SubleadPhoton_pixelSeed",&SubleadPhoton_pixelSeed,"SubleadPhoton_pixelSeed/B");  
  tout->Branch("SubleadPhoton_r9",&SubleadPhoton_r9,"SubleadPhoton_r9/F");
  tout->Branch("SubleadPhoton_sieie",&SubleadPhoton_sieie,"SubleadPhoton_sieie/F");
  tout->Branch("SubleadPhoton_pfPhoIso03",&SubleadPhoton_pfPhoIso03,"SubleadPhoton_pfPhoIso03/F");
  tout->Branch("SubleadPhoton_trkSumPtHollowConeDR03",&SubleadPhoton_trkSumPtHollowConeDR03,"SubleadPhoton_trkSumPtHollowConeDR03/F");
  tout->Branch("SubleadPhoton_chargedHadronIso",&SubleadPhoton_chargedHadronIso,"SubleadPhoton_chargedHadronIso/F");
  tout->Branch("SubleadPhoton_mvaID",&SubleadPhoton_mvaID,"SubleadPhoton_mvaID/F");

  tout->Branch("Diphoton_pt",&Diphoton_pt,"Diphoton_pt/F");
  tout->Branch("Diphoton_eta",&Diphoton_eta,"Diphoton_eta/F");
  tout->Branch("Diphoton_phi",&Diphoton_phi,"Diphoton_phi/F");
  tout->Branch("Diphoton_mass",&Diphoton_mass,"Diphoton_mass/F");
  tout->Branch("Diphoton_pt_mgg",&Diphoton_pt_mgg,"Diphoton_pt_mgg/F");
  tout->Branch("Diphoton_dR",&Diphoton_dR,"Diphoton_dR/F");

  tout->Branch("n_jets",&n_jets,"n_jets/I");  
  tout->Branch("dijet_lead_pt",&dijet_lead_pt,"dijet_lead_pt/F");
  tout->Branch("dijet_lead_eta",&dijet_lead_eta,"dijet_lead_eta/F");
  tout->Branch("dijet_lead_phi",&dijet_lead_phi,"dijet_lead_phi/F");
  tout->Branch("dijet_lead_mass",&dijet_lead_mass,"dijet_lead_mass/F");
  tout->Branch("dijet_lead_btagDeepFlavB",&dijet_lead_btagDeepFlavB,"dijet_lead_btagDeepFlavB/F");

  tout->Branch("dijet_sublead_pt",&dijet_sublead_pt,"dijet_sublead_pt/F");
  tout->Branch("dijet_sublead_eta",&dijet_sublead_eta,"dijet_sublead_eta/F");
  tout->Branch("dijet_sublead_phi",&dijet_sublead_phi,"dijet_sublead_phi/F");
  tout->Branch("dijet_sublead_mass",&dijet_sublead_mass,"dijet_sublead_mass/F");
  tout->Branch("dijet_sublead_btagDeepFlavB",&dijet_sublead_btagDeepFlavB,"dijet_sublead_btagDeepFlavB/F");

  tout->Branch("dijet_pt",&dijet_pt,"dijet_pt/F");
  tout->Branch("dijet_eta",&dijet_eta,"dijet_eta/F");  
  tout->Branch("dijet_phi",&dijet_phi,"dijet_phi/F");  
  tout->Branch("dijet_mass",&dijet_mass,"dijet_mass/F");  
  tout->Branch("dijet_dR",&dijet_dR,"dijet_dR/F"); 

  tout->Branch("pfmet_pt",&pfmet_pt,"pfmet_pt/F"); 
  tout->Branch("puppimet_pt",&puppimet_pt,"puppimet_pt/F"); 

  tout->Branch("year",&year_out,"year/I");
  tout->Branch("weight_central",&weight_central,"weight_central/F");
  tout->Branch("weight_central_initial",&weight_central_initial,"weight_central_initial/F");
  tout->Branch("weight_central_no_lumi",&weight_central_no_lumi,"weight_central_no_lumi/F");
  tout->Branch("weight_beforeBTagSF",&weight_beforeBTagSF,"weight_beforeBTagSF/F");
  tout->Branch("weight_afterBTagSF",&weight_afterBTagSF,"weight_afterBTagSF/F");
  tout->Branch("event",&eventNum,"event/I");
  tout->Branch("process_id",&process_id,"process_id/I");

  if (year=="2016nonAPV" || year=="2016APV") year_out = 2016;
  else if (year=="2017") year_out = 2017;
  else if (year=="2018") year_out = 2018;
  else year_out = 0;
  
  tout->Branch("LeadPhoton_genPartFlav",&LeadPhoton_genPartFlav,"LeadPhoton_genPartFlav/I");
  tout->Branch("SubleadPhoton_genPartFlav",&SubleadPhoton_genPartFlav,"SubleadPhoton_genPartFlav/I");
  tout->Branch("n_gen_matched_jets",&n_gen_matched_jets,"n_gen_matched_jets/I");
  tout->Branch("n_gen_matched_in_dijet",&n_gen_matched_in_dijet,"n_gen_matched_in_dijet/I");
  tout->Branch("dijet_lead_gen_match",&dijet_lead_gen_match,"dijet_lead_gen_match/B");
  tout->Branch("dijet_sublead_gen_match",&dijet_sublead_gen_match,"dijet_sublead_gen_match/B");
  tout->Branch("GenHiggs_pt",&GenHiggs_pt,"GenHiggs_pt/F");
  tout->Branch("GenHiggs_eta",&GenHiggs_eta,"GenHiggs_eta/F");
  tout->Branch("GenHiggs_phi",&GenHiggs_phi,"GenHiggs_phi/F");
  tout->Branch("GenHiggs_mass",&GenHiggs_mass,"GenHiggs_mass/F");
  tout->Branch("GenHiggs_dR",&GenHiggs_dR,"GenHiggs_dR/F");
  tout->Branch("GenY_pt",&GenY_pt,"GenY_pt/F");
  tout->Branch("GenY_eta",&GenY_eta,"GenY_eta/F");
  tout->Branch("GenY_phi",&GenY_phi,"GenY_phi/F");
  tout->Branch("GenY_mass",&GenY_mass,"GenY_mass/F");
  tout->Branch("GenY_dR",&GenY_dR,"GenY_dR/F");
  tout->Branch("GenX_pt",&GenX_pt,"GenX_pt/F");
  tout->Branch("GenX_eta",&GenX_eta,"GenX_eta/F");
  tout->Branch("GenX_phi",&GenX_phi,"GenX_phi/F");
  tout->Branch("GenX_mass",&GenX_mass,"GenX_mass/F");
  tout->Branch("GenX_dR",&GenX_dR,"GenX_dR/F");
  tout->Branch("GenBFromHiggs_1_pt",&GenBFromHiggs_1_pt,"GenBFromHiggs_1_pt/F");
  tout->Branch("GenBFromHiggs_1_eta",&GenBFromHiggs_1_eta,"GenBFromHiggs_1_eta/F");
  tout->Branch("GenBFromHiggs_1_phi",&GenBFromHiggs_1_phi,"GenBFromHiggs_1_phi/F");
  tout->Branch("GenBFromHiggs_1_mass",&GenBFromHiggs_1_mass,"GenBFromHiggs_1_mass/F");
  tout->Branch("GenBFromHiggs_2_pt",&GenBFromHiggs_2_pt,"GenBFromHiggs_2_pt/F");
  tout->Branch("GenBFromHiggs_2_eta",&GenBFromHiggs_2_eta,"GenBFromHiggs_2_eta/F");
  tout->Branch("GenBFromHiggs_2_phi",&GenBFromHiggs_2_phi,"GenBFromHiggs_2_phi/F");
  tout->Branch("GenBFromHiggs_2_mass",&GenBFromHiggs_2_mass,"GenBFromHiggs_2_mass/F");


  // Define histo info maps
  map<TString, int> nbins { };
  map<TString, float> low { };
  map<TString, float> high { };
  map<TString, vector<float>> binsx { };
  map<TString, TString> title { };

  // Define histos
  H1(cutflow,20,0,20,"");
  H1(weight,1,0,1,"");
  H1(weight_full,1,0,1,"");
  H1(weight_beforeBTagSF,1,0,1,"");
  H1(weight_afterBTagSF,1,0,1,"");

  TF1 *fakePhotonID_shape_barrel = get_fakePhotonIDShape(year,/*inclusive*/false,/*isEndcap*/false);
  TF1 *fakePhotonID_shape_endcap = get_fakePhotonIDShape(year,/*inclusive*/false,/*isEndcap*/true);
  if ( electronVetoSF != 0) electronVetoSF::set_electronVetoSF();
  if ( triggerSF != 0) (lowMassMode ? lowMassHggTriggerSF::set_lowMassHggTriggerSF() : highMassHggTriggerSF::set_highMassHggTriggerSF() );
  if ( preselSF != 0) lowMassHggPreselSF::set_lowMassHggPreselSF();
  if ( phoMVAIDWP90SF != 0) phoMVAIDWP90SF::set_phoMVAIDWP90SF();

  int nEventsTotal = 0;
  int nDuplicates = 0;
  int nEventsChain = ch->GetEntries();
  TFile *currentFile = 0;
  TObjArray *listOfFiles = ch->GetListOfFiles();
  TIter fileIter(listOfFiles);
  tqdm bar;

  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TFile *file = TFile::Open( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    TString filename(currentFile->GetTitle());

    tree->SetCacheSize(128*1024*1024);
    tree->SetCacheLearnEntries(100);

    nt.Init(tree);

    for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
      nt.GetEntry(event);
      tree->LoadTree(event);

      nEventsTotal++;
      bar.progress(nEventsTotal, nEventsChain);

      //initialize variables in each event loop
      xcand_pt=-999, xcand_eta=-999, xcand_phi=-999, xcand_mass=-999;

      LeadPhoton_pt=-999, LeadPhoton_eta=-999, LeadPhoton_phi=-999, LeadPhoton_mvaID=-999;
      SubleadPhoton_pt=-999, SubleadPhoton_eta=-999, SubleadPhoton_phi=-999, SubleadPhoton_mvaID=-999;
      Diphoton_pt=-999, Diphoton_eta=-999, Diphoton_phi=-999, Diphoton_mass=-999, Diphoton_pt_mgg=-999, Diphoton_dR=-999;
      LeadPhoton_sieie=-999, LeadPhoton_pfPhoIso03=-999, LeadPhoton_chargedHadronIso=-999, LeadPhoton_r9=-999, LeadPhoton_trkSumPtHollowConeDR03=-999;
      SubleadPhoton_sieie=-999, SubleadPhoton_pfPhoIso03=-999, SubleadPhoton_chargedHadronIso=-999, SubleadPhoton_r9=-999, SubleadPhoton_trkSumPtHollowConeDR03=-999;
      LeadPhoton_pixelSeed=true, SubleadPhoton_pixelSeed=true;

      n_jets=-1;
      dijet_lead_pt=-999, dijet_lead_eta=-999, dijet_lead_phi=-999, dijet_lead_mass=-999, dijet_lead_btagDeepFlavB=-999;
      dijet_sublead_pt=-999, dijet_sublead_eta=-999, dijet_sublead_phi=-999, dijet_sublead_mass=-999, dijet_sublead_btagDeepFlavB=-999;
      dijet_pt=-999, dijet_eta=-999, dijet_phi=-999, dijet_mass=-999, dijet_dR=-999;
      pfmet_pt=-999, puppimet_pt=-999;
      eventNum=0;
      weight_central=1.0, weight_central_initial=1.0, weight_central_no_lumi=1.0, weight_beforeBTagSF=1.0, weight_afterBTagSF=1.0;

      LeadPhoton_genPartFlav=0; SubleadPhoton_genPartFlav=0;
      n_gen_matched_jets=0; n_gen_matched_in_dijet=0;
      dijet_lead_gen_match=false; dijet_sublead_gen_match=false;
      GenHiggs_pt=-999; GenHiggs_eta=-999; GenHiggs_phi=-999; GenHiggs_mass=-999; GenHiggs_dR=-999;
      GenY_pt=-999; GenY_eta=-999; GenY_phi=-999; GenY_mass=-999; GenY_dR=-999;
      GenX_pt=-999; GenX_eta=-999; GenX_phi=-999; GenX_mass=-999; GenX_dR=-999;
      GenBFromHiggs_1_pt=-999; GenBFromHiggs_1_eta=-999; GenBFromHiggs_1_phi=-999; GenBFromHiggs_1_mass=-999;
      GenBFromHiggs_2_pt=-999; GenBFromHiggs_2_eta=-999; GenBFromHiggs_2_phi=-999; GenBFromHiggs_2_mass=-999;

      float weight = 1.0;
      if ( isMC ) {
        weight = nt.genWeight();

        // Apply PreFiring weight
        if ( prefireWeight!=0 ) {
          if ( prefireWeight==1  ) weight *= nt.L1PreFiringWeight_Nom();
          if ( prefireWeight==2  ) weight *= nt.L1PreFiringWeight_Up();
          if ( prefireWeight==-2 ) weight *= nt.L1PreFiringWeight_Dn();
        }
        //if(removeSpikes && weight*factor>1e2) continue; //comment out for synchronizing

        // Apply PU weight
        if ( PUWeight!=0 ) {
          if ( PUWeight==1  ) weight *= nt.puWeight();
          if ( PUWeight==2  ) weight *= nt.puWeightUp();
          if ( PUWeight==-2 ) weight *= nt.puWeightDown();
        }
      }

      h_weight_full->Fill(0.5, weight*factor);

      unsigned int runnb = nt.run();
      unsigned int lumiblock = nt.luminosityBlock();
      unsigned long int evtnb = nt.event();
      int npv = nt.PV_npvs();


      // Apply Golden JSON
      if ( !isMC ) {
        if ( !(goodrun(runnb, lumiblock)) )
          continue;
        if ( removeDataDuplicates ) {
          DorkyEventIdentifier id(runnb, evtnb, lumiblock);
          if ( is_duplicate(id) ) {
            ++nDuplicates;
            continue;
          }
        }
      }

      // HLT selection
      if ( !isMC ) {
        if ( enrichDY ) {
          // HLT selection Ele path (inverted)
          if ( (year=="2016nonAPV" || year=="2016APV") &&
              !( (tree->GetBranch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ? nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() : 0)
                || (tree->GetBranch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ? nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() : 0) ) ) continue;
          if ( (year=="2017") &&
              !( (tree->GetBranch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ? nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() : 0)  )  ) continue;
          if ( (year=="2018") &&
              !( (tree->GetBranch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ? nt.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() : 0) ) ) continue;
        }
        else{
            if ( lowMassMode ) {
              if ( (year=="2016nonAPV" || year=="2016APV") &&
                  !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0)
                    || (tree->GetBranch("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55") ? nt.HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55() : 0) ) ) continue;
              if ( (year=="2017") &&
                  !( (tree->GetBranch("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55") ? nt.HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55() : 0)  )  ) continue;
              if ( (year=="2018") &&
                  !( (tree->GetBranch("HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto") ? nt.HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto() : 0) ) ) continue;
            }
            else {
              if ( (year=="2016nonAPV" || year=="2016APV") &&
                  !( (tree->GetBranch("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90") ? nt.HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() : 0) ) ) continue;
              if ( (year=="2017") &&
                  !( (tree->GetBranch("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90") ? nt.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() : 0) ) ) continue;
              if ( (year=="2018") &&
                  !( (tree->GetBranch("HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90") ? nt.HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() : 0) ) ) continue;
            }
        }

      }

      // MET filters
      if ( !( nt.Flag_goodVertices()>=1 &&
              nt.Flag_globalSuperTightHalo2016Filter()>=1 &&
              nt.Flag_HBHENoiseFilter()>=1 &&
              nt.Flag_HBHENoiseIsoFilter()>=1 &&
              nt.Flag_EcalDeadCellTriggerPrimitiveFilter()>=1 &&
              nt.Flag_BadPFMuonFilter()>=1 &&
              nt.Flag_BadPFMuonDzFilter()>=1 &&
              nt.Flag_hfNoisyHitsFilter()>=1 &&
              nt.Flag_eeBadScFilter()>=1 &&
              ( year.Contains("2016") ? 1 : nt.Flag_ecalBadCalibFilter()>=1 ) ) )
      continue;

      // Object selection
      Photons photons;
      if (isMC) {
        photons = getPhotons(year, lowMassMode, fnufUnc, materialUnc, PhoScaleUnc, PhoSmearUnc, enrichDY);
      }
      else {
        photons = getPhotons(year, lowMassMode, 0, 0, 0, 0, enrichDY);
      }
      DiPhotons diphotons = DiPhotonPreselection(photons, lowMassMode);

      if (diphotons.size() == 0 ) continue; 

      DiPhoton selectedDiPhoton = diphotons[0];
      Photons selectedPhotons={selectedDiPhoton.leadPho, selectedDiPhoton.subleadPho};

      // Get DDQCDGJets from data
      float minPhotonID_cut = -1.0, maxPhotonID_cut = -1.0;
      float minPhotonID = -1.0, maxPhotonID = 1.0;
      if ( selectedDiPhoton.leadPho.mvaID() > selectedDiPhoton.subleadPho.mvaID() ) {
        minPhotonID = selectedDiPhoton.subleadPho.mvaID();
        minPhotonID_cut = ( fabs(selectedDiPhoton.subleadPho.eta())<1.442 ? -0.02 : -0.26 );

        maxPhotonID = selectedDiPhoton.leadPho.mvaID();
        maxPhotonID_cut = ( fabs(selectedDiPhoton.leadPho.eta())<1.442 ? -0.02 : -0.26 );
      }
      else {
        minPhotonID = selectedDiPhoton.leadPho.mvaID();
        minPhotonID_cut = ( fabs(selectedDiPhoton.leadPho.eta())<1.442 ? -0.02 : -0.26 );

        maxPhotonID = selectedDiPhoton.subleadPho.mvaID();
        maxPhotonID_cut = ( fabs(selectedDiPhoton.subleadPho.eta())<1.442 ? -0.02 : -0.26 );
      }
      if ( process == "DDQCDGJets" ) {
        if ( !(maxPhotonID > maxPhotonID_cut && minPhotonID < minPhotonID_cut) ) continue; // Select events in which only one photon fails the photon MVA selection
        if ( maxPhotonID <= minPhotonID_cut ) continue; // Incompatible photon ID configuration, where the max and min ID photon would be swapped
        
        if ( selectedDiPhoton.leadPho.mvaID() > selectedDiPhoton.subleadPho.mvaID() ) {
          if ( fabs(selectedDiPhoton.subleadPho.eta()) < 1.5 ) imputePhotonID(minPhotonID_cut, maxPhotonID, fakePhotonID_shape_barrel, minPhotonID, weight);
          else imputePhotonID(minPhotonID_cut, maxPhotonID, fakePhotonID_shape_endcap, minPhotonID, weight);
          selectedDiPhoton.subleadPho.setMvaID(minPhotonID);
        }
        else {
          if ( fabs(selectedDiPhoton.leadPho.eta()) < 1.5 ) imputePhotonID(minPhotonID_cut, maxPhotonID, fakePhotonID_shape_barrel, minPhotonID, weight);
          else imputePhotonID(minPhotonID_cut, maxPhotonID, fakePhotonID_shape_endcap, minPhotonID, weight);
          selectedDiPhoton.leadPho.setMvaID(minPhotonID);
        }
      }

      // Select photons passing the photon MVA ID WP 90
      if ( !( (fabs(selectedDiPhoton.leadPho.eta())<1.442 ? selectedDiPhoton.leadPho.mvaID()>-0.02 : selectedDiPhoton.leadPho.mvaID()>-0.26) && \
              (fabs(selectedDiPhoton.subleadPho.eta())<1.442 ? selectedDiPhoton.subleadPho.mvaID()>-0.02 : selectedDiPhoton.subleadPho.mvaID()>-0.26) ) ) continue;

      Electrons electrons = getElectrons(selectedPhotons);
      Muons muons = getMuons(selectedPhotons);
      if ( enrichDY ) { if (electrons.size() == 0) continue; }
      else { if (electrons.size() != 0 ) continue; }
      if (muons.size() != 0 ) continue;

      Jets jets;
      if (isMC) jets = getJets(selectedPhotons, JESUnc, JERUnc, HEMCheck);
      else jets = getJets(selectedPhotons, 0, 0, 0);
      if (jets.size() < 2) continue; 

      DiJets dijets = DiJetPreselection(jets);
      DiJet selectedDiJet = dijets[0];

      if (dijets[0].p4.M()<50) continue;

      if (isMC) {
        // Apply electron veto SF
        if ( electronVetoSF!=0 ) {
          if ( electronVetoSF==1  ) weight *= electronVetoSF::get_electronVetoSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "central")*electronVetoSF::get_electronVetoSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "central");
          if ( electronVetoSF==2  ) weight *= electronVetoSF::get_electronVetoSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "up")*electronVetoSF::get_electronVetoSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "up");
          if ( electronVetoSF==-2 ) weight *= electronVetoSF::get_electronVetoSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "down")*electronVetoSF::get_electronVetoSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "down");
        }

        // Apply low mass trigger SF
        if ( triggerSF!=0 ) {
          if ( lowMassMode ) {
            bool leadEBHiR9=false, subleadEBHiR9=false;
            if ( year=="2016nonAPV" || year=="2016APV" ) {
              leadEBHiR9 = fabs(selectedDiPhoton.leadPho.eta()) < 1.5 && selectedDiPhoton.leadPho.r9() > 0.85 && fabs(selectedDiPhoton.subleadPho.eta()) < 1.5;
              subleadEBHiR9 = fabs(selectedDiPhoton.subleadPho.eta()) < 1.5 && selectedDiPhoton.subleadPho.r9() > 0.85 && fabs(selectedDiPhoton.leadPho.eta()) < 1.5;
            }
            if ( triggerSF==1  ) weight *= lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", leadEBHiR9, year, "central")*lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", subleadEBHiR9, year, "central");
            if ( triggerSF==2  ) weight *= lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", leadEBHiR9, year, "up")*lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", subleadEBHiR9, year, "up");
            if ( triggerSF==-2 ) weight *= lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", leadEBHiR9, year, "down")*lowMassHggTriggerSF::get_lowMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", subleadEBHiR9, year, "down");
          }
          else {
            if ( triggerSF==1  ) weight *= highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", year, "central")*highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", year, "central");
            if ( triggerSF==2  ) weight *= highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", year, "up")*highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", year, "up");
            if ( triggerSF==-2 ) weight *= highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), "Lead", year, "down")*highMassHggTriggerSF::get_highMassHggTriggerSF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), "Sublead", year, "down");
          }
        }

        // Apply low mass preselection SF
        if ( preselSF!=0 ) {
          if ( preselSF==1  ) weight *= lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "central")*lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "central");
          if ( preselSF==2  ) weight *= lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "up")*lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "up");
          if ( preselSF==-2 ) weight *= lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.leadPho.eta(), selectedDiPhoton.leadPho.r9(), year, "down")*lowMassHggPreselSF::get_lowMassHggPreselSF(selectedDiPhoton.subleadPho.eta(), selectedDiPhoton.subleadPho.r9(), year, "down");
        }

        // Apply photon MVA ID WP 90 SF
        if ( phoMVAIDWP90SF!=0 ) {
          if ( phoMVAIDWP90SF==1  ) weight *= phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), year, "central")*phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), year, "central");
          if ( phoMVAIDWP90SF==2  ) weight *= phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), year, "up")*phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), year, "up");
          if ( phoMVAIDWP90SF==-2 ) weight *= phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.leadPho.pt(), selectedDiPhoton.leadPho.eta(), year, "down")*phoMVAIDWP90SF::get_phoMVAIDWP90SF(selectedDiPhoton.subleadPho.pt(), selectedDiPhoton.subleadPho.eta(), year, "down");
        }

        // Apply bTagSF and get the weights before and after the application to normalize post-preselection. electedDiJet.leadJet. Links:
        // - https://btv-wiki.docs.cern.ch/PerformanceCalibration/SFUncertaintiesAndCorrelations/#working-point-based-sfs-fixedwp-sfs
        // - https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration#Correlation_across_years_2016_20
        // Jet_hadronFlavour = 0 for light jets, 4 for c jets, and 5 for b jets
        if ( bTagSF!=0 ) {
          float leadJetBTagSF, subleadJetBTagSF;
          if ( JESUnc!=0 ) { // btag JES correlated with jet JES no matter what.
            if ( JESUnc>1 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_jes();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_jes();
            }
            if ( JESUnc<-1 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_jes();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_jes();
            }
          }
          else {
            if ( bTagSF==1 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape();
            }
            if ( bTagSF==2 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_hf();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_hf();
            }
            if ( bTagSF==-2 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_hf();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_hf();
            }
            if ( bTagSF==3 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_hfstats1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_hfstats1();
            }
            if ( bTagSF==-3 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_hfstats1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_hfstats1();
            }
            if ( bTagSF==4 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_hfstats2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_hfstats2();
            }
            if ( bTagSF==-4 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_hfstats2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_hfstats2();
            }
            if ( bTagSF==5 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_lf();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_lf();
            }
            if ( bTagSF==-5 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_lf();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_lf();
            }
            if ( bTagSF==6 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_lfstats1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_lfstats1();
            }
            if ( bTagSF==-6 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_lfstats1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_lfstats1();
            }
            if ( bTagSF==7 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_lfstats2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_lfstats2();
            }
            if ( bTagSF==-7 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==5 || selectedDiJet.leadJet.hadronFlavour()==0 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_lfstats2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==5 || selectedDiJet.subleadJet.hadronFlavour()==0 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_lfstats2();
            }
            if ( bTagSF==8 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==4 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_cferr1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==4 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_cferr1();
            }
            if ( bTagSF==-8 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==4 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_cferr1();
              if ( selectedDiJet.subleadJet.hadronFlavour()==4 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_cferr1();
            }
            if ( bTagSF==9 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==4 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_up_cferr2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==4 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_up_cferr2();
            }
            if ( bTagSF==-9 ) {
              if ( selectedDiJet.leadJet.hadronFlavour()==4 ) leadJetBTagSF = selectedDiJet.leadJet.btagSF_deepjet_shape_down_cferr2();
              if ( selectedDiJet.subleadJet.hadronFlavour()==4 ) subleadJetBTagSF = selectedDiJet.subleadJet.btagSF_deepjet_shape_down_cferr2();
            }
          }
          weight_beforeBTagSF = weight;
          h_weight_beforeBTagSF->Fill(0.5, weight);
          weight *= leadJetBTagSF * subleadJetBTagSF;
          weight_afterBTagSF = weight;
          h_weight_afterBTagSF->Fill(0.5, weight);
        }
      }

      // Setting output variables
      TLorentzVector x_cand = selectedDiPhoton.p4 + selectedDiJet.p4;
      xcand_pt = x_cand.Pt();
      xcand_eta = x_cand.Eta();
      xcand_phi = x_cand.Phi();
      xcand_mass = x_cand.M();

      LeadPhoton_pt = selectedDiPhoton.leadPho.pt();
      LeadPhoton_eta = selectedDiPhoton.leadPho.eta();
      LeadPhoton_phi = selectedDiPhoton.leadPho.phi();
      LeadPhoton_r9 = selectedDiPhoton.leadPho.r9();
      LeadPhoton_pixelSeed = selectedDiPhoton.leadPho.pixelSeed();
      LeadPhoton_sieie = selectedDiPhoton.leadPho.sieie();
      LeadPhoton_pfPhoIso03 = selectedDiPhoton.leadPho.phoIso();
      LeadPhoton_trkSumPtHollowConeDR03 = selectedDiPhoton.leadPho.trkIso();
      LeadPhoton_chargedHadronIso = selectedDiPhoton.leadPho.chargedHadIso();
      LeadPhoton_mvaID = selectedDiPhoton.leadPho.mvaID();

      SubleadPhoton_pt = selectedDiPhoton.subleadPho.pt();
      SubleadPhoton_eta = selectedDiPhoton.subleadPho.eta();
      SubleadPhoton_phi = selectedDiPhoton.subleadPho.phi();
      SubleadPhoton_r9 = selectedDiPhoton.subleadPho.r9();
      SubleadPhoton_pixelSeed = selectedDiPhoton.subleadPho.pixelSeed();
      SubleadPhoton_sieie = selectedDiPhoton.subleadPho.sieie();
      SubleadPhoton_pfPhoIso03 = selectedDiPhoton.subleadPho.phoIso();
      SubleadPhoton_trkSumPtHollowConeDR03 = selectedDiPhoton.subleadPho.trkIso();
      SubleadPhoton_chargedHadronIso = selectedDiPhoton.subleadPho.chargedHadIso();
      SubleadPhoton_mvaID = selectedDiPhoton.subleadPho.mvaID();

      Diphoton_pt = selectedDiPhoton.p4.Pt();
      Diphoton_eta = selectedDiPhoton.p4.Eta();
      Diphoton_phi = selectedDiPhoton.p4.Phi();
      Diphoton_mass = selectedDiPhoton.p4.M();
      Diphoton_pt_mgg = Diphoton_pt/Diphoton_mass;
      Diphoton_dR = selectedDiPhoton.dR;

      n_jets = jets.size();
      dijet_lead_pt = selectedDiJet.leadJet.pt();
      dijet_lead_eta = selectedDiJet.leadJet.eta();
      dijet_lead_phi = selectedDiJet.leadJet.phi();
      dijet_lead_mass = selectedDiJet.leadJet.mass();
      dijet_lead_btagDeepFlavB = selectedDiJet.leadJet.btagDeepFlavB();
      dijet_sublead_pt = selectedDiJet.subleadJet.pt();
      dijet_sublead_eta = selectedDiJet.subleadJet.eta();
      dijet_sublead_phi = selectedDiJet.subleadJet.phi();
      dijet_sublead_mass = selectedDiJet.subleadJet.mass();
      dijet_sublead_btagDeepFlavB = selectedDiJet.subleadJet.btagDeepFlavB();
      dijet_pt = selectedDiJet.p4.Pt();
      dijet_eta = selectedDiJet.p4.Eta();
      dijet_phi = selectedDiJet.p4.Phi();
      dijet_mass = selectedDiJet.p4.M();
      dijet_dR = selectedDiJet.dR;
      pfmet_pt = isMC ? nt.MET_T1_pt() : nt.MET_pt();
      puppimet_pt = nt.PuppiMET_pt();

      weight_central = weight*factor;
      weight_central_initial = weight;
      weight_central_no_lumi = weight*factor/lumi;
      eventNum = evtnb;

      count_test++;

      if (isMC){
        LeadPhoton_genPartFlav = int(selectedDiPhoton.leadPho.genPartFlav());
        SubleadPhoton_genPartFlav = int(selectedDiPhoton.subleadPho.genPartFlav());
        TLorentzVector GenHiggs, GenX, GenY;
        GenParts genparts = getGenParts();
        vector<TLorentzVector> gen_child_xyh, gen_child_ygg, gen_child_hbb;
        for (int igenpart=0; igenpart<genparts.size(); igenpart++)
          {
            if (genparts[igenpart].isxyh()) {
              GenX = genparts[igenpart].mother_p4(); 
              gen_child_xyh.push_back(genparts[igenpart].p4());
              GenX_pt = GenX.Pt();
              GenX_eta = GenX.Eta();
              GenX_phi = GenX.Phi();
              GenX_mass = GenX.M();
            }
            if (genparts[igenpart].isygg()) {
              GenY = genparts[igenpart].mother_p4(); 
              gen_child_ygg.push_back(genparts[igenpart].p4());
              GenY_pt = GenY.Pt();
              GenY_eta = GenY.Eta();
              GenY_phi = GenY.Phi();
              GenY_mass = GenY.M();
            }
            if (genparts[igenpart].ishbb()) {
              GenHiggs = genparts[igenpart].mother_p4();
              gen_child_hbb.push_back(genparts[igenpart].p4());
              GenHiggs_pt = GenHiggs.Pt();
              GenHiggs_eta = GenHiggs.Eta();
              GenHiggs_phi = GenHiggs.Phi();
              GenHiggs_mass = GenHiggs.M();
            }
          }
        if (!abs(GenX_pt+999)<0.0001){
          GenX_dR = gen_child_xyh[0].DeltaR(gen_child_xyh[1]);
        }
        if (!abs(GenY_pt+999)<0.0001){
          GenY_dR = gen_child_ygg[0].DeltaR(gen_child_ygg[1]);
        }
        if (!abs(GenHiggs_pt+999)<0.0001){
          GenHiggs_dR = gen_child_hbb[0].DeltaR(gen_child_hbb[1]);
          TLorentzVector gen_child_hbb_temp;
          if (gen_child_hbb[0].Pt()<gen_child_hbb[1].Pt()) {gen_child_hbb_temp = gen_child_hbb[0]; gen_child_hbb[0]= gen_child_hbb[1]; gen_child_hbb[1] = gen_child_hbb_temp;}
          GenBFromHiggs_1_pt = gen_child_hbb[0].Pt();
          GenBFromHiggs_1_eta = gen_child_hbb[0].Eta();
          GenBFromHiggs_1_phi = gen_child_hbb[0].Phi();
          GenBFromHiggs_1_mass = gen_child_hbb[0].M();
          GenBFromHiggs_2_pt = gen_child_hbb[1].Pt();
          GenBFromHiggs_2_eta = gen_child_hbb[1].Eta();
          GenBFromHiggs_2_phi = gen_child_hbb[1].Phi();
          GenBFromHiggs_2_mass = gen_child_hbb[1].M();
          if (selectedDiJet.leadJet.p4().DeltaR(gen_child_hbb[0])<=0.4 || selectedDiJet.leadJet.p4().DeltaR(gen_child_hbb[1])<=0.4) {dijet_lead_gen_match=true; n_gen_matched_in_dijet++;}
          if (selectedDiJet.subleadJet.p4().DeltaR(gen_child_hbb[0])<=0.4 || selectedDiJet.subleadJet.p4().DeltaR(gen_child_hbb[1])<=0.4) {dijet_sublead_gen_match=true; n_gen_matched_in_dijet++;}
          for (int ijet=0; ijet<jets.size(); ijet++)
            if (jets[ijet].p4().DeltaR(gen_child_hbb[0])<=0.4 || jets[ijet].p4().DeltaR(gen_child_hbb[1])<=0.4)
              n_gen_matched_jets++;
        }
      }

      // Histo filling
      h_LeadPhoton_sieie->Fill(LeadPhoton_sieie);
      h_LeadPhoton_pfPhoIso03->Fill(LeadPhoton_pfPhoIso03);
      h_LeadPhoton_chargedHadronIso->Fill(LeadPhoton_chargedHadronIso);
      h_LeadPhoton_trkSumPtHollowConeDR03->Fill(LeadPhoton_trkSumPtHollowConeDR03);
      h_SubleadPhoton_sieie->Fill(SubleadPhoton_sieie);
      h_SubleadPhoton_pfPhoIso03->Fill(SubleadPhoton_pfPhoIso03);
      h_SubleadPhoton_chargedHadronIso->Fill(SubleadPhoton_chargedHadronIso);
      h_SubleadPhoton_trkSumPtHollowConeDR03->Fill(SubleadPhoton_trkSumPtHollowConeDR03);
      h_weight->Fill(0.5, weight*factor);

      tout->Fill();
    } // Event loop
    delete file;
  } // File loop
    
  bar.finish();
  if ( electronVetoSF != 0) electronVetoSF::reset_electronVetoSF();
  if ( triggerSF != 0) ( lowMassMode ? lowMassHggTriggerSF::reset_lowMassHggTriggerSF() : highMassHggTriggerSF::reset_highMassHggTriggerSF() );
  if ( preselSF != 0) lowMassHggPreselSF::reset_lowMassHggPreselSF();
  if ( phoMVAIDWP90SF != 0) phoMVAIDWP90SF::reset_phoMVAIDWP90SF();
  cout << "nTotal: " << h_weight_full->GetBinContent(1) << ", nPass: " << h_weight->GetBinContent(1) << ", eff: " << h_weight->GetBinContent(1)/h_weight_full->GetBinContent(1) << endl;
  cout << endl;

  if ( removeDataDuplicates )
    cout << "Number of duplicates found: " << nDuplicates << endl;

  txtout.close();
  fout->Write();
  fout->Close();

  return 0;
}
