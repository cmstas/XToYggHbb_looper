#ifndef LOWMASSHGGTRIGGERSF_H
#define LOWMASSHGGTRIGGERSF_H

#include <map>
#include <stdexcept>

using namespace std;

namespace lowMassHggTriggerSF {
  map<TString,map<TString,map<TString,map<TString,float>>>> idsf_lead = { };
  map<TString,map<TString,map<TString,map<TString,float>>>> idsfunc_lead = { };
  map<TString,map<TString,map<TString,map<TString,float>>>> idsf_sublead = { };
  map<TString,map<TString,map<TString,map<TString,float>>>> idsfunc_sublead = { };

  constexpr int netabins = 3;

  float *thresholds_eta = new float[netabins];
  TString *etabins = new TString[netabins];
  inline void set_etabins() {
    thresholds_eta[0] = 0.0;
    etabins[0] = "eta0";
    thresholds_eta[1] = 1.5;
    etabins[1] = "eta1";
    thresholds_eta[2] = 3.0;
    etabins[2] = "eta2";
  }

  TString get_etaBin(const float abseta) {
    for ( unsigned int b=netabins-1; b>=1; b-- ) {
      if ( abseta > thresholds_eta[b] ) {
        return etabins[b];
      }
    }
    return etabins[0];
  }

  constexpr int nr9bins = 5;

  float *thresholds_r9 = new float[nr9bins];
  TString *r9bins = new TString[nr9bins];
  inline void set_r9bins(const float abseta, const TString photon, const bool EBHiR9, const TString year) {
    if ( abseta <= 1.5 ) {
      if ( year == "2018" ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.5;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.54;
        r9bins[2] = "r92";
        thresholds_r9[3] = 0.85; // The following bins are dummy
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0;
        r9bins[4] = "r94";
      }
      else if ( year == "2017" && photon == "Lead" ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.85;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.88;
        r9bins[2] = "r92";
        thresholds_r9[3] = 999.0; // The following bins are dummy
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0;
        r9bins[4] = "r94";
      }
      else if ( year == "2017" && photon == "Sublead" ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.5;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.53;
        r9bins[2] = "r92";
        thresholds_r9[3] = 0.85;
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0; // The following bins are dummy
        r9bins[4] = "r94";
      }
      else if ( year == "2016" && EBHiR9 ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.85;
        r9bins[1] = "r91";
        thresholds_r9[2] = 999.0; // The following bins are dummy
        r9bins[2] = "r92";
        thresholds_r9[3] = 999.0;
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0;
        r9bins[4] = "r94";
      }
      else { // 2016
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.5;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.53;
        r9bins[2] = "r92";
        thresholds_r9[3] = 0.85;
        r9bins[3] = "r93";
        thresholds_r9[4] = 0.88;
        r9bins[4] = "r94";
      }
    }
    else if ( abseta <= 3.0 ) {
      if ( year == "2017" && photon == "Sublead" ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.8;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.83;
        r9bins[2] = "r92";
        thresholds_r9[3] = 0.9;
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0; // The following bins are dummy
        r9bins[4] = "r94";
      }
      else if ( year == "2016" && EBHiR9 ) {
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 999.0; // The following bins are dummy
        r9bins[1] = "r91";
        thresholds_r9[2] = 999.0;
        r9bins[2] = "r92";
        thresholds_r9[3] = 999.0;
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0;
        r9bins[4] = "r94";
      }
      else { // 2018, 2017Lead, 2016
        thresholds_r9[0] = 0.0;
        r9bins[0] = "r90";
        thresholds_r9[1] = 0.9;
        r9bins[1] = "r91";
        thresholds_r9[2] = 0.93;
        r9bins[2] = "r92";
        thresholds_r9[3] = 999.0; // The following bins are dummy
        r9bins[3] = "r93";
        thresholds_r9[4] = 999.0;
        r9bins[4] = "r94";
      }
    }
    else {
      thresholds_r9[0] = 0.0;
      r9bins[0] = "r90";
      thresholds_r9[1] = 999.0; // The following bins are dummy
      r9bins[1] = "r91";
      thresholds_r9[2] = 999.0;
      r9bins[2] = "r92";
      thresholds_r9[3] = 999.0;
      r9bins[3] = "r93";
      thresholds_r9[4] = 999.0;
      r9bins[4] = "r94";
    }
  }

  TString get_r9Bin(const float r9) {
    for ( unsigned int b=nr9bins-1; b>=1; b-- ) {
      if ( r9 > thresholds_r9[b] ) {
        return r9bins[b];
      }
    }
    return r9bins[0];
  }

  constexpr int nptbins = 12;

  float *thresholds_pt = new float[nptbins];
  TString *ptbins = new TString[nptbins];
  inline void set_ptbins(const float abseta, const float r9, const TString photon, const bool EBHiR9, const TString year) {
    if ( abseta > 3.0 ) {
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else if ( abseta > 1.5 && EBHiR9 ) {
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else if ( (year == "2017" && photon == "Sublead") && abseta > 1.5 && r9 <= 0.8 ) {
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else if ( !(year == "2017" && photon == "Sublead")  && abseta > 1.5 && r9 <= 0.9 ) { // 2016, 2017Lead, 2018
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else if ( (year == "2017" && photon == "Lead") && r9 <= 0.85 ) {
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else if ( !(year == "2017" && photon == "Lead") && r9 <= 0.5 ) { // 2016, 2017Sublead, 2018
      thresholds_pt[0] = 0.0;
      ptbins[0] = "pt0";
      thresholds_pt[1] = 999999.0; // The following bins are dummy
      ptbins[1] = "pt1";
      thresholds_pt[2] = 999999.0;
      ptbins[2] = "pt2";
      thresholds_pt[3] = 999999.0;
      ptbins[3] = "pt3";
      thresholds_pt[4] = 999999.0;
      ptbins[4] = "pt4";
      thresholds_pt[5] = 999999.0;
      ptbins[5] = "pt5";
      thresholds_pt[6] = 999999.0;
      ptbins[6] = "pt6";
      thresholds_pt[7] = 999999.0;
      ptbins[7] = "pt7";
      thresholds_pt[8] = 999999.0;
      ptbins[8] = "pt8";
      thresholds_pt[9] = 999999.0;
      ptbins[9] = "pt9";
      thresholds_pt[10] = 999999.0;
      ptbins[10] = "pt10";
      thresholds_pt[11] = 999999.0;
      ptbins[11] = "pt11";
    }
    else {
      if ( photon == "Lead" ) {
        if ( year == "2018" ) {
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 35.0;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 37.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 40.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 45.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 50.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 60.0;
          ptbins[6] = "pt6";
          thresholds_pt[7] = 70.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 90.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 999999.0; // The following bins are dummy
          ptbins[9] = "pt9";
          thresholds_pt[10] = 999999.0;
          ptbins[10] = "pt10";
          thresholds_pt[11] = 999999.0;
          ptbins[11] = "pt11";
        }
        else if ( year == "2017" ) {
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 33.0;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 35.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 37.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 40.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 45.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 50.0;
          ptbins[6] = "pt6";
          thresholds_pt[7] = 60.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 70.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 90.0;
          ptbins[9] = "pt9";
          thresholds_pt[10] = 999999.0; // The following bins are dummy
          ptbins[10] = "pt10";
          thresholds_pt[11] = 999999.0;
          ptbins[11] = "pt11";
        }
        else { // 2016
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 33.3333;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 35.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 40.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 45.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 50.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 60.0;
          ptbins[6] = "pt6";
          thresholds_pt[7] = 70.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 90.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 999999.0; // The following bins are dummy
          ptbins[9] = "pt9";
          thresholds_pt[10] = 999999.0;
          ptbins[10] = "pt10";
          thresholds_pt[11] = 999999.0;
          ptbins[11] = "pt11";
        }
      }
      else if ( photon == "Sublead" ) {
        if ( year == "2018" ) {
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 25.0;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 28.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 31.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 35.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 40.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 45.0;
          ptbins[6] = "pt6";
          thresholds_pt[7] = 50.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 60.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 70.0;
          ptbins[9] = "pt9";
          thresholds_pt[10] = 90.0;
          ptbins[10] = "pt10";
          thresholds_pt[11] = 999999.0; // The following bins are dummy
          ptbins[11] = "pt11";
        }
        else if ( year == "2017" ) {
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 25.0;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 28.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 31.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 40.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 70.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 999999.0; // The following bins are dummy
          ptbins[6] = "pt6";
          thresholds_pt[7] = 999999.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 999999.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 999999.0;
          ptbins[9] = "pt9";
          thresholds_pt[10] = 999999.0;
          ptbins[10] = "pt10";
          thresholds_pt[11] = 999999.0;
          ptbins[11] = "pt11";
        }
        else { // 2016
          thresholds_pt[0] = 0.0;
          ptbins[0] = "pt0";
          thresholds_pt[1] = 32.0;
          ptbins[1] = "pt1";
          thresholds_pt[2] = 34.0;
          ptbins[2] = "pt2";
          thresholds_pt[3] = 36.0;
          ptbins[3] = "pt3";
          thresholds_pt[4] = 38.0;
          ptbins[4] = "pt4";
          thresholds_pt[5] = 40.0;
          ptbins[5] = "pt5";
          thresholds_pt[6] = 45.0;
          ptbins[6] = "pt6";
          thresholds_pt[7] = 50.0;
          ptbins[7] = "pt7";
          thresholds_pt[8] = 55.0;
          ptbins[8] = "pt8";
          thresholds_pt[9] = 65.0;
          ptbins[9] = "pt9";
          thresholds_pt[10] = 70.0;
          ptbins[10] = "pt10";
          thresholds_pt[11] = 90.0;
          ptbins[11] = "pt11";
        }
      }
      else {
        throw std::invalid_argument( "Argument can only be 'Lead' or 'Sublead'." );
      }
    }
  }

  TString get_ptBin(const float pt) {
    for ( unsigned int b=nptbins-1; b>=1; b-- ) {
      if ( pt > thresholds_pt[b] ) {
        return ptbins[b];
      }
    }
    return ptbins[0];
  }

  inline void reset_lowMassHggTriggerSF() {
    idsf_lead.clear();
    idsfunc_lead.clear();
    idsf_sublead.clear();
    idsfunc_sublead.clear();
  }

  inline void set_lowMassHggTriggerSF() {
    idsf_lead.insert({"2016", { }});
    idsfunc_lead.insert({"2016", { }});
      idsf_lead["2016"].insert({"eta0", { }});
      idsfunc_lead["2016"].insert({"eta0", { }});
        idsf_lead["2016"]["eta0"].insert({"r90", { }});
        idsfunc_lead["2016"]["eta0"].insert({"r90", { }});
          idsf_lead["2016"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2016"]["eta0"].insert({"r91", { }});
        idsfunc_lead["2016"]["eta0"].insert({"r91", { }});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt0", 0.594793000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt0", 0.006127890});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt1", 0.662945000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt1", 0.007832360});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt2", 0.692324000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt2", 0.003986420});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt3", 0.708503000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt3", 0.003432880});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt4", 0.732292000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt4", 0.004387480});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt5", 0.744355000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt5", 0.006468870});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt6", 0.792566000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt6", 0.015080600});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt7", 0.808000000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt7", 0.028349700});
          idsf_lead["2016"]["eta0"]["r91"].insert({"pt8", 0.783784000});
          idsfunc_lead["2016"]["eta0"]["r91"].insert({"pt8", 0.089627400});
        idsf_lead["2016"]["eta0"].insert({"r92", { }});
        idsfunc_lead["2016"]["eta0"].insert({"r92", { }});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt0", 0.880773000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt0", 0.001000070});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt1", 0.951792000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt1", 0.000842411});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt2", 0.961186000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt2", 0.000371599});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt3", 0.967805000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt3", 0.000267569});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt4", 0.971748000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt4", 0.000303140});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt5", 0.974850000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt5", 0.000389931});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt6", 0.979271000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt6", 0.000769622});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt7", 0.983438000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt7", 0.001079670});
          idsf_lead["2016"]["eta0"]["r92"].insert({"pt8", 0.990463000});
          idsfunc_lead["2016"]["eta0"]["r92"].insert({"pt8", 0.001437670});
        idsf_lead["2016"]["eta0"].insert({"r93", { }});
        idsfunc_lead["2016"]["eta0"].insert({"r93", { }});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt0", 0.480261000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt0", 0.004790900});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt1", 0.544853000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt1", 0.005953340});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt2", 0.564528000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt2", 0.002857180});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt3", 0.584602000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt3", 0.002182620});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt4", 0.610664000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt4", 0.002551310});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt5", 0.638541000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt5", 0.003218980});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt6", 0.679178000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt6", 0.006107500});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt7", 0.719895000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt7", 0.007941300});
          idsf_lead["2016"]["eta0"]["r93"].insert({"pt8", 0.782112000});
          idsfunc_lead["2016"]["eta0"]["r93"].insert({"pt8", 0.009429590});
        idsf_lead["2016"]["eta0"].insert({"r94", { }});
        idsfunc_lead["2016"]["eta0"].insert({"r94", { }});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt0", 0.911612000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt0", 0.000551594});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt1", 0.964408000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt1", 0.000451792});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt2", 0.973082000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt2", 0.000187761});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt3", 0.978494000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt3", 0.000126388});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt4", 0.981905000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt4", 0.000137401});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt5", 0.982553000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt5", 0.000174033});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt6", 0.983744000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt6", 0.000329198});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt7", 0.986766000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt7", 0.000405361});
          idsf_lead["2016"]["eta0"]["r94"].insert({"pt8", 0.990753000});
          idsfunc_lead["2016"]["eta0"]["r94"].insert({"pt8", 0.000407389});
      idsf_lead["2016"].insert({"eta1", { }});
      idsfunc_lead["2016"].insert({"eta1", { }});
        idsf_lead["2016"]["eta1"].insert({"r90", { }});
        idsfunc_lead["2016"]["eta1"].insert({"r90", { }});
          idsf_lead["2016"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2016"]["eta1"].insert({"r91", { }});
        idsfunc_lead["2016"]["eta1"].insert({"r91", { }});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt0", 0.595211000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt0", 0.005601170});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt1", 0.702697000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt1", 0.006740310});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt2", 0.734526000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt2", 0.003284520});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt3", 0.760539000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt3", 0.002529220});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt4", 0.778439000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt4", 0.003045640});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt5", 0.791388000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt5", 0.003802740});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt6", 0.814570000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt6", 0.007199710});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt7", 0.837738000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt7", 0.008972370});
          idsf_lead["2016"]["eta1"]["r91"].insert({"pt8", 0.884554000});
          idsfunc_lead["2016"]["eta1"]["r91"].insert({"pt8", 0.009787100});
        idsf_lead["2016"]["eta1"].insert({"r92", { }});
        idsfunc_lead["2016"]["eta1"].insert({"r92", { }});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt0", 0.853852000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt0", 0.001546490});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt1", 0.950491000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt1", 0.001278020});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt2", 0.965209000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt2", 0.000558641});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt3", 0.970964000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt3", 0.000416357});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt4", 0.973743000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt4", 0.000499808});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt5", 0.977689000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt5", 0.000597825});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt6", 0.981832000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt6", 0.001076080});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt7", 0.989488000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt7", 0.001100480});
          idsf_lead["2016"]["eta1"]["r92"].insert({"pt8", 0.993663000});
          idsfunc_lead["2016"]["eta1"]["r92"].insert({"pt8", 0.001054880});
      idsf_lead["2016"].insert({"eta2", { }});
      idsfunc_lead["2016"].insert({"eta2", { }});
        idsf_lead["2016"]["eta2"].insert({"r90", { }});
        idsfunc_lead["2016"]["eta2"].insert({"r90", { }});
          idsf_lead["2016"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016"]["eta2"]["r90"].insert({"pt0", 0.0});

    idsf_lead.insert({"2016EBHiR9", { }});
    idsfunc_lead.insert({"2016EBHiR9", { }});
      idsf_lead["2016EBHiR9"].insert({"eta0", { }});
      idsfunc_lead["2016EBHiR9"].insert({"eta0", { }});
        idsf_lead["2016EBHiR9"]["eta0"].insert({"r90", { }});
        idsfunc_lead["2016EBHiR9"]["eta0"].insert({"r90", { }});
          idsf_lead["2016EBHiR9"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2016EBHiR9"]["eta0"].insert({"r91", { }});
        idsfunc_lead["2016EBHiR9"]["eta0"].insert({"r91", { }});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt0", 0.917161000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt0", 0.000453874});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt1", 0.969697000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt1", 0.000353493});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt2", 0.977014000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt2", 0.000146821});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt3", 0.980632000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt3", 0.000100407});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt4", 0.983649000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt4", 0.000108830});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt5", 0.984346000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt5", 0.000136420});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt6", 0.985730000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt6", 0.000254647});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt7", 0.989070000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt7", 0.000308570});
          idsf_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt8", 0.992703000});
          idsfunc_lead["2016EBHiR9"]["eta0"]["r91"].insert({"pt8", 0.000311862});
      idsf_lead["2016EBHiR9"].insert({"eta1", { }});
      idsfunc_lead["2016EBHiR9"].insert({"eta1", { }});
        idsf_lead["2016EBHiR9"]["eta1"].insert({"r90", { }});
        idsfunc_lead["2016EBHiR9"]["eta1"].insert({"r90", { }});
          idsf_lead["2016EBHiR9"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016EBHiR9"]["eta1"]["r90"].insert({"pt0", 0.0});
      idsf_lead["2016EBHiR9"].insert({"eta2", { }});
      idsfunc_lead["2016EBHiR9"].insert({"eta2", { }});
        idsf_lead["2016EBHiR9"]["eta2"].insert({"r90", { }});
        idsfunc_lead["2016EBHiR9"]["eta2"].insert({"r90", { }});
          idsf_lead["2016EBHiR9"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2016EBHiR9"]["eta2"]["r90"].insert({"pt0", 0.0});

    idsf_sublead.insert({"2016", { }});
    idsfunc_sublead.insert({"2016", { }});
      idsf_sublead["2016"].insert({"eta0", { }});
      idsfunc_sublead["2016"].insert({"eta0", { }});
        idsf_sublead["2016"]["eta0"].insert({"r90", { }});
        idsfunc_sublead["2016"]["eta0"].insert({"r90", { }});
          idsf_sublead["2016"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2016"]["eta0"].insert({"r91", { }});
        idsfunc_sublead["2016"]["eta0"].insert({"r91", { }});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt0", 0.641895000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt0", 0.005655140});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt1", 0.682347000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt1", 0.007698890});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt2", 0.695856000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt2", 0.007184140});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt3", 0.695441000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt3", 0.006825560});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt4", 0.714604000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt4", 0.006521380});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt5", 0.719518000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt5", 0.004430170});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt6", 0.743710000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt6", 0.005164890});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt7", 0.754574000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt7", 0.007924110});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt8", 0.763377000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt8", 0.010285900});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt9", 0.796185000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt9", 0.027643500});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt10", 0.802318000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt10", 0.027871300});
          idsf_sublead["2016"]["eta0"]["r91"].insert({"pt11", 0.769673000});
          idsfunc_sublead["2016"]["eta0"]["r91"].insert({"pt11", 0.088070700});
        idsf_sublead["2016"]["eta0"].insert({"r92", { }});
        idsfunc_sublead["2016"]["eta0"].insert({"r92", { }});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt0", 0.963248000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt0", 0.005819560});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt1", 0.970865000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt1", 0.004565600});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt2", 0.973962000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt2", 0.004128940});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt3", 0.975698000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt3", 0.004075500});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt4", 0.976490000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt4", 0.004026880});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt5", 0.978155000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt5", 0.004048440});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt6", 0.978768000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt6", 0.004000310});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt7", 0.979168000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt7", 0.003997050});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt8", 0.979229000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt8", 0.003999800});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt9", 0.978977000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt9", 0.004035520});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt10", 0.978874000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt10", 0.004047160});
          idsf_sublead["2016"]["eta0"]["r92"].insert({"pt11", 0.979177000});
          idsfunc_sublead["2016"]["eta0"]["r92"].insert({"pt11", 0.004121280});
        idsf_sublead["2016"]["eta0"].insert({"r93", { }});
        idsfunc_sublead["2016"]["eta0"].insert({"r93", { }});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt0", 0.520538000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt0", 0.004318930});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt1", 0.544783000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt1", 0.006225360});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt2", 0.552613000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt2", 0.005689080});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt3", 0.569686000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt3", 0.005211650});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt4", 0.575110000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt4", 0.004813300});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt5", 0.586177000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt5", 0.003234180});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt6", 0.612741000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt6", 0.003568640});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt7", 0.635098000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt7", 0.004720430});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt8", 0.656744000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt8", 0.005286250});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt9", 0.686169000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt9", 0.010787500});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt10", 0.716083000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt10", 0.008521880});
          idsf_sublead["2016"]["eta0"]["r93"].insert({"pt11", 0.771437000});
          idsfunc_sublead["2016"]["eta0"]["r93"].insert({"pt11", 0.010028300});
        idsf_sublead["2016"]["eta0"].insert({"r94", { }});
        idsfunc_sublead["2016"]["eta0"].insert({"r94", { }});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt0", 0.957582000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt0", 0.004087330});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt1", 0.969330000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt1", 0.003986260});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt2", 0.971038000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt2", 0.003953990});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt3", 0.973016000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt3", 0.003954450});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt4", 0.975327000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt4", 0.003958410});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt5", 0.977557000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt5", 0.003965220});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt6", 0.978469000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt6", 0.003966840});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt7", 0.979083000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt7", 0.003970480});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt8", 0.979269000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt8", 0.003972170});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt9", 0.979605000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt9", 0.003988670});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt10", 0.978773000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt10", 0.003980140});
          idsf_sublead["2016"]["eta0"]["r94"].insert({"pt11", 0.980216000});
          idsfunc_sublead["2016"]["eta0"]["r94"].insert({"pt11", 0.003988400});
      idsf_sublead["2016"].insert({"eta1", { }});
      idsfunc_sublead["2016"].insert({"eta1", { }});
        idsf_sublead["2016"]["eta1"].insert({"r90", { }});
        idsfunc_sublead["2016"]["eta1"].insert({"r90", { }});
          idsf_sublead["2016"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2016"]["eta1"].insert({"r91", { }});
        idsfunc_sublead["2016"]["eta1"].insert({"r91", { }});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt0", 0.685109000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt0", 0.014049000});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt1", 0.697586000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt1", 0.015049100});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt2", 0.708077000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt2", 0.015013700});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt3", 0.714550000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt3", 0.014968400});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt4", 0.728204000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt4", 0.015043100});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt5", 0.745900000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt5", 0.014774000});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt6", 0.761913000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt6", 0.015180200});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt7", 0.769585000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt7", 0.015741700});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt8", 0.784002000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt8", 0.016251100});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt9", 0.796937000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt9", 0.019829800});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt10", 0.815029000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt10", 0.018342300});
          idsf_sublead["2016"]["eta1"]["r91"].insert({"pt11", 0.858506000});
          idsfunc_sublead["2016"]["eta1"]["r91"].insert({"pt11", 0.019501300});
        idsf_sublead["2016"]["eta1"].insert({"r92", { }});
        idsfunc_sublead["2016"]["eta1"].insert({"r92", { }});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt0", 0.952053000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt0", 0.018743300});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt1", 0.960736000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt1", 0.018768800});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt2", 0.963011000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt2", 0.018795100});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt3", 0.962529000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt3", 0.018784900});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt4", 0.965013000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt4", 0.018830500});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt5", 0.965534000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt5", 0.018837100});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt6", 0.966361000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt6", 0.018853200});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt7", 0.965761000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt7", 0.018845200});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt8", 0.966527000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt8", 0.018861900});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt9", 0.964668000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt9", 0.018870900});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt10", 0.966477000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt10", 0.018878600});
          idsf_sublead["2016"]["eta1"]["r92"].insert({"pt11", 0.966970000});
          idsfunc_sublead["2016"]["eta1"]["r92"].insert({"pt11", 0.018897900});
      idsf_sublead["2016"].insert({"eta2", { }});
      idsfunc_sublead["2016"].insert({"eta2", { }});
        idsf_sublead["2016"]["eta2"].insert({"r90", { }});
        idsfunc_sublead["2016"]["eta2"].insert({"r90", { }});
          idsf_sublead["2016"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016"]["eta2"]["r90"].insert({"pt0", 0.0});

    idsf_sublead.insert({"2016EBHiR9", { }});
    idsfunc_sublead.insert({"2016EBHiR9", { }});
      idsf_sublead["2016EBHiR9"].insert({"eta0", { }});
      idsfunc_sublead["2016EBHiR9"].insert({"eta0", { }});
        idsf_sublead["2016EBHiR9"]["eta0"].insert({"r90", { }});
        idsfunc_sublead["2016EBHiR9"]["eta0"].insert({"r90", { }});
          idsf_sublead["2016EBHiR9"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2016EBHiR9"]["eta0"].insert({"r91", { }});
        idsfunc_sublead["2016EBHiR9"]["eta0"].insert({"r91", { }});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt0", 0.973687000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt0", 0.005079770});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt1", 0.978161000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt1", 0.004996890});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt2", 0.979632000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt2", 0.004981490});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt3", 0.980581000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt3", 0.004981970});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt4", 0.981697000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt4", 0.004984960});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt5", 0.982505000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt5", 0.004989110});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt6", 0.983331000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt6", 0.004991810});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt7", 0.984008000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt7", 0.004995070});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt8", 0.984377000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt8", 0.004996930});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt9", 0.984618000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt9", 0.004998240});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt10", 0.984768000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt10", 0.004999020});
          idsf_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt11", 0.984870000});
          idsfunc_sublead["2016EBHiR9"]["eta0"]["r91"].insert({"pt11", 0.004999510});
      idsf_sublead["2016EBHiR9"].insert({"eta1", { }});
      idsfunc_sublead["2016EBHiR9"].insert({"eta1", { }});
        idsf_sublead["2016EBHiR9"]["eta1"].insert({"r90", { }});
        idsfunc_sublead["2016EBHiR9"]["eta1"].insert({"r90", { }});
          idsf_sublead["2016EBHiR9"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016EBHiR9"]["eta1"]["r90"].insert({"pt0", 0.0});
      idsf_sublead["2016EBHiR9"].insert({"eta2", { }});
      idsfunc_sublead["2016EBHiR9"].insert({"eta2", { }});
        idsf_sublead["2016EBHiR9"]["eta2"].insert({"r90", { }});
        idsfunc_sublead["2016EBHiR9"]["eta2"].insert({"r90", { }});
          idsf_sublead["2016EBHiR9"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2016EBHiR9"]["eta2"]["r90"].insert({"pt0", 0.0});


    idsf_lead.insert({"2017", { }});
    idsfunc_lead.insert({"2017", { }});
      idsf_lead["2017"].insert({"eta0", { }});
      idsfunc_lead["2017"].insert({"eta0", { }});
        idsf_lead["2017"]["eta0"].insert({"r90", { }});
        idsfunc_lead["2017"]["eta0"].insert({"r90", { }});
          idsf_lead["2017"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2017"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2017"]["eta0"].insert({"r91", { }});
        idsfunc_lead["2017"]["eta0"].insert({"r91", { }});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt0", 0.697427000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt0", 0.001378320});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt1", 0.882715000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt1", 0.001075570});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt2", 0.894138000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt2", 0.001007880});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt3", 0.902827000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt3", 0.001008490});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt4", 0.912114000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt4", 0.001009140});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt5", 0.915718000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt5", 0.001009390});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt6", 0.920406000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt6", 0.001009720});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt7", 0.924396000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt7", 0.001339350});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt8", 0.928790000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt8", 0.001556060});
          idsf_lead["2017"]["eta0"]["r91"].insert({"pt9", 0.939501000});
          idsfunc_lead["2017"]["eta0"]["r91"].insert({"pt9", 0.001672750});
        idsf_lead["2017"]["eta0"].insert({"r92", { }});
        idsfunc_lead["2017"]["eta0"].insert({"r92", { }});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt0", 0.764902000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt0", 0.000999564});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt1", 0.930524000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt1", 0.001010450});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt2", 0.938207000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt2", 0.001011000});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt3", 0.943564000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt3", 0.001011390});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt4", 0.947708000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt4", 0.001011690});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt5", 0.949471000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt5", 0.001011820});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt6", 0.949001000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt6", 0.001011780});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt7", 0.950820000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt7", 0.001011920});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt8", 0.957201000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt8", 0.001012390});
          idsf_lead["2017"]["eta0"]["r92"].insert({"pt9", 0.963203000});
          idsfunc_lead["2017"]["eta0"]["r92"].insert({"pt9", 0.001012830});
      idsf_lead["2017"].insert({"eta1", { }});
      idsfunc_lead["2017"].insert({"eta1", { }});
        idsf_lead["2017"]["eta1"].insert({"r90", { }});
        idsfunc_lead["2017"]["eta1"].insert({"r90", { }});
          idsf_lead["2017"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2017"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2017"]["eta1"].insert({"r91", { }});
        idsfunc_lead["2017"]["eta1"].insert({"r91", { }});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt0", 0.533179000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt0", 0.001663490});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt1", 0.781670000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt1", 0.001519800});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt2", 0.853343000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt2", 0.001217370});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt3", 0.882869000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt3", 0.001049660});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt4", 0.901054000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt4", 0.001052390});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt5", 0.914075000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt5", 0.001054380});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt6", 0.921084000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt6", 0.001055460});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt7", 0.933343000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt7", 0.001613840});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt8", 0.943441000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt8", 0.001773070});
          idsf_lead["2017"]["eta1"]["r91"].insert({"pt9", 0.958013000});
          idsfunc_lead["2017"]["eta1"]["r91"].insert({"pt9", 0.001830420});
        idsf_lead["2017"]["eta1"].insert({"r92", { }});
        idsfunc_lead["2017"]["eta1"].insert({"r92", { }});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt0", 0.587245000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt0", 0.001296540});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt1", 0.860594000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt1", 0.001047770});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt2", 0.904747000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt2", 0.001052960});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt3", 0.926744000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt3", 0.001056340});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt4", 0.942814000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt4", 0.001058860});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt5", 0.951324000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt5", 0.001060210});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt6", 0.953806000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt6", 0.001060600});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt7", 0.958042000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt7", 0.001068090});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt8", 0.964110000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt8", 0.001083310});
          idsf_lead["2017"]["eta1"]["r92"].insert({"pt9", 0.966957000});
          idsfunc_lead["2017"]["eta1"]["r92"].insert({"pt9", 0.001232760});
      idsf_lead["2017"].insert({"eta2", { }});
      idsfunc_lead["2017"].insert({"eta2", { }});
        idsf_lead["2017"]["eta2"].insert({"r90", { }});
        idsfunc_lead["2017"]["eta2"].insert({"r90", { }});
          idsf_lead["2017"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2017"]["eta2"]["r90"].insert({"pt0", 0.0});

    idsf_sublead.insert({"2017", { }});
    idsfunc_sublead.insert({"2017", { }});
      idsf_sublead["2017"].insert({"eta0", { }});
      idsfunc_sublead["2017"].insert({"eta0", { }});
        idsf_sublead["2017"]["eta0"].insert({"r90", { }});
        idsfunc_sublead["2017"]["eta0"].insert({"r90", { }});
          idsf_sublead["2017"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2017"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2017"]["eta0"].insert({"r91", { }});
        idsfunc_sublead["2017"]["eta0"].insert({"r91", { }});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt0", 0.762829000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt0", 0.027059500});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt1", 0.881282000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt1", 0.043552900});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt2", 0.926675000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt2", 0.041405600});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt3", 0.902538000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt3", 0.022605200});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt4", 0.881184000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt4", 0.025446900});
          idsf_sublead["2017"]["eta0"]["r91"].insert({"pt5", 0.821656000});
          idsfunc_sublead["2017"]["eta0"]["r91"].insert({"pt5", 0.178344000});
        idsf_sublead["2017"]["eta0"].insert({"r92", { }});
        idsfunc_sublead["2017"]["eta0"].insert({"r92", { }});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt0", 0.841898000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt0", 0.011523100});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt1", 0.925735000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt1", 0.014695000});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt2", 0.929883000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt2", 0.013892700});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt3", 0.935106000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt3", 0.011281100});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt4", 0.940301000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt4", 0.011331300});
          idsf_sublead["2017"]["eta0"]["r92"].insert({"pt5", 0.938715000});
          idsfunc_sublead["2017"]["eta0"]["r92"].insert({"pt5", 0.025564900});
        idsf_sublead["2017"]["eta0"].insert({"r93", { }});
        idsfunc_sublead["2017"]["eta0"].insert({"r93", { }});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt0", 0.931975000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt0", 0.015274700});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt1", 0.959115000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt1", 0.017123400});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt2", 0.968810000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt2", 0.012925400});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt3", 0.975143000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt3", 0.011415200});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt4", 0.970157000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt4", 0.011705000});
          idsf_sublead["2017"]["eta0"]["r93"].insert({"pt5", 0.975143000});
          idsfunc_sublead["2017"]["eta0"]["r93"].insert({"pt5", 0.011415200});
      idsf_sublead["2017"].insert({"eta1", { }});
      idsfunc_sublead["2017"].insert({"eta1", { }});
        idsf_sublead["2017"]["eta1"].insert({"r90", { }});
        idsfunc_sublead["2017"]["eta1"].insert({"r90", { }});
          idsf_sublead["2017"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2017"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2017"]["eta1"].insert({"r91", { }});
        idsfunc_sublead["2017"]["eta1"].insert({"r91", { }});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt0", 0.709933000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt0", 0.030332300});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt1", 0.815994000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt1", 0.045079200});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt2", 0.832926000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt2", 0.038794700});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt3", 0.850107000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt3", 0.021142000});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt4", 0.882189000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt4", 0.024057300});
          idsf_sublead["2017"]["eta1"]["r91"].insert({"pt5", 0.695921000});
          idsfunc_sublead["2017"]["eta1"]["r91"].insert({"pt5", 0.086399500});
        idsf_sublead["2017"]["eta1"].insert({"r92", { }});
        idsfunc_sublead["2017"]["eta1"].insert({"r92", { }});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt0", 0.803149000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt0", 0.018451300});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt1", 0.934628000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt1", 0.024966200});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt2", 0.932669000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt2", 0.022900700});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt3", 0.939581000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt3", 0.014787700});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt4", 0.933654000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt4", 0.014820300});
          idsf_sublead["2017"]["eta1"]["r92"].insert({"pt5", 0.965191000});
          idsfunc_sublead["2017"]["eta1"]["r92"].insert({"pt5", 0.034809000});
        idsf_sublead["2017"]["eta1"].insert({"r93", { }});
        idsfunc_sublead["2017"]["eta1"].insert({"r93", { }});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt0", 0.842708000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt0", 0.017774600});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt1", 0.945462000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt1", 0.023360200});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt2", 0.962701000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt2", 0.018884600});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt3", 0.965889000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt3", 0.013939300});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt4", 0.969626000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt4", 0.013400600});
          idsf_sublead["2017"]["eta1"]["r93"].insert({"pt5", 0.977952000});
          idsfunc_sublead["2017"]["eta1"]["r93"].insert({"pt5", 0.011491900});
      idsf_sublead["2017"].insert({"eta2", { }});
      idsfunc_sublead["2017"].insert({"eta2", { }});
        idsf_sublead["2017"]["eta2"].insert({"r90", { }});
        idsfunc_sublead["2017"]["eta2"].insert({"r90", { }});
          idsf_sublead["2017"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2017"]["eta2"]["r90"].insert({"pt0", 0.0});


    idsf_lead.insert({"2018", { }});
    idsfunc_lead.insert({"2018", { }});
      idsf_lead["2018"].insert({"eta0", { }});
      idsfunc_lead["2018"].insert({"eta0", { }});
        idsf_lead["2018"]["eta0"].insert({"r90", { }});
        idsfunc_lead["2018"]["eta0"].insert({"r90", { }});
          idsf_lead["2018"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2018"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2018"]["eta0"].insert({"r91", { }});
        idsfunc_lead["2018"]["eta0"].insert({"r91", { }});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt0", 0.780938869});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt0", 0.004864648});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt1", 0.919135282});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt1", 0.003581765});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt2", 0.930928443});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt2", 0.004336357});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt3", 0.928327221});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt3", 0.002356278});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt4", 0.931627238});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt4", 0.004075929});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt5", 0.944861906});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt5", 0.015164412});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt6", 0.947722049});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt6", 0.025845703});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt7", 0.945577779});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt7", 0.013170956});
          idsf_lead["2018"]["eta0"]["r91"].insert({"pt8", 0.928267404});
          idsfunc_lead["2018"]["eta0"]["r91"].insert({"pt8", 0.035055827});
        idsf_lead["2018"]["eta0"].insert({"r92", { }});
        idsfunc_lead["2018"]["eta0"].insert({"r92", { }});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt0", 0.820914293});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt0", 0.001011564});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt1", 0.947347227});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt1", 0.001568624});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt2", 0.954592064});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt2", 0.003151450});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt3", 0.956382009});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt3", 0.001001644});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt4", 0.959609780});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt4", 0.001166888});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt5", 0.961361418});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt5", 0.002349434});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt6", 0.968326999});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt6", 0.001620777});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt7", 0.975550504});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt7", 0.002212991});
          idsf_lead["2018"]["eta0"]["r92"].insert({"pt8", 0.972812037});
          idsfunc_lead["2018"]["eta0"]["r92"].insert({"pt8", 0.006461995});
        idsf_lead["2018"]["eta0"].insert({"r93", { }});
        idsfunc_lead["2018"]["eta0"].insert({"r93", { }});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt0", 0.861269899});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt0", 0.003239360});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt1", 0.959156600});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt1", 0.001134487});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt2", 0.966150439});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt2", 0.001121107});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt3", 0.971342941});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt3", 0.001000526});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt4", 0.976606940});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt4", 0.001059896});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt5", 0.978579373});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt5", 0.001928226});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt6", 0.978487233});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt6", 0.003315547});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt7", 0.980704339});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt7", 0.001127807});
          idsf_lead["2018"]["eta0"]["r93"].insert({"pt8", 0.987216810});
          idsfunc_lead["2018"]["eta0"]["r93"].insert({"pt8", 0.001017551});
      idsf_lead["2018"].insert({"eta1", { }});
      idsfunc_lead["2018"].insert({"eta1", { }});
        idsf_lead["2018"]["eta1"].insert({"r90", { }});
        idsfunc_lead["2018"]["eta1"].insert({"r90", { }});
          idsf_lead["2018"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2018"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_lead["2018"]["eta1"].insert({"r91", { }});
        idsfunc_lead["2018"]["eta1"].insert({"r91", { }});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt0", 0.833663217});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt0", 0.005815814});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt1", 0.959330511});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt1", 0.001682000});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt2", 0.963932447});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt2", 0.001762105});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt3", 0.967334793});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt3", 0.001016837});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt4", 0.973084597});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt4", 0.001212774});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt5", 0.970071362});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt5", 0.001893430});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt6", 0.975354118});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt6", 0.002494532});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt7", 0.977838111});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt7", 0.002707403});
          idsf_lead["2018"]["eta1"]["r91"].insert({"pt8", 0.984013392});
          idsfunc_lead["2018"]["eta1"]["r91"].insert({"pt8", 0.003265289});
        idsf_lead["2018"]["eta1"].insert({"r92", { }});
        idsfunc_lead["2018"]["eta1"].insert({"r92", { }});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt0", 0.845371172});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt0", 0.002063441});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt1", 0.971218700});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt1", 0.001102060});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt2", 0.976836001});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt2", 0.001190716});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt3", 0.980361970});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt3", 0.001222205});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt4", 0.983873247});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt4", 0.004337875});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt5", 0.984967520});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt5", 0.009859004});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt6", 0.990564194});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt6", 0.002968815});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt7", 0.989140823});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt7", 0.003954124});
          idsf_lead["2018"]["eta1"]["r92"].insert({"pt8", 0.993155891});
          idsfunc_lead["2018"]["eta1"]["r92"].insert({"pt8", 0.002126684});
      idsf_lead["2018"].insert({"eta2", { }});
      idsfunc_lead["2018"].insert({"eta2", { }});
        idsf_lead["2018"]["eta2"].insert({"r90", { }});
        idsfunc_lead["2018"]["eta2"].insert({"r90", { }});
          idsf_lead["2018"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_lead["2018"]["eta2"]["r90"].insert({"pt0", 0.0});

    idsf_sublead.insert({"2018", { }});
    idsfunc_sublead.insert({"2018", { }});
      idsf_sublead["2018"].insert({"eta0", { }});
      idsfunc_sublead["2018"].insert({"eta0", { }});
        idsf_sublead["2018"]["eta0"].insert({"r90", { }});
        idsfunc_sublead["2018"]["eta0"].insert({"r90", { }});
          idsf_sublead["2018"]["eta0"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2018"]["eta0"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2018"]["eta0"].insert({"r91", { }});
        idsfunc_sublead["2018"]["eta0"].insert({"r91", { }});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt0", 0.913884494});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt0", 0.009517289});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt1", 0.939115197});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt1", 0.009146070});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt2", 0.947174456});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt2", 0.009407467});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt3", 0.965752087});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt3", 0.003902668});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt4", 0.972694647});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt4", 0.004627728});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt5", 0.973594440});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt5", 0.001845442});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt6", 0.973623474});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt6", 0.001954076});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt7", 0.978601841});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt7", 0.014369003});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt8", 0.970405929});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt8", 0.018451648});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt9", 0.966664489});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt9", 0.017007281});
          idsf_sublead["2018"]["eta0"]["r91"].insert({"pt10", 0.947327704});
          idsfunc_sublead["2018"]["eta0"]["r91"].insert({"pt10", 0.047975961});
        idsf_sublead["2018"]["eta0"].insert({"r92", { }});
        idsfunc_sublead["2018"]["eta0"].insert({"r92", { }});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt0", 0.943798395});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt0", 0.003316689});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt1", 0.975876626});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt1", 0.003744328});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt2", 0.982889507});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt2", 0.001468185});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt3", 0.987609719});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt3", 0.001825274});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt4", 0.989295959});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt4", 0.001088810});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt5", 0.989700703});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt5", 0.001093196});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt6", 0.988282541});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt6", 0.001045119});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt7", 0.987912109});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt7", 0.001058935});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt8", 0.991951085});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt8", 0.002761054});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt9", 0.991374293});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt9", 0.001572445});
          idsf_sublead["2018"]["eta0"]["r92"].insert({"pt10", 0.985228217});
          idsfunc_sublead["2018"]["eta0"]["r92"].insert({"pt10", 0.004677243});
        idsf_sublead["2018"]["eta0"].insert({"r93", { }});
        idsfunc_sublead["2018"]["eta0"].insert({"r93", { }});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt0", 0.967229704});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt0", 0.001583783});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt1", 0.980568240});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt1", 0.003240185});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt2", 0.985786488});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt2", 0.003539124});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt3", 0.989965533});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt3", 0.001430515});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt4", 0.992684722});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt4", 0.001003685});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt5", 0.993636346});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt5", 0.001001831});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt6", 0.993855036});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt6", 0.002266144});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt7", 0.994253518});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt7", 0.001086973});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt8", 0.994926328});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt8", 0.001139417});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt9", 0.994882077});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt9", 0.001725602});
          idsf_sublead["2018"]["eta0"]["r93"].insert({"pt10", 0.994925740});
          idsfunc_sublead["2018"]["eta0"]["r93"].insert({"pt10", 0.001273526});
      idsf_sublead["2018"].insert({"eta1", { }});
      idsfunc_sublead["2018"].insert({"eta1", { }});
        idsf_sublead["2018"]["eta1"].insert({"r90", { }});
        idsfunc_sublead["2018"]["eta1"].insert({"r90", { }});
          idsf_sublead["2018"]["eta1"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2018"]["eta1"]["r90"].insert({"pt0", 0.0});
        idsf_sublead["2018"]["eta1"].insert({"r91", { }});
        idsfunc_sublead["2018"]["eta1"].insert({"r91", { }});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt0", 0.912117162});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt0", 0.002233670});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt1", 0.954305588});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt1", 0.006573312});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt2", 0.958487902});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt2", 0.003897311});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt3", 0.966529284});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt3", 0.001350754});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt4", 0.972487735});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt4", 0.001195978});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt5", 0.975625653});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt5", 0.001021597});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt6", 0.979438512});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt6", 0.001099126});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt7", 0.978590269});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt7", 0.001301435});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt8", 0.977142848});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt8", 0.002335084});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt9", 0.978083049});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt9", 0.002683357});
          idsf_sublead["2018"]["eta1"]["r91"].insert({"pt10", 0.981152020});
          idsfunc_sublead["2018"]["eta1"]["r91"].insert({"pt10", 0.002623317});
        idsf_sublead["2018"]["eta1"].insert({"r92", { }});
        idsfunc_sublead["2018"]["eta1"].insert({"r92", { }});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt0", 0.943331900});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt0", 0.003617830});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt1", 0.983648649});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt1", 0.002260656});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt2", 0.984981002});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt2", 0.003790529});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt3", 0.989851679});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt3", 0.001246198});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt4", 0.992568512});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt4", 0.001056425});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt5", 0.993433762});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt5", 0.001000833});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt6", 0.993232320});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt6", 0.001056853});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt7", 0.993395672});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt7", 0.001039927});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt8", 0.995547720});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt8", 0.001048809});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt9", 0.993596642});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt9", 0.001035713});
          idsf_sublead["2018"]["eta1"]["r92"].insert({"pt10", 0.996633194});
          idsfunc_sublead["2018"]["eta1"]["r92"].insert({"pt10", 0.002228740});
      idsf_sublead["2018"].insert({"eta2", { }});
      idsfunc_sublead["2018"].insert({"eta2", { }});
        idsf_sublead["2018"]["eta2"].insert({"r90", { }});
        idsfunc_sublead["2018"]["eta2"].insert({"r90", { }});
          idsf_sublead["2018"]["eta2"]["r90"].insert({"pt0", 1.0});
          idsfunc_sublead["2018"]["eta2"]["r90"].insert({"pt0", 0.0});
  }

  inline float get_lowMassHggTriggerSF(const float pt, const float eta, const float r9, const TString photon, const bool EBHiR9, const TString year, const TString variation ) {
    if ( year!="2016nonAPV" && year!="2016APV" && year!="2017" && year!="2018" ) {
      std::cout << "WARNING: unknown year, returning unity low mass trigger SF!" << std::endl;
      return 1.0;
    }
    TString actualYear;
    if (year=="2016nonAPV" || year=="2016APV") actualYear = "2016";
    else actualYear = year;

    if ( photon!="Lead" && photon!="Sublead" ) {
      std::cout << "WARNING: unknown photon, returning unity low mass trigger SF!" << std::endl;
      return 1.0;
    }

    float abseta = fabs(eta);
    set_etabins();
    set_r9bins(abseta, photon, EBHiR9, actualYear);
    set_ptbins(abseta, r9, photon, EBHiR9, actualYear);
    TString etabin = get_etaBin(abseta);
    TString r9bin = get_r9Bin(r9);
    TString ptbin = get_ptBin(pt);

    if ( actualYear == "2016" && EBHiR9 ) actualYear = "2016EBHiR9";
    if ( variation == "central" ) return (photon=="Lead" ? idsf_lead[actualYear][etabin][r9bin][ptbin] : idsf_sublead[actualYear][etabin][r9bin][ptbin]);
    else if ( variation == "up" ) return (photon=="Lead" ? idsf_lead[actualYear][etabin][r9bin][ptbin]+idsfunc_lead[actualYear][etabin][r9bin][ptbin] : idsf_sublead[actualYear][etabin][r9bin][ptbin]+idsfunc_sublead[actualYear][etabin][r9bin][ptbin]);
    else if ( variation == "down" ) return (photon=="Lead" ? idsf_lead[actualYear][etabin][r9bin][ptbin]-idsfunc_lead[actualYear][etabin][r9bin][ptbin] : idsf_sublead[actualYear][etabin][r9bin][ptbin]-idsfunc_sublead[actualYear][etabin][r9bin][ptbin]);
    else {
      std::cout << "WARNING: unknown variation (central, up, down), returning unity scale factor!" << std::endl;
      return 1.0;
    }
  }
}

#endif

