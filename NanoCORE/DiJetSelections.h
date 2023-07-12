#ifndef DiJetSELECTIONS_H
#define DiJetSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"
#include "DiPhotonSelections.h"

struct Jet {
    Jet(unsigned int idx = 0) : idx_(idx) {
        try { pt_ = nt.Jet_pt_nom()[idx_]; }
        catch(const std::exception& e) { pt_ = nt.Jet_pt()[idx_]; }
        try { pt_jerUp_ = nt.Jet_pt_jerUp()[idx_]; }
        catch(const std::exception& e) { pt_jerUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jerDown_ = nt.Jet_pt_jerDown()[idx_]; }
        catch(const std::exception& e) { pt_jerDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesTotalUp_ = nt.Jet_pt_jesTotalUp()[idx_]; }
        catch(const std::exception& e) { pt_jesTotalUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesTotalDown_ = nt.Jet_pt_jesTotalDown()[idx_]; }
        catch(const std::exception& e) { pt_jesTotalDown_ = nt.Jet_pt()[idx_]; }
        mass_ = nt.Jet_mass()[idx_];
        eta_ = nt.Jet_eta()[idx_];
        phi_ = nt.Jet_phi()[idx_];
        jetId_ = nt.Jet_jetId()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.Jet_eta()[idx_], nt.Jet_phi()[idx_], nt.Jet_mass()[idx_]);
        btagDeepFlavB_ = nt.Jet_btagDeepFlavB()[idx_];
        try {btagSF_deepjet_shape_ = nt.Jet_btagSF_deepjet_shape()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_=1;}
        try {btagSF_deepjet_shape_up_hf_ = nt.Jet_btagSF_deepjet_shape_up_hf()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_hf_=1;}
        try {btagSF_deepjet_shape_down_hf_ = nt.Jet_btagSF_deepjet_shape_down_hf()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_hf_=1;}
        try {btagSF_deepjet_shape_up_hfstats1_ = nt.Jet_btagSF_deepjet_shape_up_hfstats1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_hfstats1_=1;}
        try {btagSF_deepjet_shape_down_hfstats1_ = nt.Jet_btagSF_deepjet_shape_down_hfstats1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_hfstats1_=1;}
        try {btagSF_deepjet_shape_up_hfstats2_ = nt.Jet_btagSF_deepjet_shape_up_hfstats2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_hfstats2_=1;}
        try {btagSF_deepjet_shape_down_hfstats2_ = nt.Jet_btagSF_deepjet_shape_down_hfstats2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_hfstats2_=1;}
        try {btagSF_deepjet_shape_up_lf_ = nt.Jet_btagSF_deepjet_shape_up_lf()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_lf_=1;}
        try {btagSF_deepjet_shape_down_lf_ = nt.Jet_btagSF_deepjet_shape_down_lf()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_lf_=1;}
        try {btagSF_deepjet_shape_up_lfstats1_ = nt.Jet_btagSF_deepjet_shape_up_lfstats1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_lfstats1_=1;}
        try {btagSF_deepjet_shape_down_lfstats1_ = nt.Jet_btagSF_deepjet_shape_down_lfstats1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_lfstats1_=1;}
        try {btagSF_deepjet_shape_up_lfstats2_ = nt.Jet_btagSF_deepjet_shape_up_lfstats2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_lfstats2_=1;}
        try {btagSF_deepjet_shape_down_lfstats2_ = nt.Jet_btagSF_deepjet_shape_down_lfstats2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_lfstats2_=1;}
        try {btagSF_deepjet_shape_up_jes_ = nt.Jet_btagSF_deepjet_shape_up_jes()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_jes_=1;}
        try {btagSF_deepjet_shape_down_jes_ = nt.Jet_btagSF_deepjet_shape_down_jes()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_jes_=1;}
        try {btagSF_deepjet_shape_up_cferr1_ = nt.Jet_btagSF_deepjet_shape_up_cferr1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_cferr1_=1;}
        try {btagSF_deepjet_shape_down_cferr1_ = nt.Jet_btagSF_deepjet_shape_down_cferr1()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_cferr1_=1;}
        try {btagSF_deepjet_shape_up_cferr2_ = nt.Jet_btagSF_deepjet_shape_up_cferr2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_up_cferr2_=1;}
        try {btagSF_deepjet_shape_down_cferr2_ = nt.Jet_btagSF_deepjet_shape_down_cferr2()[idx_];}
        catch(const std::exception& e) {btagSF_deepjet_shape_down_cferr2_=1;}
    }
    void setPt(float pt) { pt_ = pt; }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float pt_jerUp() { return pt_jerUp_; }
    float pt_jerDown() { return pt_jerDown_; }
    float pt_jesTotalUp() { return pt_jesTotalUp_; }
    float pt_jesTotalDown() { return pt_jesTotalDown_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int jetId() { return jetId_; }
    float btagDeepFlavB() { return btagDeepFlavB_; }
    float btagSF_deepjet_shape() { return btagSF_deepjet_shape_; }
    float btagSF_deepjet_shape_up_hf() { return btagSF_deepjet_shape_up_hf_; }
    float btagSF_deepjet_shape_down_hf() { return btagSF_deepjet_shape_down_hf_; }
    float btagSF_deepjet_shape_up_hfstats1() { return btagSF_deepjet_shape_up_hfstats1_; }
    float btagSF_deepjet_shape_down_hfstats1() { return btagSF_deepjet_shape_down_hfstats1_; }
    float btagSF_deepjet_shape_up_hfstats2() { return btagSF_deepjet_shape_up_hfstats2_; }
    float btagSF_deepjet_shape_down_hfstats2() { return btagSF_deepjet_shape_down_hfstats2_; }
    float btagSF_deepjet_shape_up_lf() { return btagSF_deepjet_shape_up_lf_; }
    float btagSF_deepjet_shape_down_lf() { return btagSF_deepjet_shape_down_lf_; }
    float btagSF_deepjet_shape_up_lfstats1() { return btagSF_deepjet_shape_up_lfstats1_; }
    float btagSF_deepjet_shape_down_lfstats1() { return btagSF_deepjet_shape_down_lfstats1_; }
    float btagSF_deepjet_shape_up_lfstats2() { return btagSF_deepjet_shape_up_lfstats2_; }
    float btagSF_deepjet_shape_down_lfstats2() { return btagSF_deepjet_shape_down_lfstats2_; }
    float btagSF_deepjet_shape_up_jes() { return btagSF_deepjet_shape_up_jes_; }
    float btagSF_deepjet_shape_down_jes() { return btagSF_deepjet_shape_down_jes_; }
    float btagSF_deepjet_shape_up_cferr1() { return btagSF_deepjet_shape_up_cferr1_; }
    float btagSF_deepjet_shape_down_cferr1() { return btagSF_deepjet_shape_down_cferr1_; }
    float btagSF_deepjet_shape_up_cferr2() { return btagSF_deepjet_shape_up_cferr2_; }
    float btagSF_deepjet_shape_down_cferr2() { return btagSF_deepjet_shape_down_cferr2_; }

  private:
    int id_;
    float pt_ = 0.;
    float pt_jerUp_ = 0.;
    float pt_jerDown_ = 0.;
    float pt_jesTotalUp_ = 0.;
    float pt_jesTotalDown_ = 0.;
    float eta_ = 0.;
    float mass_ = 0;
    TLorentzVector p4_;
    float phi_ = 0.;
    int jetId_ = 0;
    unsigned int idx_;
    float btagDeepFlavB_ = 0.;
    float btagSF_deepjet_shape_ = 1.0;
    float btagSF_deepjet_shape_up_hf_ = 1.0;
    float btagSF_deepjet_shape_down_hf_ = 1.0;
    float btagSF_deepjet_shape_up_hfstats1_ = 1.0;
    float btagSF_deepjet_shape_down_hfstats1_ = 1.0;
    float btagSF_deepjet_shape_up_hfstats2_ = 1.0;
    float btagSF_deepjet_shape_down_hfstats2_ = 1.0;
    float btagSF_deepjet_shape_up_lf_ = 1.0;
    float btagSF_deepjet_shape_down_lf_ = 1.0;
    float btagSF_deepjet_shape_up_lfstats1_ = 1.0;
    float btagSF_deepjet_shape_down_lfstats1_ = 1.0;
    float btagSF_deepjet_shape_up_lfstats2_ = 1.0;
    float btagSF_deepjet_shape_down_lfstats2_ = 1.0;
    float btagSF_deepjet_shape_up_jes_ = 1.0;
    float btagSF_deepjet_shape_down_jes_ = 1.0;
    float btagSF_deepjet_shape_up_cferr1_ = 1.0;
    float btagSF_deepjet_shape_down_cferr1_ = 1.0;
    float btagSF_deepjet_shape_up_cferr2_ = 1.0;
    float btagSF_deepjet_shape_down_cferr2_ = 1.0;
};

struct GenJet {
    GenJet(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.GenJet_pt()[idx_];
        mass_ = nt.GenJet_mass()[idx_];
        eta_ = nt.GenJet_eta()[idx_];
        phi_ = nt.GenJet_phi()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.GenJet_eta()[idx_], nt.GenJet_phi()[idx_], nt.GenJet_mass()[idx_]);
    }
    TLorentzVector p4() { return p4_; }
    unsigned int idx() { return idx_; }
    float pt() { return pt_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }

  private:
    float pt_ = 0.;
    float eta_ = 0.;
    float mass_ = 0.;
    float phi_ = 0.;
    TLorentzVector p4_;
    unsigned int idx_;
};

struct FatJet {
    FatJet(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.FatJet_pt()[idx_];
        mass_ = nt.FatJet_msoftdrop()[idx_];
        eta_ = nt.FatJet_eta()[idx_];
        phi_ = nt.FatJet_phi()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.FatJet_eta()[idx_], nt.FatJet_phi()[idx_], nt.FatJet_msoftdrop()[idx_]);
        jetId_ = nt.Jet_jetId()[idx_];
        Hbb_score_ = nt.FatJet_particleNetMD_Xbb()[idx_]/(nt.FatJet_particleNetMD_Xbb()[idx_]+nt.FatJet_particleNetMD_QCD()[idx_]);
        subJetIdx1_ = nt.FatJet_subJetIdx1()[idx_];
        subJetIdx2_ = nt.FatJet_subJetIdx2()[idx_];
        nSubJet_ = nt.nSubJet();
        if (subJetIdx1_!=-1 && subJetIdx2_!=-1)
        {
            if (nt.SubJet_btagDeepB()[subJetIdx1_]>nt.SubJet_btagDeepB()[subJetIdx2_])
            {
                SubJet1_p4_.SetPtEtaPhiM(nt.SubJet_pt()[subJetIdx1_], nt.SubJet_eta()[subJetIdx1_], nt.SubJet_phi()[subJetIdx1_], nt.SubJet_mass()[subJetIdx1_]);
                SubJet2_p4_.SetPtEtaPhiM(nt.SubJet_pt()[subJetIdx2_], nt.SubJet_eta()[subJetIdx2_], nt.SubJet_phi()[subJetIdx2_], nt.SubJet_mass()[subJetIdx2_]);
                subjet1_bscore_ = nt.SubJet_btagDeepB()[subJetIdx1_];
                subjet2_bscore_ = nt.SubJet_btagDeepB()[subJetIdx2_];
            }
            if (nt.SubJet_btagDeepB()[subJetIdx2_]>nt.SubJet_btagDeepB()[subJetIdx1_])
            {
                SubJet2_p4_.SetPtEtaPhiM(nt.SubJet_pt()[subJetIdx1_], nt.SubJet_eta()[subJetIdx1_], nt.SubJet_phi()[subJetIdx1_], nt.SubJet_mass()[subJetIdx1_]);
                SubJet1_p4_.SetPtEtaPhiM(nt.SubJet_pt()[subJetIdx2_], nt.SubJet_eta()[subJetIdx2_], nt.SubJet_phi()[subJetIdx2_], nt.SubJet_mass()[subJetIdx2_]);
                subjet1_bscore_ = nt.SubJet_btagDeepB()[subJetIdx2_];
                subjet2_bscore_ = nt.SubJet_btagDeepB()[subJetIdx1_];
            }
        }
    }
    TLorentzVector p4() { return p4_; }
    TLorentzVector SubJet1_p4() { return SubJet1_p4_; }
    TLorentzVector SubJet2_p4() { return SubJet2_p4_; }
    unsigned int idx() { return idx_; }
    float pt() { return pt_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int jetId() { return jetId_; }
    float Hbb_score() {return Hbb_score_;}
    float subjet1_bscore() {return subjet1_bscore_;}
    float subjet2_bscore() {return subjet2_bscore_;}
    int subJetIdx1() {return subJetIdx1_;}
    int subJetIdx2() {return subJetIdx2_;}
    int nSubJet() {return nSubJet_;}

  private:
    float pt_ = 0.;
    float eta_ = 0.;
    float mass_ = 0.;
    float phi_ = 0.;
    TLorentzVector p4_;
    TLorentzVector SubJet1_p4_;
    TLorentzVector SubJet2_p4_;
    unsigned int idx_;
    int jetId_ = 0;
    float Hbb_score_ = 0;
    int subJetIdx1_ = 0;
    int subJetIdx2_ = 0;
    int nSubJet_ = 0;
    float subjet1_bscore_ = 0;
    float subjet2_bscore_ = 0;
};

vector<Jet> getJets(Photons photons, const int JESUnc, const int JERUnc);
typedef std::vector<Jet> Jets;
vector<GenJet> getGenJets();
typedef std::vector<GenJet> GenJets;
vector<FatJet> getFatJets(Photons photons);
typedef std::vector<FatJet> FatJets;

struct DiJet{
    Jet leadJet;
    Jet subleadJet;
    TLorentzVector p4;
    float dR;
    DiJet(Jet p1, Jet p2)
    {
        if (p1.pt() > p2.pt()) {
            leadJet = p1;
            subleadJet = p2;
        } else {
            leadJet = p2;
            subleadJet = p1;
        }
        TLorentzVector leadJetp4 = leadJet.p4();
        TLorentzVector subleadJetp4 = subleadJet.p4();
        p4 = leadJetp4 + subleadJetp4;
        dR = leadJetp4.DeltaR(subleadJetp4);
    }
};

typedef std::vector<DiJet> DiJets;

struct GenDiJet{
    GenJet leadGenJet;
    GenJet subleadGenJet;
    TLorentzVector p4;
    float dR;
    GenDiJet(GenJet p1, GenJet p2)
    {
        if (p1.pt() > p2.pt()) {
            leadGenJet = p1;
            subleadGenJet = p2;
        } else {
            leadGenJet = p2;
            subleadGenJet = p1;
        }
        TLorentzVector leadGenJetp4 = leadGenJet.p4();
        TLorentzVector subleadGenJetp4 = subleadGenJet.p4();
        p4 = leadGenJetp4 + subleadGenJetp4;
        dR = leadGenJetp4.DeltaR(subleadGenJetp4);
    }
};

typedef std::vector<GenDiJet> GenDiJets;

inline bool sortBybscore(Jet &p1, Jet &p2)
{
        return p1.btagDeepFlavB() > p2.btagDeepFlavB();    
}

DiJets DiJetPreselection(Jets &jets); 
GenDiJets GenDiJetPreselection(GenJets &Genjets); 

#endif
