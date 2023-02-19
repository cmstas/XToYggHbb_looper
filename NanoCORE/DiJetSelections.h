#ifndef DiJetSELECTIONS_H
#define DiJetSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"
#include "DiPhotonSelections.h"

struct Jet {
    Jet(unsigned int idx = 0) : idx_(idx) {
        try
        {
            pt_ = nt.Jet_pt_nom()[idx_];
        }
        catch(const std::exception& e)
        {
            pt_ = nt.Jet_pt()[idx_];
        }
        
        mass_ = nt.Jet_mass()[idx_];
        eta_ = nt.Jet_eta()[idx_];
        phi_ = nt.Jet_phi()[idx_];
        jetId_ = nt.Jet_jetId()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.Jet_eta()[idx_], nt.Jet_phi()[idx_], nt.Jet_mass()[idx_]);
        btagDeepFlavB_ = nt.Jet_btagDeepFlavB()[idx_];
        try {btagSF_deepjet_shape_ = nt.Jet_btagSF_deepjet_shape()[idx_];} catch(const std::exception& e) {btagSF_deepjet_shape_=1;}
    }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    int id() { return id_; }
    unsigned int idx() { return idx_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int jetId() { return jetId_; }
    float btagDeepFlavB() { return btagDeepFlavB_; }
    float btagSF_deepjet_shape() { return btagSF_deepjet_shape_; }

  private:
    int id_;
    float pt_ = 0.;
    float eta_ = 0.;
    float mass_ = 0;
    TLorentzVector p4_;
    float phi_ = 0.;
    int jetId_ = 0;
    unsigned int idx_;
    float btagDeepFlavB_ = 0;
    float btagSF_deepjet_shape_ = 1;
};

vector<Jet> getJets(Photons photons);
typedef std::vector<Jet> Jets;

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

inline bool sortBybscore(Jet &p1, Jet &p2)
{
        return p1.btagDeepFlavB() > p2.btagDeepFlavB();    
}

DiJets DiJetPreselection(Jets &jets); 

#endif
