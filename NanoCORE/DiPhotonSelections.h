#ifndef DiPhotonSELECTIONS_H
#define DiPhotonSELECTIONS_H
#include "Nano.h"
#include "Base.h"
#include "TLorentzVector.h"

struct Photon {
    Photon(unsigned int idx = 0) : idx_(idx) {
        pt_ = nt.Photon_pt()[idx_];
        eta_ = nt.Photon_eta()[idx_];
        phi_ = nt.Photon_phi()[idx_];
        //mass_ = nt.Photon_mass()[idx_]; // Not in newest custom NanoAOD
//        p4_ = nt.Photon_p4()[idx_]; need to convert to TLorentzVector to compute DeltaR later
        p4_.SetPtEtaPhiM(pt_, eta_, phi_, 0 );
        //id_ = nt.Photon_pdgId()[idx_]; // Not in newest custom NanoAOD
        r9_ = nt.Photon_r9()[idx_];
        chargedHadIso_ = nt.Photon_pfChargedIsoPFPV()[idx_]; // Equivalent to chargedHadronIso, as noted here: https://github.com/cms-sw/cmssw/pull/36526/#discussion_r772295841
        hoe_ = nt.Photon_hoe()[idx_];
        phoIso_ = nt.Photon_pfPhoIso03()[idx_];
        sieie_ = nt.Photon_sieie()[idx_];
        eveto_ = nt.Photon_electronVeto()[idx_];
        pixelSeed_ = nt.Photon_pixelSeed()[idx_];
        mvaID_ = nt.Photon_mvaID()[idx_];
        isScEtaEB_ = nt.Photon_isScEtaEB()[idx_];
        isScEtaEE_ = nt.Photon_isScEtaEE()[idx_];
        //trkSumPtHollowConeDR03_ = nt.Photon_trkSumPtHollowConeDR03()[idx_]; // To be readded
        //idlevel_ = whichPhotonLevel(id_, idx_);
        try { fixedGridRhoFastjetAll_ = nt.Rho_fixedGridRhoFastjetAll(); } catch(const std::exception& e) { fixedGridRhoFastjetAll_ = nt.fixedGridRhoFastjetAll(); }
        try { genPartFlav_ = nt.Photon_genPartFlav()[idx_]; } catch(const std::exception& e) { genPartFlav_ = -1; }
    }
//    void setGenPartFlav(unsigned int idx) { genPartFlav_ = nt.Photon_genPartFlav()[idx_]; }
    //void set_idlevel(int idlevel) { idlevel_ = idlevel; }
    //int id() { return id_; } // Not in newest custom NanoAOD
    unsigned int idx() { return idx_; }
    //int idlevel() { return idlevel_; }
    TLorentzVector p4() { return p4_; }
    float pt() { return pt_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    //float mass() { return mass_; } // Not in newest custom NanoAOD
    float r9() { return r9_; }
    float chargedHadIso() { return chargedHadIso_; }
    float hoe() { return hoe_; }
    float phoIso() { return phoIso_; }
    float trkIso() { return trkIso_; }
    float sieie() { return sieie_; }
    bool eveto() { return eveto_; }
    bool pixelSeed() { return pixelSeed_;}
    float mvaID() { return mvaID_; }
    float perEvtRho() { return fixedGridRhoFastjetAll_; }
    float isScEtaEE() { return isScEtaEE_; }
    float isScEtaEB() { return isScEtaEB_; }
    //float trkSumPtHollowConeDR03() { return trkSumPtHollowConeDR03_; } // To be readded
    unsigned char genPartFlav() { return genPartFlav_; }

  private:
    //int id_; // Not in newest custom NanoAOD
    float pt_ = 0.;
    float eta_ = 0.;
    float phi_ = 0.;
    //float mass_ = 0.; // Not in newest custom NanoAOD
    TLorentzVector p4_;
    unsigned int idx_;
    float r9_ = 0.;
    float chargedHadIso_ = 0.;
    float hoe_ = 0.;
    float phoIso_ = 0.;
    float trkIso_ = 0.;
    float sieie_ = 0.;
    bool eveto_ = 0.;
    bool pixelSeed_ = 0.;
    float mvaID_ = 0.;
    float fixedGridRhoFastjetAll_ = 0.; // this variable is the same for each event
    float isScEtaEB_ = 0;
    float isScEtaEE_ = 0;
    //float trkSumPtHollowConeDR03_ = 0; // To be readded
    //int idlevel_ = SS::IDdefault;
    unsigned char genPartFlav_ = -1;
};

vector<Photon> getPhotons();
typedef std::vector<Photon> Photons;

struct DiPhoton{
    Photon leadPho;
    Photon subleadPho;
    TLorentzVector p4;
    float dR;
    DiPhoton(Photon p1, Photon p2)
    {
        leadPho = p1;
        subleadPho = p2;
        TLorentzVector leadPhop4 = TLorentzVector();
        TLorentzVector subleadPhop4 = TLorentzVector();
        leadPhop4.SetPtEtaPhiM(p1.pt(), p1.eta(), p1.phi(), 0);
        subleadPhop4.SetPtEtaPhiM(p2.pt(), p2.eta(), p2.phi(), 0);
        p4 = leadPhop4 + subleadPhop4;
        dR = leadPhop4.DeltaR(subleadPhop4);
    }
};

//typedef std::pair<Photon, Photon> DiPhoton;
typedef std::vector<DiPhoton> DiPhotons;

inline bool sortByPt(Photon &p1, Photon &p2)
{
        return p1.pt() > p2.pt();    
}

DiPhotons DiPhotonPreselection(Photons &photons); 
bool UseLowR9Photon(Photon pho, bool isEB);
/*
SS::IDLevel whichLeptonLevel(int id, int idx);

typedef std::pair<Lepton, Lepton> Hyp;
typedef std::vector<Lepton> Leptons;

std::ostream &operator<<(std::ostream &os, Lepton &lep) {
    std::string lepstr = (abs(lep.id()) == 11) ? "Electron" : "Muon";
    return os << "<" << lepstr << " id=" << std::showpos << setw(3) << lep.id() << std::noshowpos << ", idx=" << setw(2)
              << lep.idx() << ", level=" << lep.idlevel() << ", (pT,eta)="
              << "(" << lep.pt() << "," << lep.eta() << ")>";
}
template <typename T1, typename T2> std::ostream &operator<<(std::ostream &os, std::pair<T1, T2> &p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

vector<Lepton> getLeptons();
std::tuple<int, int, float> getJetInfo(vector<Lepton> &leps, int variation = 0);
std::pair<int, int> makesResonance(Leptons &leps, Lepton lep1, Lepton lep2, float mass, float window);
std::pair<int, Hyp> getBestHyp(vector<Lepton> &leptons, bool verbose);
bool isLeptonLevel(SS::IDLevel idlevel, int id, int idx);
void dumpLeptonProperties(Lepton lep);
*/
#endif
