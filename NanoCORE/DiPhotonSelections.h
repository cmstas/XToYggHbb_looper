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
        //p4_ = nt.Photon_p4()[idx_]; need to convert to TLorentzVector to compute DeltaR later
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
        trkSumPtHollowConeDR03_ = nt.Photon_trkSumPtHollowConeDR03()[idx_];
        //idlevel_ = whichPhotonLevel(id_, idx_);
        try { pt_ScaleUp_ = nt.Photon_pt_ScaleUp()[idx_]; }
        catch(const std::exception& e) { pt_ScaleUp_ = 1.0; }
        try { pt_ScaleDown_ = nt.Photon_pt_ScaleDown()[idx_]; }
        catch(const std::exception& e) { pt_ScaleDown_ = 1.0; }
        dEsigmaUp_ = nt.Photon_dEsigmaUp()[idx_];
        dEsigmaDown_ = nt.Photon_dEsigmaDown()[idx_];
        fixedGridRhoAll_ = nt.Rho_fixedGridRhoAll();
        try { genPartFlav_ = nt.Photon_genPartFlav()[idx_]; } catch(const std::exception& e) { genPartFlav_ = 0; }
    }
    void setPt(float pt) { pt_ = pt; }
    void setMvaID(float mvaID) { mvaID_ = mvaID; }
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
    float perEvtRho() { return fixedGridRhoAll_; }
    float isScEtaEE() { return isScEtaEE_; }
    float isScEtaEB() { return isScEtaEB_; }
    float trkSumPtHollowConeDR03() { return trkSumPtHollowConeDR03_; }
    float pt_ScaleUp() { return pt_ScaleUp_; }
    float pt_ScaleDown() { return pt_ScaleDown_; }
    float dEsigmaUp() { return dEsigmaUp_; }
    float dEsigmaDown() { return dEsigmaDown_; }
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
    float fixedGridRhoAll_ = 0.;
    float isScEtaEB_ = 0;
    float isScEtaEE_ = 0;
    float trkSumPtHollowConeDR03_ = 0;
    float pt_ScaleUp_ = 0.;
    float pt_ScaleDown_ = 0.;
    float dEsigmaUp_ = 0.;
    float dEsigmaDown_ = 0.;
    //int idlevel_ = SS::IDdefault;
    unsigned char genPartFlav_ = 0;
};

vector<Photon> getPhotons(const TString year, const int lowMassMode, const int fnufUnc, const int materialUnc, const int PhoScaleUnc, const int PhoSmearUnc);
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

DiPhotons DiPhotonPreselection(Photons &photons, const int lowMassMode); 
bool UseLowR9Photon(Photon pho, bool isEB);

#endif
