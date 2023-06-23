#include "DiJetSelections.h"

using namespace tas;

Jets getJets(Photons photons, const int JESUnc, const int JERUnc) {
    Jets jets;
    for (unsigned int ijet = 0; ijet < nt.nJet(); ijet++) {
        Jet cand_jet = Jet(ijet);

        // Apply JES unc
        if ( JESUnc!=0 ) {
          // JESUnc is only an uncertainty => JESUnc==1 does not change the nominal value
          if ( JESUnc==2  ) cand_jet.setPt( cand_jet.pt_jesTotalUp() );
          if ( JESUnc==-2 ) cand_jet.setPt( cand_jet.pt_jesTotalDown() );
        }
        // Apply JER unc
        if ( JERUnc!=0 ) {
          // JERUnc is only an uncertainty => JERUnc==1 does not change the nominal value
          if ( JERUnc==2  ) cand_jet.setPt( cand_jet.pt_jerUp() );
          if ( JERUnc==-2 ) cand_jet.setPt( cand_jet.pt_jerDown() );
        }

        if ( !(abs(cand_jet.eta()) < 2.4) ) continue;
        if ( !(cand_jet.pt() > 25) ) continue;
        if ( !(cand_jet.jetId()>=1) ) continue;
        bool clean_with_photon = true;
        for (unsigned int iphoton = 0; iphoton < photons.size(); iphoton++){
            if (cand_jet.p4().DeltaR(photons.at(iphoton).p4())<0.4) clean_with_photon = false;
        }
        if (clean_with_photon == false) continue;
        jets.push_back(cand_jet);
    }
    sort(jets.begin(), jets.end(), sortBybscore);

    return jets;
}

DiJets DiJetPreselection(Jets &jets) {
    DiJets dijets;   
    if (jets.size() > 1) {
        DiJet dijet = DiJet(jets[0], jets[1]);
        dijets.push_back(dijet);
    }
    return dijets;    
}

GenJets getGenJets() {
    GenJets Genjets;
    for (unsigned int ijet = 0; ijet < nt.nGenJet(); ijet++) {
        GenJet cand_genjet = GenJet(ijet);
        if ( !(abs(cand_genjet.eta()) < 2.4) ) continue;
        if ( !(cand_genjet.pt() > 25) ) continue;
        Genjets.push_back(cand_genjet);
    }
    return Genjets;
}

GenDiJets GenDiJetPreselection(GenJets &Genjets) {
    GenDiJets Gendijets;   
    if (Genjets.size() > 1) {
        GenDiJet Gendijet = GenDiJet(Genjets[0], Genjets[1]);
        Gendijets.push_back(Gendijet);
    }
    return Gendijets;    
}

FatJets getFatJets(Photons photons) {
    FatJets Fatjets;
    FatJet lead_FatJet;
    float lead_Hbb_score=0;
    if (nt.nFatJet()==0) return Fatjets;
    for (unsigned int ijet = 0; ijet < nt.nFatJet(); ijet++) {
        FatJet cand_fatjet = FatJet(ijet);
        if ( !(abs(cand_fatjet.eta()) < 2.5) ) continue;
        if ( !(cand_fatjet.pt() > 200) ) continue;
        if ( !(cand_fatjet.mass() > 40) ) continue;
        if ( !(cand_fatjet.jetId() < 1) ) continue;
        bool clean_with_photon = true;
        for (unsigned int iphoton = 0; iphoton < photons.size(); iphoton++){
            if (cand_fatjet.p4().DeltaR(photons.at(iphoton).p4())<0.8) clean_with_photon = false;
        }

        if (clean_with_photon == false) continue;
        if (cand_fatjet.Hbb_score()>lead_Hbb_score) {
            lead_Hbb_score = cand_fatjet.Hbb_score();
            lead_FatJet = cand_fatjet;
        }
    }
    //Fatjets are naturally sorted by pt in NanoAOD
    if (lead_Hbb_score>0) Fatjets.push_back(lead_FatJet);
    return Fatjets;
}