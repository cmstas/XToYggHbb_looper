#include "DiJetSelections.h"

using namespace tas;

Jets getJets(Photons photons, const int JESUnc, const int JERUnc, const int HEMCheck) {
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
        // Check HEM issue: https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
        if ( HEMCheck!=0 ) {
          if ( cand_jet.phi() > -1.57 && cand_jet.phi() < -0.87 ) {
            if ( cand_jet.eta() >= -2.5 && cand_jet.eta() < -1.3 ) cand_jet.setPt( cand_jet.pt() * 0.80 );
            if ( cand_jet.eta() >= -3.0 && cand_jet.eta() < -2.5 ) cand_jet.setPt( cand_jet.pt() * 0.65 );
          }
        }

        if ( !(abs(cand_jet.eta()) < 2.4) ) continue;
        if ( !(cand_jet.pt() > 25) ) continue;
        if ( !(cand_jet.jetId()>1) ) continue; // Apply the Jet ID (at least tight WP), as per recommendations (https://twiki.cern.ch/twiki/bin/view/CMS/JetID)
        if ( cand_jet.pt() < 50 && !(cand_jet.puId()>5) ) continue; // Apply the Jet PU ID (at least medium WP) to jets below 50 GeV, as per recommendations (https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID)
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
