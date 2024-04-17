#include "DiJetSelections.h"

using namespace tas;

Jets getJets(Photons photons, const int JESUnc, const int JERUnc, const int HEMCheck) {
    Jets jets;
    for (unsigned int ijet = 0; ijet < nt.nJet(); ijet++) {
        Jet cand_jet = Jet(ijet);

        // Apply JES unc
        if ( JESUnc!=0 ) {
          // JESUnc is only an uncertainty => JESUnc==1 does not change the nominal value
          // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources#Run_2_reduced_set_of_uncertainty
          // https://docs.google.com/spreadsheets/d/1Feuj1n0MdotcPq19Mht7SUIgvkXkA4hiB0BxEuBShLw/edit#gid=1345121349
          // Absolute up
          if ( JESUnc==2   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesAbsoluteMPFBiasUp() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesAbsoluteScaleUp() - cand_jet.pt(),    2) +
                                                                  pow(cand_jet.pt_jesFragmentationUp() - cand_jet.pt(),    2) +
                                                                  pow(cand_jet.pt_jesPileUpDataMCUp() - cand_jet.pt(),     2) +
                                                                  pow(cand_jet.pt_jesPileUpPtRefUp() - cand_jet.pt(),      2) +
                                                                  pow(cand_jet.pt_jesRelativeFSRUp() - cand_jet.pt(),      2) +
                                                                  pow(cand_jet.pt_jesSinglePionECALUp() - cand_jet.pt(),   2) +
                                                                  pow(cand_jet.pt_jesSinglePionHCALUp() - cand_jet.pt(),   2)));
          // Absolute down
          if ( JESUnc==-2  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesAbsoluteMPFBiasDown() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesAbsoluteScaleDown() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesFragmentationDown() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesPileUpDataMCDown() - cand_jet.pt(),   2) +
                                                                  pow(cand_jet.pt_jesPileUpPtRefDown() - cand_jet.pt(),    2) +
                                                                  pow(cand_jet.pt_jesRelativeFSRDown() - cand_jet.pt(),    2) +
                                                                  pow(cand_jet.pt_jesSinglePionECALDown() - cand_jet.pt(), 2) +
                                                                  pow(cand_jet.pt_jesSinglePionHCALDown() - cand_jet.pt(), 2)));

          // Absolute_year up
          if ( JESUnc==3   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesAbsoluteStatUp() - cand_jet.pt(),     2) +
                                                                  pow(cand_jet.pt_jesRelativeStatFSRUp() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesTimePtEtaUp() - cand_jet.pt(),        2)));
          // Absolute_year down
          if ( JESUnc==-3  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesAbsoluteStatDown() - cand_jet.pt(),   2) +
                                                                  pow(cand_jet.pt_jesRelativeStatFSRDown() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesTimePtEtaDown() - cand_jet.pt(),      2)));

          // FlavorQCD up
          if ( JESUnc==4   ) cand_jet.setPt( cand_jet.pt_jesFlavorQCDUp()  );
          // FlavorQCD down
          if ( JESUnc==-4  ) cand_jet.setPt( cand_jet.pt_jesFlavorQCDDown());

          // BBEC1 up
          if ( JESUnc==5   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesPileUpPtBBUp() - cand_jet.pt(),    2) +
                                                                  pow(cand_jet.pt_jesPileUpPtEC1Up() - cand_jet.pt(),   2) +
                                                                  pow(cand_jet.pt_jesRelativePtBBUp() - cand_jet.pt(),  2)));
          // BBEC1 down
          if ( JESUnc==-5  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesPileUpPtBBDown() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesPileUpPtEC1Down() - cand_jet.pt(), 2) +
                                                                  pow(cand_jet.pt_jesRelativePtBBDown() - cand_jet.pt(),2)));

          // BBEC1_year up
          if ( JESUnc==6   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesRelativeJEREC1Up() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesRelativeJEREC2Up() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesRelativeStatECUp() - cand_jet.pt(),  2)));
          // BBEC1_year down
          if ( JESUnc==-6  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesRelativeJEREC1Down() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesRelativeJEREC2Down() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesRelativeStatECDown() - cand_jet.pt(),2)));

          // EC2 up
          if ( JESUnc==7   ) cand_jet.setPt( cand_jet.pt_jesPileUpPtEC2Up()  );
          // EC2 down
          if ( JESUnc==-7  ) cand_jet.setPt( cand_jet.pt_jesPileUpPtEC2Down());

          // EC2_year up
          if ( JESUnc==8   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesRelativeJEREC2Up() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesRelativePtEC2Up() - cand_jet.pt(),   2)));
          // EC2_year down
          if ( JESUnc==-8  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesRelativeJEREC2Down() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesRelativePtEC2Down() - cand_jet.pt(), 2)));

          // HF up
          if ( JESUnc==9   ) cand_jet.setPt( cand_jet.pt() + sqrt(pow(cand_jet.pt_jesPileUpPtHFUp() - cand_jet.pt(),     2) +
                                                                  pow(cand_jet.pt_jesRelativeJERHFUp() - cand_jet.pt(),  2) +
                                                                  pow(cand_jet.pt_jesRelativePtHFUp() - cand_jet.pt(),   2)));
          // HF down
          if ( JESUnc==-9  ) cand_jet.setPt( cand_jet.pt() - sqrt(pow(cand_jet.pt_jesPileUpPtHFDown() - cand_jet.pt(),   2) +
                                                                  pow(cand_jet.pt_jesRelativeJERHFDown() - cand_jet.pt(),2) +
                                                                  pow(cand_jet.pt_jesRelativePtHFDown() - cand_jet.pt(), 2)));

          // HF_year up
          if ( JESUnc==10  ) cand_jet.setPt( cand_jet.pt_jesRelativeStatHFUp()  );
          // HF_year down
          if ( JESUnc==-10 ) cand_jet.setPt( cand_jet.pt_jesRelativeStatHFDown());

          // RelativeBal up
          if ( JESUnc==11  ) cand_jet.setPt( cand_jet.pt_jesRelativeBalUp()  );
          // RelativeBal down
          if ( JESUnc==-11 ) cand_jet.setPt( cand_jet.pt_jesRelativeBalDown());

          // RelativeSample_year up
          if ( JESUnc==12  ) cand_jet.setPt( cand_jet.pt_jesRelativeSampleUp()  );
          // RelativeSample_year down
          if ( JESUnc==-12 ) cand_jet.setPt( cand_jet.pt_jesRelativeSampleDown());
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
