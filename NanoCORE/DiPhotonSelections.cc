#include "DiPhotonSelections.h"
#include "Tools/fnufUnc.h"
#include "Tools/materialUnc.h"

using namespace tas;

bool UseLowR9Photon(Photon pho, bool isEB) {
    bool useThisPhoton = false;
    bool loweta = abs(pho.eta())<1.5;
    if (isEB) {
        if ( !(pho.sieie() < 0.015) ) return useThisPhoton;       
        if ( !(pho.trkSumPtHollowConeDR03() < 6.0) ) return useThisPhoton;
        if ( loweta && !(pho.phoIso() - 0.16544*pho.perEvtRho() < 4.0) ) return useThisPhoton;       
        if ( !(loweta) && !(pho.phoIso() - 0.13212*pho.perEvtRho() < 4.0) ) return useThisPhoton;
    } else {
        if ( !(pho.sieie() < 0.035) ) return useThisPhoton;       
        if ( !(pho.trkSumPtHollowConeDR03() < 6.0) ) return useThisPhoton;
        if ( loweta && !(pho.phoIso() - 0.16544*pho.perEvtRho() < 4.0) ) return useThisPhoton;       
        if ( !(loweta) && !(pho.phoIso() - 0.13212*pho.perEvtRho() < 4.0) ) return useThisPhoton;      
    }
    // 0.16544 and 0.13212 are copied from flashggPreselectedDiPhotons_cfi.py
    useThisPhoton = true;
    return useThisPhoton;
}

Photons getPhotons(const TString year, const int fnufUnc, const int materialUnc, const int PhoScaleUnc, const int PhoSmearUnc) {
    Photons photons;

    if ( fnufUnc != 0) fnufUnc::set_fnufUnc();
    if ( materialUnc != 0) materialUnc::set_materialUnc();

    for (unsigned int ipho = 0; ipho < nt.nPhoton(); ipho++) {
        Photon pho = Photon(ipho);

        // Apply fnuf unc
        if ( fnufUnc!=0 ) {
          // fnufUnc is only an uncertainty => fnufUnc==1 leads to dummy multiplication by unity
          if ( fnufUnc==2  ) pho.setPt( pho.pt()*fnufUnc::get_fnufUnc(pho.eta(), pho.r9(), year, "up") );
          if ( fnufUnc==-2 ) pho.setPt( pho.pt()*fnufUnc::get_fnufUnc(pho.eta(), pho.r9(), year, "down") );
        }
        // Apply material unc
        if ( materialUnc!=0 ) {
          // materialUnc is only an uncertainty => materialUnc==1 leads to dummy multiplication by unity
          if ( materialUnc==2  ) pho.setPt( pho.pt()*materialUnc::get_materialUnc(pho.eta(), pho.r9(), year, "up") );
          if ( materialUnc==-2 ) pho.setPt( pho.pt()*materialUnc::get_materialUnc(pho.eta(), pho.r9(), year, "down") );
        }
        // Apply photon scale unc
        if ( PhoScaleUnc!=0 ) {
          // PhoScaleUnc is only an uncertainty => PhoScaleUnc==1 leads to dummy multiplication by unity
          if ( PhoScaleUnc==2  ) pho.setPt( pho.pt_ScaleUp() );
          if ( PhoScaleUnc==-2 ) pho.setPt( pho.pt_ScaleDown() );
        }
        // Apply photon smear unc
        if ( PhoSmearUnc!=0 ) {
          // PhoSmearUnc is only an uncertainty => PhoSmearUnc==1 leads to dummy multiplication by unity
          if ( PhoSmearUnc==2  ) pho.setPt( pho.pt()+pho.dEsigmaUp() );
          if ( PhoSmearUnc==-2 ) pho.setPt( pho.pt()+pho.dEsigmaDown() );
        }

        if ( !(pho.pt()>18) ) continue;
        if ( !(pho.isScEtaEB() || pho.isScEtaEE()) ) continue;
        if ( !(pho.hoe()<0.08) ) continue;
        //if ( pho.pixelSeed() > 0.5 ) continue; // Not standard photon selection, but used to suppress electrons on the DY peak
        // inverted 
        if ( !(pho.pixelSeed() <= 0.5) ) continue;
        if ( pho.eveto() < 0.5 ) continue;

        if ( !(pho.r9() > 0.8 || pho.chargedHadIso() < 20 || pho.chargedHadIso()/pho.pt() < 0.3) ) continue;

        bool pho_EB_highR9 = pho.isScEtaEB() && pho.r9() > 0.85; 
        bool pho_EB_lowR9 = pho.isScEtaEB() && pho.r9() < 0.85 && pho.r9() > 0.50 && UseLowR9Photon(pho, true);
        bool pho_EE_highR9 = pho.isScEtaEE() && pho.r9() > 0.90; 
        bool pho_EE_lowR9 = pho.isScEtaEE() && pho.r9() < 0.90 && pho.r9() > 0.80 && UseLowR9Photon(pho, false);
        if ( !(pho_EB_highR9 || pho_EE_highR9 || pho_EB_lowR9 || pho_EE_lowR9) ) continue;
        photons.push_back(pho);
    }

    sort(photons.begin(), photons.end(), sortByPt);

    if ( fnufUnc != 0) fnufUnc::reset_fnufUnc();
    if ( materialUnc != 0) materialUnc::reset_materialUnc();

    return photons;
}

DiPhotons DiPhotonPreselection(Photons &photons) {
    DiPhotons diphotons; 
    float maxSumDiphoPt = 0;    
    for (unsigned int i1 = 0; i1 < photons.size(); i1++) {
        for (unsigned int i2 = i1+1; i2 < photons.size(); i2++) {
            Photon pho1 = photons[i1];
            Photon pho2 = photons[i2];
            DiPhoton dipho = DiPhoton(pho1, pho2);
            if ( !(dipho.leadPho.pt() >= 30.0 && dipho.subleadPho.pt() > 18.0) ) continue;        
//            if ( !(dipho.leadPho.pt()/dipho.p4.M() > 0.33 && dipho.subleadPho.pt()/dipho.p4.M() > 0.25) ) continue; 
//it's a historical cut ,when searching for Higgs, people don't know the dipho mass, then for higher mass they require tighter pt
//also, trigger cut around 20/30 , bkg mgg shape is not a smooth falling distribution ~100GeV (need to be tested)

            if ( !(dipho.p4.M() >= 55 && dipho.p4.M() <= 999999) ) continue;

            float sumDiPhoPt = dipho.leadPho.pt() + dipho.subleadPho.pt();

            if ( diphotons.size() == 0 )
            {
                diphotons.push_back(dipho);
                maxSumDiphoPt = sumDiPhoPt;
            } 
            else if (sumDiPhoPt > maxSumDiphoPt) 
            {
                diphotons.pop_back();
                diphotons.push_back(dipho);
                maxSumDiphoPt = sumDiPhoPt;
            }
        }
    } 

    // thre could be either 1 or 0 diphoton in the output
    return diphotons;    
}
