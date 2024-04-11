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
        try { pt_jesAbsoluteMPFBiasUp_ = nt.Jet_pt_jesAbsoluteMPFBiasUp()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteMPFBiasUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesAbsoluteMPFBiasDown_ = nt.Jet_pt_jesAbsoluteMPFBiasDown()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteMPFBiasDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesAbsoluteScaleUp_ = nt.Jet_pt_jesAbsoluteScaleUp()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteScaleUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesAbsoluteScaleDown_ = nt.Jet_pt_jesAbsoluteScaleDown()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteScaleDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesAbsoluteStatUp_ = nt.Jet_pt_jesAbsoluteStatUp()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteStatUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesAbsoluteStatDown_ = nt.Jet_pt_jesAbsoluteStatDown()[idx_]; }
        catch(const std::exception& e) { pt_jesAbsoluteStatDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesFlavorQCDUp_ = nt.Jet_pt_jesFlavorQCDUp()[idx_]; }
        catch(const std::exception& e) { pt_jesFlavorQCDUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesFlavorQCDDown_ = nt.Jet_pt_jesFlavorQCDDown()[idx_]; }
        catch(const std::exception& e) { pt_jesFlavorQCDDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesFragmentationUp_ = nt.Jet_pt_jesFragmentationUp()[idx_]; }
        catch(const std::exception& e) { pt_jesFragmentationUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesFragmentationDown_ = nt.Jet_pt_jesFragmentationDown()[idx_]; }
        catch(const std::exception& e) { pt_jesFragmentationDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpDataMCUp_ = nt.Jet_pt_jesPileUpDataMCUp()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpDataMCUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpDataMCDown_ = nt.Jet_pt_jesPileUpDataMCDown()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpDataMCDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtBBUp_ = nt.Jet_pt_jesPileUpPtBBUp()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtBBUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtBBDown_ = nt.Jet_pt_jesPileUpPtBBDown()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtBBDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtEC1Up_ = nt.Jet_pt_jesPileUpPtEC1Up()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtEC1Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtEC1Down_ = nt.Jet_pt_jesPileUpPtEC1Down()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtEC1Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtEC2Up_ = nt.Jet_pt_jesPileUpPtEC2Up()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtEC2Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtEC2Down_ = nt.Jet_pt_jesPileUpPtEC2Down()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtEC2Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtHFUp_ = nt.Jet_pt_jesPileUpPtHFUp()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtHFUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtHFDown_ = nt.Jet_pt_jesPileUpPtHFDown()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtHFDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtRefUp_ = nt.Jet_pt_jesPileUpPtRefUp()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtRefUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesPileUpPtRefDown_ = nt.Jet_pt_jesPileUpPtRefDown()[idx_]; }
        catch(const std::exception& e) { pt_jesPileUpPtRefDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeFSRUp_ = nt.Jet_pt_jesRelativeFSRUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeFSRUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeFSRDown_ = nt.Jet_pt_jesRelativeFSRDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeFSRDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJEREC1Up_ = nt.Jet_pt_jesRelativeJEREC1Up()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJEREC1Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJEREC1Down_ = nt.Jet_pt_jesRelativeJEREC1Down()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJEREC1Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJEREC2Up_ = nt.Jet_pt_jesRelativeJEREC2Up()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJEREC2Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJEREC2Down_ = nt.Jet_pt_jesRelativeJEREC2Down()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJEREC2Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJERHFUp_ = nt.Jet_pt_jesRelativeJERHFUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJERHFUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeJERHFDown_ = nt.Jet_pt_jesRelativeJERHFDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeJERHFDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtBBUp_ = nt.Jet_pt_jesRelativePtBBUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtBBUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtBBDown_ = nt.Jet_pt_jesRelativePtBBDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtBBDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtEC1Up_ = nt.Jet_pt_jesRelativePtEC1Up()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtEC1Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtEC1Down_ = nt.Jet_pt_jesRelativePtEC1Down()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtEC1Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtEC2Up_ = nt.Jet_pt_jesRelativePtEC2Up()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtEC2Up_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtEC2Down_ = nt.Jet_pt_jesRelativePtEC2Down()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtEC2Down_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtHFUp_ = nt.Jet_pt_jesRelativePtHFUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtHFUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativePtHFDown_ = nt.Jet_pt_jesRelativePtHFDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativePtHFDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeBalUp_ = nt.Jet_pt_jesRelativeBalUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeBalUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeBalDown_ = nt.Jet_pt_jesRelativeBalDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeBalDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeSampleUp_ = nt.Jet_pt_jesRelativeSampleUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeSampleUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeSampleDown_ = nt.Jet_pt_jesRelativeSampleDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeSampleDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatECUp_ = nt.Jet_pt_jesRelativeStatECUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatECUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatECDown_ = nt.Jet_pt_jesRelativeStatECDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatECDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatFSRUp_ = nt.Jet_pt_jesRelativeStatFSRUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatFSRUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatFSRDown_ = nt.Jet_pt_jesRelativeStatFSRDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatFSRDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatHFUp_ = nt.Jet_pt_jesRelativeStatHFUp()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatHFUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesRelativeStatHFDown_ = nt.Jet_pt_jesRelativeStatHFDown()[idx_]; }
        catch(const std::exception& e) { pt_jesRelativeStatHFDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesSinglePionECALUp_ = nt.Jet_pt_jesSinglePionECALUp()[idx_]; }
        catch(const std::exception& e) { pt_jesSinglePionECALUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesSinglePionECALDown_ = nt.Jet_pt_jesSinglePionECALDown()[idx_]; }
        catch(const std::exception& e) { pt_jesSinglePionECALDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesSinglePionHCALUp_ = nt.Jet_pt_jesSinglePionHCALUp()[idx_]; }
        catch(const std::exception& e) { pt_jesSinglePionHCALUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesSinglePionHCALDown_ = nt.Jet_pt_jesSinglePionHCALDown()[idx_]; }
        catch(const std::exception& e) { pt_jesSinglePionHCALDown_ = nt.Jet_pt()[idx_]; }
        try { pt_jesTimePtEtaUp_ = nt.Jet_pt_jesTimePtEtaUp()[idx_]; }
        catch(const std::exception& e) { pt_jesTimePtEtaUp_ = nt.Jet_pt()[idx_]; }
        try { pt_jesTimePtEtaDown_ = nt.Jet_pt_jesTimePtEtaDown()[idx_]; }
        catch(const std::exception& e) { pt_jesTimePtEtaDown_ = nt.Jet_pt()[idx_]; }
        mass_ = nt.Jet_mass()[idx_];
        eta_ = nt.Jet_eta()[idx_];
        phi_ = nt.Jet_phi()[idx_];
        jetId_ = nt.Jet_jetId()[idx_];
        puId_ = nt.Jet_puId()[idx_];
        p4_.SetPtEtaPhiM(pt_, nt.Jet_eta()[idx_], nt.Jet_phi()[idx_], nt.Jet_mass()[idx_]);
        try {hadronFlavour_ = nt.Jet_hadronFlavour()[idx_];}
        catch(const std::exception& e) {hadronFlavour_=-1;}
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
    float pt_jesAbsoluteMPFBiasUp() { return pt_jesAbsoluteMPFBiasUp_; }
    float pt_jesAbsoluteMPFBiasDown() { return pt_jesAbsoluteMPFBiasDown_; }
    float pt_jesAbsoluteScaleUp() { return pt_jesAbsoluteScaleUp_; }
    float pt_jesAbsoluteScaleDown() { return pt_jesAbsoluteScaleDown_; }
    float pt_jesAbsoluteStatUp() { return pt_jesAbsoluteStatUp_; }
    float pt_jesAbsoluteStatDown() { return pt_jesAbsoluteStatDown_; }
    float pt_jesFlavorQCDUp() { return pt_jesFlavorQCDUp_; }
    float pt_jesFlavorQCDDown() { return pt_jesFlavorQCDDown_; }
    float pt_jesFragmentationUp() { return pt_jesFragmentationUp_; }
    float pt_jesFragmentationDown() { return pt_jesFragmentationDown_; }
    float pt_jesPileUpDataMCUp() { return pt_jesPileUpDataMCUp_; }
    float pt_jesPileUpDataMCDown() { return pt_jesPileUpDataMCDown_; }
    float pt_jesPileUpPtBBUp() { return pt_jesPileUpPtBBUp_; }
    float pt_jesPileUpPtBBDown() { return pt_jesPileUpPtBBDown_; }
    float pt_jesPileUpPtEC1Up() { return pt_jesPileUpPtEC1Up_; }
    float pt_jesPileUpPtEC1Down() { return pt_jesPileUpPtEC1Down_; }
    float pt_jesPileUpPtEC2Up() { return pt_jesPileUpPtEC2Up_; }
    float pt_jesPileUpPtEC2Down() { return pt_jesPileUpPtEC2Down_; }
    float pt_jesPileUpPtHFUp() { return pt_jesPileUpPtHFUp_; }
    float pt_jesPileUpPtHFDown() { return pt_jesPileUpPtHFDown_; }
    float pt_jesPileUpPtRefUp() { return pt_jesPileUpPtRefUp_; }
    float pt_jesPileUpPtRefDown() { return pt_jesPileUpPtRefDown_; }
    float pt_jesRelativeFSRUp() { return pt_jesRelativeFSRUp_; }
    float pt_jesRelativeFSRDown() { return pt_jesRelativeFSRDown_; }
    float pt_jesRelativeJEREC1Up() { return pt_jesRelativeJEREC1Up_; }
    float pt_jesRelativeJEREC1Down() { return pt_jesRelativeJEREC1Down_; }
    float pt_jesRelativeJEREC2Up() { return pt_jesRelativeJEREC2Up_; }
    float pt_jesRelativeJEREC2Down() { return pt_jesRelativeJEREC2Down_; }
    float pt_jesRelativeJERHFUp() { return pt_jesRelativeJERHFUp_; }
    float pt_jesRelativeJERHFDown() { return pt_jesRelativeJERHFDown_; }
    float pt_jesRelativePtBBUp() { return pt_jesRelativePtBBUp_; }
    float pt_jesRelativePtBBDown() { return pt_jesRelativePtBBDown_; }
    float pt_jesRelativePtEC1Up() { return pt_jesRelativePtEC1Up_; }
    float pt_jesRelativePtEC1Down() { return pt_jesRelativePtEC1Down_; }
    float pt_jesRelativePtEC2Up() { return pt_jesRelativePtEC2Up_; }
    float pt_jesRelativePtEC2Down() { return pt_jesRelativePtEC2Down_; }
    float pt_jesRelativePtHFUp() { return pt_jesRelativePtHFUp_; }
    float pt_jesRelativePtHFDown() { return pt_jesRelativePtHFDown_; }
    float pt_jesRelativeBalUp() { return pt_jesRelativeBalUp_; }
    float pt_jesRelativeBalDown() { return pt_jesRelativeBalDown_; }
    float pt_jesRelativeSampleUp() { return pt_jesRelativeSampleUp_; }
    float pt_jesRelativeSampleDown() { return pt_jesRelativeSampleDown_; }
    float pt_jesRelativeStatECUp() { return pt_jesRelativeStatECUp_; }
    float pt_jesRelativeStatECDown() { return pt_jesRelativeStatECDown_; }
    float pt_jesRelativeStatFSRUp() { return pt_jesRelativeStatFSRUp_; }
    float pt_jesRelativeStatFSRDown() { return pt_jesRelativeStatFSRDown_; }
    float pt_jesRelativeStatHFUp() { return pt_jesRelativeStatHFUp_; }
    float pt_jesRelativeStatHFDown() { return pt_jesRelativeStatHFDown_; }
    float pt_jesSinglePionECALUp() { return pt_jesSinglePionECALUp_; }
    float pt_jesSinglePionECALDown() { return pt_jesSinglePionECALDown_; }
    float pt_jesSinglePionHCALUp() { return pt_jesSinglePionHCALUp_; }
    float pt_jesSinglePionHCALDown() { return pt_jesSinglePionHCALDown_; }
    float pt_jesTimePtEtaUp() { return pt_jesTimePtEtaUp_; }
    float pt_jesTimePtEtaDown() { return pt_jesTimePtEtaDown_; }
    float mass() { return mass_; }
    float eta() { return eta_; }
    float phi() { return phi_; }
    int jetId() { return jetId_; }
    int puId() { return puId_; }
    int hadronFlavour() { return  hadronFlavour_; }
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
    float pt_jesAbsoluteMPFBiasUp_ = 0.;
    float pt_jesAbsoluteMPFBiasDown_ = 0.;
    float pt_jesAbsoluteScaleUp_ = 0.;
    float pt_jesAbsoluteScaleDown_ = 0.;
    float pt_jesAbsoluteStatUp_ = 0.;
    float pt_jesAbsoluteStatDown_ = 0.;
    float pt_jesFlavorQCDUp_ = 0.;
    float pt_jesFlavorQCDDown_ = 0.;
    float pt_jesFragmentationUp_ = 0.;
    float pt_jesFragmentationDown_ = 0.;
    float pt_jesPileUpDataMCUp_ = 0.;
    float pt_jesPileUpDataMCDown_ = 0.;
    float pt_jesPileUpPtBBUp_ = 0.;
    float pt_jesPileUpPtBBDown_ = 0.;
    float pt_jesPileUpPtEC1Up_ = 0.;
    float pt_jesPileUpPtEC1Down_ = 0.;
    float pt_jesPileUpPtEC2Up_ = 0.;
    float pt_jesPileUpPtEC2Down_ = 0.;
    float pt_jesPileUpPtHFUp_ = 0.;
    float pt_jesPileUpPtHFDown_ = 0.;
    float pt_jesPileUpPtRefUp_ = 0.;
    float pt_jesPileUpPtRefDown_ = 0.;
    float pt_jesRelativeFSRUp_ = 0.;
    float pt_jesRelativeFSRDown_ = 0.;
    float pt_jesRelativeJEREC1Up_ = 0.;
    float pt_jesRelativeJEREC1Down_ = 0.;
    float pt_jesRelativeJEREC2Up_ = 0.;
    float pt_jesRelativeJEREC2Down_ = 0.;
    float pt_jesRelativeJERHFUp_ = 0.;
    float pt_jesRelativeJERHFDown_ = 0.;
    float pt_jesRelativePtBBUp_ = 0.;
    float pt_jesRelativePtBBDown_ = 0.;
    float pt_jesRelativePtEC1Up_ = 0.;
    float pt_jesRelativePtEC1Down_ = 0.;
    float pt_jesRelativePtEC2Up_ = 0.;
    float pt_jesRelativePtEC2Down_ = 0.;
    float pt_jesRelativePtHFUp_ = 0.;
    float pt_jesRelativePtHFDown_ = 0.;
    float pt_jesRelativeBalUp_ = 0.;
    float pt_jesRelativeBalDown_ = 0.;
    float pt_jesRelativeSampleUp_ = 0.;
    float pt_jesRelativeSampleDown_ = 0.;
    float pt_jesRelativeStatECUp_ = 0.;
    float pt_jesRelativeStatECDown_ = 0.;
    float pt_jesRelativeStatFSRUp_ = 0.;
    float pt_jesRelativeStatFSRDown_ = 0.;
    float pt_jesRelativeStatHFUp_ = 0.;
    float pt_jesRelativeStatHFDown_ = 0.;
    float pt_jesSinglePionECALUp_ = 0.;
    float pt_jesSinglePionECALDown_ = 0.;
    float pt_jesSinglePionHCALUp_ = 0.;
    float pt_jesSinglePionHCALDown_ = 0.;
    float pt_jesTimePtEtaUp_ = 0.;
    float pt_jesTimePtEtaDown_ = 0.;
    float eta_ = 0.;
    float mass_ = 0;
    TLorentzVector p4_;
    float phi_ = 0.;
    int jetId_ = 0;
    int puId_ = 0;
    unsigned int idx_;
    float btagDeepFlavB_ = 0.;
    int  hadronFlavour_ = -1;
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

vector<Jet> getJets(Photons photons, const int JESUnc, const int JERUnc, const int HEMCheck);
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
