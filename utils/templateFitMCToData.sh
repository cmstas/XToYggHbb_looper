#!/bin/bash

IN=$1
OUT=$2
while ! [ -z "$3" ]; do
    FLAGS="$FLAGS $3"; shift;
done

# Data
python python/tree_plotting.py --inDir $IN --outDir $OUT --dataOnly --outFormat root --label _data_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --dataOnly --outFormat root --label _data_weighted --DDHistos
# QCD
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples QCD_Pt-30to40_MGG-80toInf QCD_Pt-30toInf_MGG-40to80 QCD_Pt-40toInf_MGG-80toInf --outFormat root --label _ff_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples QCD_Pt-30to40_MGG-80toInf QCD_Pt-30toInf_MGG-40to80 QCD_Pt-40toInf_MGG-80toInf --outFormat root --label _ff_weighted --DDHistos
# GJets
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples GJets --outFormat root --label _fp_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples GJets --outFormat root --label _fp_weighted --DDHistos
# DD QCD+GJets
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DDQCDGJets --outFormat root --label _fffp_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DDQCDGJets --outFormat root --label _fffp_weighted --DDHistos
# DiPhoton+DiPhotonLow
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhoton DiPhotonLow --outFormat root --label _pp_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhoton DiPhotonLow --outFormat root --label _pp_weighted --DDHistos
# DiPhotonLow
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhotonLow --outFormat root --label _ppL_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhotonLow --outFormat root --label _ppL_weighted --DDHistos
# DiPhoton
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhoton --outFormat root --label _ppH_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples DiPhoton --outFormat root --label _ppH_weighted --DDHistos
# Bkg
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples TTGG TTGJets TTJets DY WG ZG WW WZ ZZ VBFH_M125 VH_M125 ggHToDiPhoM125 ttH_M125 HHbbgg --outFormat root --label _bkg_unweighted --unweighted --DDHistos
python python/tree_plotting.py --inDir $IN --outDir $OUT --noSignal --samples TTGG TTGJets TTJets DY WG ZG WW WZ ZZ VBFH_M125 VH_M125 ggHToDiPhoM125 ttH_M125 HHbbgg --outFormat root --label _bkg_weighted --DDHistos

python python/do_fits_qcd.py --inputDir $OUT $FLAGS
