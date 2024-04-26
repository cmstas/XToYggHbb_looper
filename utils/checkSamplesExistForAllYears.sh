#!/bin/bash

DIR=$1; shift;
if [ -z $1 ]; then
  samples=(
    "Data"
    "DiPhoton"
    "TTGG"
    "TTGJets"
    "TTJets"
    "VH_M125"
    "VBFH_M125"
    "ttH_M125"
    "ggHToDiPhoM125"
    "GJets_HT-40To100"
    "GJets_HT-100To200"
    "GJets_HT-200To400"
    "GJets_HT-400To600"
    "GJets_HT-600ToInf"
    "DY"
    "WG"
    "ZG"
    "HHbbgg"
    "DiPhotonLow"
    "WW"
    "WZ"
    "ZZ"
    "QCD_Pt-30to40_MGG-80toInf"
    "QCD_Pt-30toInf_MGG-40to80"
    "QCD_Pt-40toInf_MGG-80toInf"
    "DDQCDGJets"
    "tZq"
    "TTZ"
    "TTW"
    "NMSSM_XToYHTo2G2B_MX_1000_MY_90" # Signal_low
    "NMSSM_XToYHTo2G2B_MX_500_MY_90" # Signal_med
    "NMSSM_XToYHTo2G2B_MX_900_MY_90" # Signal_high
    "NMSSM_XToYHTo2B2G_MX_650_MY_70" # Inverted signal
  );
elif [ $1 = '--reduced' ]; then
  samples=(
    "VH_M125"
    "VBFH_M125"
    "ttH_M125"
    "ggHToDiPhoM125"
    "HHbbgg"
    "NMSSM_XToYHTo2G2B_MX_1000_MY_90" # Signal_low
    "NMSSM_XToYHTo2G2B_MX_500_MY_90" # Signal_med
    "NMSSM_XToYHTo2G2B_MX_900_MY_90" # Signal_high
  );
else
  echo "Unknown option!";
fi

for s in "${samples[@]}"; do
  if [ $(ls ${DIR}"output_"${s}$"_"*".root" -1 | wc -l) -ne 4 ]; then echo "Missing file for $s"; fi;
done
