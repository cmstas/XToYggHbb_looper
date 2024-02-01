#!/bin/bash

cd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_30/ ; cmsenv ; cd - >/dev/null
source ~/miniconda3/etc/profile.d/conda.sh

echo "root_to_pickle started!"
INTERMEDIATEPKLFILES=$(python python/root_to_pickle.py $1)
echo "root_to_pickle finished!"
conda activate fastparquet
echo "pickle_to_parquet started!"
python python/pickle_to_parquet.py $1
echo "pickle_to_parquet finished!"
conda deactivate

echo "Intermediate files to be removed!"
rm $INTERMEDIATEPKLFILES
