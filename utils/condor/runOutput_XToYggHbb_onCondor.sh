#!/bin/bash

usage()
{
    echo "Usage:"
    echo ""
    echo "  sh utils/condor/runOutput_XToYggHbb_onCondor.sh output_dir"
    echo ""
    echo "The output_dir will be created in /ceph/cms/store/user/$USER/XToYggHbbOutput/"
    echo "Control the jobs to be run by editing runOutput_XToYggHbb_onCondor.sub"
    echo ""
    exit
}

if [ -z $1 ]; then usage; fi

export OUTPUTDIR=$1; shift;
export STARTDIR=$PWD

if [ -z $1 ]; then
  sh utils/condor/create_package.sh
elif [ $1 = '--soft' ]; then
  echo 'Package will not be created due to "--soft" option.'
else
  echo "Unknown option!";
fi

mkdir -p utils/condor/plotting_logs
mkdir -p /ceph/cms/store/user/$USER/XToYggHbbOutput/$OUTPUTDIR
cp cpp/summary.json /ceph/cms/store/user/$USER/XToYggHbbOutput/$OUTPUTDIR/.

condor_submit utils/condor/runOutput_XToYggHbb_onCondor.sub
