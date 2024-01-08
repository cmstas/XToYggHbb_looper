# X -> YH -> γγbb instructions

## Environment & setup

```bash
git clone git@github.com:cmstas/XToYggHbb_looper.git
cd XToYggHbb_looper/
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_2_13/src ; eval `scramv1 runtime -sh` ; cd -
cd NanoCORE
make -j8
cd -
```

**N.B.:** The compilation of the NanoCORE package needs to be repeated each time a file inside the package is modified.

## Preselection code

Edit `cpp/ScanChain_Hgg.C` and/or `cpp/main.cc` appropriately. The former contains the looper that defines the analysis preselection, while the latter contains the logic that includes the samples to be run. Before starting working with the code, the following command needs to be run:

```bash
source setup.sh
```

This makes ROOT available and initiates the grid certificate (hence the prompt for a code).

### Compilation

Any time any of the above files is edited, the code needs to be compiled. This is achieved by:

```bash
cp cpp/
make clean
make -j4
cd -
```

This creates the cpp/main.exe executable that runs the preselection code.

### Running

**To run locally:**

After the code is compiled, it can be run locally with the following command:

```bash
cd cpp
./main.exe OUTPUTDIRNAME YEAR RUNDATA RUNBKG RUNSIGNAL SAMPLE ADDITIONALBOOLEANFLAGS
cd -
```

For example, to run all data, bkg and signal, for all years and with all additional flags to their default values, the command to save the files in a folder called "temp_data" would be:

```bash
cd cpp
./main.exe temp_data all 1 1 1 all
cd -
```

To run the script to generate the Drell-Yan enriched data needed for the ABCD method use the command
```bash
./main.exe <output_dir> all 1 0 0 enrichDY
```
This will run only the "Data", but using the path to the appropriate input files

One should check `cpp/main.cc` for details on the (additional) arguments and their meaning.

This loops and creates a number of output files of the form `output_"process"_"year".root` containing histograms. 

**To run on condor:**

To run on condor, the grid certificate needs to be active and export to the correct environmental variable. This is all done by the `setup.sh` script. Then, the jobs can be sent to condor by running the following command:

```bash
sh utils/condor/runOutput_XToYggHbb_onCondor.sh FOLDER/FOR/OUTPUT/FILES
```

This script will package the current state of the repository and send it to condor jobs running the `cpp/runOutput_XToYggHbb.sh` script with the arguments included in the different lines of `utils/condor/runOutput_XToYggHbb_onCondor.sub`.

*Please edit the latter file to control what condor jobs you send.*

The output of your jobs will be found under `/ceph/cms/store/user/$USER/XToYggHbbOutput/FOLDER/FOR/OUTPUT/FILES` and the plotting logs under `utils/condor/plotting_logs`.

**To produce plots:**

The `python/tree_plotting` script has been added to produce plots from the trees containing the preselected events.
Apart from the command line options, which can be shown by running `python python/tree_plotting -h`, some parameters are controlled from within the script:
- The `samples` list contains the list of samples to run on.
- The `weight` string can be used to multiplicatively scale the events.
- The `cut` string is used to apply additional selections to the plots. The selection can be formed by (combinations of) existing branches, using C++ syntax.
- The plots, which are created using the ROOT `TTree::Draw` as backend, are defined by:
  - The `plotNames` list, which contains the (combination of) branches to be plotted. This is also the name of the plot.
  - The `plotBins` dictionary which contains the binning definition of TH1 or TH2 of ROOT, either with fixed or variable binning.
  - The `plotXTitles` dictionary which is the x(-y) axis(axes) title.

**To produce cutflow table:**

A full cutflow table is still: :construction: **WIP** :construction:

However, a final yield printer has been incorporated in the plotting script and can be run by enabling the `--yields` flag.

### Data-driven QCD+GJets background estimation

Run the preselection on GJets, disabling the MVA ID cut (comment out and recompile). After the looper has finished, move the output files in the proper directories, as per in the `cpp/getFakePhotonsFromGJetsFromPresel.cpp`. Then, run it with:
```
cd cpp
root -l -b -q getFakePhotonsFromGJetsFromPresel.cpp
cd -
```

The fit parameters are extracted by running:
```
python python/derive_impute_shape.py --input cpp/fakePhotonsFromGJetsFromPresel.root
```
and the numerical values are then inserted in the appropriate function of the looper script. Once this is done, one can run the `DDQCDGJets` sample.

An extra step (not applied in the current analysis) is the fitting of the fake-fake, fake-prompt and prompt-prompt photon background processes to the data yields. This can be done by running:
```
ls utils/templateFitMCToData.sh /dir/with/input/trees /dir/for/output/files extraFlags
```
where the extra flags are the flags of the `python/do_fits_qcd.py` script.

### Converting .root files to .parquet files

The output of the preselection code is a list of .root files. These output files are meant to be the input of the analysis pNN, which expects a single .parquet file for the merged output. This can be done by running the `root_to_parquet.sh` script. However, running this script requires setting up a virtual environment with the proper python package to read/write parquet files.

**Setup**
```bash
# download conda installer
curl -O -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b 

# add conda to the end of ~/.bashrc, so relogin after executing this line
~/miniconda3/bin/conda init

# stop conda from activating the base environment on login
conda config --set auto_activate_base false
conda config --add channels conda-forge

# create environment
conda create --name fastparquet fastparquet
```

**Running**

Once the above commands are setup correctly once, the `root_to_parquet.sh` script can be run by providing the *absolute path* to a single .root file to be converted or to a whole directory for all the .root files in it to be converted to a single, merged .parquet files. For example, in the former case:

```bash
sh utils/root_to_parquet.sh /home/users/$USER/XToYggHbb_looper/cpp/temp_data/output_DY_2018.root
```

or, the latter case:

```bash
sh utils/root_to_parquet.sh /home/users/$USER/XToYggHbb_looper/cpp/temp_data/
```

## Development

The `main` is protected from pushing commits directly to it. For new developments, a new branch needs to be created and then a PR needs to be made to the `main` branch.

