import pandas as pd
import uproot
import sys
import glob
import ROOT

output_name = sys.argv[1].replace(".root",".pkl")
file_names = [ sys.argv[1] ]
if sys.argv[1].endswith("/"):
  output_name = (sys.argv[1]+"*.root").replace(".root",".pkl")
  file_names = glob.glob(sys.argv[1]+"*.root")

df = pd.DataFrame()
for file_name in file_names:
  weight_scale = 1.0
  if ("Data" not in file_name) and ("DDQCDGJets" not in file_name):
    fileIn = ROOT.TFile.Open(file_name,"READ")
    h_beforeBTagSF = fileIn.Get("weight_beforeBTagSF")
    beforeBTagSF = h_beforeBTagSF.GetBinContent(1)
    h_afterBTagSF = fileIn.Get("weight_afterBTagSF")
    afterBTagSF = h_afterBTagSF.GetBinContent(1)
    if beforeBTagSF == 0:
      print("Zero sum of weights for file %s..."%file_name)
    else:
      weight_scale = beforeBTagSF / afterBTagSF
  tree = uproot.open(file_name)["tout"]
  df = tree.pandas.df()
  df.loc[df.process_id!=0,"weight_central"] = df.loc[df.process_id!=0,"weight_central"]*weight_scale
  df.to_pickle(file_name.replace(".root",".pkl"))

print(output_name)
