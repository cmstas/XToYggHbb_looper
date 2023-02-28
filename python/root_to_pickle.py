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
  weight_scale=1
  if "Data" not in file_name:
    file_in = ROOT.TFile.Open(file_name ,"READ")
    h_weight = file_in.Get("weight")
    h_weight_before_btagsf = file_in.Get("weight_before_btagsf")
    if (h_weight.GetBinContent(1)!=0):
      weight_scale = h_weight_before_btagsf.GetBinContent(1)/h_weight.GetBinContent(1)
  tree = uproot.open(file_name)["tout"]
  df = tree.pandas.df()
  df.loc[df.process_id!=0, "weight_central"] = df.loc[df.process_id!=0, "weight_central"]*weight_scale
  df.to_pickle(file_name.replace(".root",".pkl"))

print(output_name)
