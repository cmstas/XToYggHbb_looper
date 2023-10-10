import ROOT
import numpy
import math

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", help = "File with histogram to derive shape from", type=str)
parser.add_argument("--histName", help = "Name of histogram to fit", type=str, default = "htemp")
parser.add_argument("--twoD", help = "Fit 2D function in min/max photon MVA ID (not really tested)", default = False, action='store_true')
args = parser.parse_args()

#ROOT.gROOT.SetBatch(1)

fGJets = ROOT.TFile(args.input)

hist_name = args.histName

fakePhotonID = fGJets.Get(hist_name)


def fit_func(x, par):
  if x[1] >= x[0]:
    return 0
  else:
    return par[0] + par[1]*x[0] + par[2]*x[1] + par[3]*x[0]**2 + par[4]*x[0]*x[1] + par[5]*x[1]**2 + par[6]*x[0]**3 + par[7]*x[0]**2*x[1] + par[8]*x[0]*x[1]**2 + par[9]*x[1]**3 + par[10]*x[0]**4 + par[11]*x[0]**3*x[1] + par[12]*x[0]**2*x[1]**2 + par[13]*x[0]*x[1]**3 + par[14]*x[1]**4 + par[15]*x[0]**5 + par[16]*x[0]**4*x[1] + par[17]*x[0]**3*x[1]**2 + par[18]*x[0]**2*x[1]**3 + par[19]*x[0]*x[1]**4 + par[20]*x[1]**5 + par[21]*x[0]**6 + par[22]*x[0]**5*x[1] + par[23]*x[0]**4*x[1]**2 + par[24]*x[0]**3*x[1]**3 + par[25]*x[0]**2*x[1]**4 + par[26]*x[0]*x[1]**5 + par[27]*x[1]**6 + par[28]*x[0]**7 + par[29]*x[0]**6*x[1] + par[30]*x[0]**5*x[1]**2 + par[31]*x[0]**4*x[1]**3 + par[32]*x[0]**3*x[1]**4 + par[33]*x[0]**2*x[1]**5 + par[34]*x[0]*x[1]**6 + par[36]*x[1]**7 


lowerBound = -0.9
upperBound = 1.0
if not args.twoD:
  fakePhotonID.GetXaxis().SetRangeUser(lowerBound, upperBound)
  # From skim
  #f1 = ROOT.TF1("f1", "expo(0)+pol8(2)", lowerBound, upperBound); # from skim
  # From preselection
  f1 = ROOT.TF1("f1", "pol9", lowerBound, upperBound); # from presel
  fit = (fakePhotonID.Fit("f1")).Get()
  fakePhotonID.Draw("")
else:
  fakePhotonID.GetXaxis().SetRangeUser(lowerBound, upperBound)
  fakePhotonID.GetYaxis().SetRangeUser(lowerBound, upperBound)
  func = ROOT.TF2("func", fit_func, lowerBound, upperBound, lowerBound, upperBound, 37)
  func.SetContour(100) 
  fit = (fakePhotonID.Fit("func")).Get()
  fakePhotonID.Draw("LEGO")

input('Press ENTER to exit')
