import os, sys, copy
import ROOT
import numpy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--inputDir", required=True, help = "input directory containing root files with histograms to fit", type = str)
parser.add_argument("--jetBin", help = "n_jets bin to do fit in", type = str, default = "2+")
parser.add_argument("--DD", help = "perform fits on DD estimation", action = 'store_true', default = False)
parser.add_argument("--sepDiPhotonLow", help = "separate DiPhoton from DiPhotonLowMass contribution in the fit", action = 'store_true', default = False)
args = parser.parse_args()

def find_jet_bin(jet_bin):
  return jet_bin + 1

def combine_hists(hist1, hist2, name, jet_bin, overflow):
  #print name
  target_hist = ROOT.TH1D(name, "", hist1.GetNbinsX() + hist2.GetNbinsX(), 0, 1)
  jet_bin_number = find_jet_bin(jet_bin)
  for i in range(hist1.GetNbinsX()):
    target_hist.SetBinContent(i+1, hist1.GetBinContent(i+1, jet_bin_number))
    if overflow:
      for j in range(jet_bin_number, hist1.GetNbinsY()):
        target_hist.AddBinContent(i+1, hist1.GetBinContent(i+1, j+1))
  for i in range(hist2.GetNbinsX()):
    target_hist.SetBinContent(hist1.GetNbinsX()+i+1, hist2.GetBinContent(i+1, jet_bin_number))
    if overflow:
      for j in range(jet_bin_number, hist2.GetNbinsY()):
        target_hist.AddBinContent(hist1.GetNbinsX()+i+1, hist2.GetBinContent(i+1, j+1))
  for i in range(hist1.GetNbinsX() + hist2.GetNbinsX()):
    if target_hist.GetBinContent(i+1) < 0:
      target_hist.SetBinContent(i+1, 0)

  return target_hist


latex_dict = { "DiPhoton"     : "$ \\gamma \\gamma $ + jets",
               "DiPhotonLow"  : "$ \\gamma \\gamma $ + jets (low mass)",
               "DiPhotonHigh" : "$ \\gamma \\gamma $ + jets (high mass)",
               "GJets"        : "$ \\gamma $ + jets",
               "QCD"          : "QCD",
               "DDQCDGJets"   : "QCD, $ \\gamma $ + jets"
}

pp_template = "DiPhoton"
ppL_template = "DiPhotonLow"
ppH_template = "DiPhotonHigh"
fp_template = "GJets"
ff_template = "QCD"
fffp_template = "DDQCDGJets"

file_names = ["n_jets:Diphoton_minMvaID", "n_jets:Diphoton_maxMvaID"]

templates_dict = { "data" : "Data",
		               "pp"   : "DiPhoton",
		               "ppL"  : "DiPhotonLow",
		               "ppH"  : "DiPhotonHigh",
		               "fp"   : "GJets",
		               "ff"   : "QCD",
		               "fffp" : "DDQCDGJets",
		               "bkg"  : "OtherBkg"
}

hists = {}
hists_unweighted = {}; hists_unweighted_minMVAID = {}; hists_unweighted_maxMVAID = {}
hists_weighted = {}; hists_weighted_minMVAID = {}; hists_weighted_maxMVAID = {}
hists_weights = {}

overflow = "+" in args.jetBin
jet_bin = int(args.jetBin.replace("+",""))

for file_name in file_names:
  for template, name in templates_dict.iteritems():
    if args.DD:
      if template == "ff" or template == "fp": continue
    else:
      if template == "fffp": continue
    if args.sepDiPhotonLow:
      if template == "pp": continue
    else:
      if template == "ppL" or template == "ppH": continue
    f_unweighted = ROOT.TFile(args.inputDir+"/"+file_name+"_"+template+"_unweighted.root")
    if "min" in file_name: hists_unweighted_minMVAID[template] = copy.deepcopy(f_unweighted.Get("totalSM"))
    elif "max" in file_name: hists_unweighted_maxMVAID[template] = copy.deepcopy(f_unweighted.Get("totalSM"))

    f_weighted = ROOT.TFile(args.inputDir+"/"+file_name+"_"+template+"_weighted.root")
    if "min" in file_name: hists_weighted_minMVAID[template] = copy.deepcopy(f_weighted.Get("totalSM"))
    elif "max" in file_name: hists_weighted_maxMVAID[template] = copy.deepcopy(f_weighted.Get("totalSM"))

for template, name in templates_dict.iteritems():
  if args.DD:
    if template == "ff" or template == "fp": continue
  else:
    if template == "fffp": continue
  if args.sepDiPhotonLow:
    if template == "pp": continue
  else:
    if template == "ppL" or template == "ppH": continue
  hists_unweighted[template] = combine_hists(hists_unweighted_minMVAID[template], hists_unweighted_maxMVAID[template], template + "unweighted", jet_bin, overflow)
  hists_weighted[template] = combine_hists(hists_weighted_minMVAID[template], hists_weighted_maxMVAID[template], template + "weighted", jet_bin, overflow)
  hists_weights[template] = hists_weighted[template].Clone(template + "weights")
  hists_weights[template].Divide(hists_unweighted[template])
  for i in range(hists_weights[template].GetNbinsX()):
    if hists_weights[template].GetBinContent(i+1) <= 0:
       hists_weights[template].SetBinContent(i+1, 0.000000001)
  hists[template] = [ hists_unweighted[template], hists_weighted[template], hists_weights[template] ]


hist_unweighted_list = []; hist_weighted_list = []; hist_weights_list = []
if args.DD:
  hist_unweighted_list.append(hists["fffp"][0])
  hist_weighted_list.append(hists["fffp"][1])
  hist_weights_list.append(hists["fffp"][2])
else:
  hist_unweighted_list.append(hists["ff"][0]); hist_unweighted_list.append(hists["fp"][0])
  hist_weighted_list.append(hists["ff"][1]); hist_weighted_list.append(hists["fp"][1])
  hist_weights_list.append(hists["ff"][2]); hist_weights_list.append(hists["fp"][2])
if args.sepDiPhotonLow:
  hist_unweighted_list.append(hists["ppL"][0]); hist_unweighted_list.append(hists["ppH"][0])
  hist_weighted_list.append(hists["ppL"][1]); hist_weighted_list.append(hists["ppH"][1])
  hist_weights_list.append(hists["ppL"][2]); hist_weights_list.append(hists["ppH"][2])
else:
  hist_unweighted_list.append(hists["pp"][0])
  hist_weighted_list.append(hists["pp"][1])
  hist_weights_list.append(hists["pp"][2])
hist_unweighted_list.append(hists["bkg"][0])
hist_weighted_list.append(hists["bkg"][1])
hist_weights_list.append(hists["bkg"][2])


h_data =  hists["data"][1]
initial_fracs = []
for hist in hist_weighted_list:
  initial_fracs.append(hist.Integral() / h_data.Integral())

# Fit
mc = ROOT.TObjArray(4)
if args.DD:
  mc.Add(hists["fffp"][0])
else:
  mc.Add(hists["ff"][0])
  mc.Add(hists["fp"][0])
if args.sepDiPhotonLow:
  mc.Add(hists["ppL"][0])
  mc.Add(hists["ppH"][0])
else:
  mc.Add(hists["pp"][0])
mc.Add(hists["bkg"][0])

fit = ROOT.TFractionFitter(h_data, mc, "Q")
i_range = 3 if args.DD else 4
if args.sepDiPhotonLow: i_range += 1
for i in range(i_range):
  fit.SetWeight(i, hist_weights_list[i]) # set bin-by-bin weights for raw MC counts 
  # Other bkgs, fixed in fit
  i_bkg = 2 if args.DD else 3
  if args.sepDiPhotonLow: i_bkg += 1
  if i == i_bkg:
    fit.Constrain(i, initial_fracs[i]-0.000000001, initial_fracs[i]+0.000000001)
  else:
    fit.Constrain(i, 0.0, 1.0)

fit.Fit()

fracs = []
frac_errs = []
scales = []
for i in range(i_range):
  f = ROOT.Double()
  ferr = ROOT.Double()
  fit.GetResult(i, f, ferr)
  fracs.append(f)
  frac_errs.append(ferr)
  scales.append(f/initial_fracs[i])

# Print fit results
print "N_jets == %s" % args.jetBin 
print "Initial fractions", initial_fracs
print "Fitted fractions", fracs
print "Errors", frac_errs
print "Scales" , scales

# Make pretty table
print "\\begin{center} \\Fontvi \\begin{tabular}{|l|| r| r| r|} \\hline"
print "Template & Initial Fraction & Fitted Fraction & Scale \\\\ \\hline"
if args.DD:
  print "%s & %.2f & %.2f $ \\pm $ %.2f & %.2f \\\\" % (latex_dict[fffp_template] + " (DD fake/fake+fake/prompt)", initial_fracs[0], fracs[0], frac_errs[0], scales[0])
  print "%s & %.2f & %.2f $ \\pm $ %.2f & %.2f \\\\ \\hline" % (latex_dict[pp_template] + " (prompt/prompt)", initial_fracs[1], fracs[1], frac_errs[1], scales[1])
else:
  print "%s & %.2f & %.2f $ \\pm $ %.2f & %.2f \\\\" % (latex_dict[ff_template] + " (fake/fake)", initial_fracs[0], fracs[0], frac_errs[0], scales[0])
  print "%s & %.2f & %.2f $ \\pm $ %.2f & %.2f \\\\" % (latex_dict[fp_template] + " (fake/prompt)", initial_fracs[1], fracs[1], frac_errs[1], scales[1])
  print "%s & %.2f & %.2f $ \\pm $ %.2f & %.2f \\\\ \\hline" % (latex_dict[pp_template] + " (prompt/prompt)", initial_fracs[2], fracs[2], frac_errs[2], scales[2])
print "\\end{tabular} \\end{center}"

