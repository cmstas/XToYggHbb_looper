import ROOT
from array import array
import argparse
from collections import OrderedDict
import copy
from datetime import date    
import glob
import numpy
import os
import plotUtils


sampleFillColor=dict()
sampleFillColor["Data"]              = None
sampleFillColor["QCD"]               = ROOT.kMagenta
sampleFillColor["DDQCDGJets"]        = ROOT.kMagenta
sampleFillColor["GJets"]             = ROOT.kCyan-7
sampleFillColor["Diphoton"]          = ROOT.kSpring+5
sampleFillColor["tt+X"]              = ROOT.kPink-4
sampleFillColor["DY"]                = ROOT.kOrange+3
sampleFillColor["VG"]                = ROOT.kMagenta+2
sampleFillColor["Diboson"]           = ROOT.kPink
sampleFillColor["Single H"]          = ROOT.kGreen+2
sampleFillColor["ttH_M125"]          = ROOT.kOrange+7
sampleFillColor["HHbbgg"]            = ROOT.kBlue-4
sampleFillColor["NMSSM_XToYHTo2G2B"] = None
sampleFillColor["NMSSM_XToYHTo2B2G"] = None

sampleLineColor=dict()
sampleLineColor["Data"]              = ROOT.kBlack
sampleLineColor["QCD"]               = None
sampleLineColor["DDQCDGJets"]        = None
sampleLineColor["GJets"]             = None
sampleLineColor["Diphoton"]          = None
sampleLineColor["tt+X"]              = None
sampleLineColor["DY"]                = None
sampleLineColor["VG"]                = None
sampleLineColor["Diboson"]           = None
sampleLineColor["Single H"]          = None
sampleLineColor["ttH_M125"]          = None
sampleLineColor["HHbbgg"]            = None
sampleLineColor["NMSSM_XToYHTo2G2B"] = ROOT.kRed
sampleLineColor["NMSSM_XToYHTo2B2G"] = ROOT.kRed

sampleLineWidth=dict()
sampleLineWidth["Data"]              = 1
sampleLineWidth["QCD"]               = 0
sampleLineWidth["DDQCDGJets"]        = 0
sampleLineWidth["GJets"]             = 0
sampleLineWidth["Diphoton"]          = 0
sampleLineWidth["tt+X"]              = 0
sampleLineWidth["DY"]                = 0
sampleLineWidth["VG"]                = 0
sampleLineWidth["Diboson"]           = 0
sampleLineWidth["Single H"]          = 0
sampleLineWidth["ttH_M125"]          = 0
sampleLineWidth["HHbbgg"]            = 0
sampleLineWidth["NMSSM_XToYHTo2G2B"] = 2
sampleLineWidth["NMSSM_XToYHTo2B2G"] = 2

sampleMarkerStyle=dict()
sampleMarkerStyle["Data"]              = 20
sampleMarkerStyle["QCD"]               = None
sampleMarkerStyle["DDQCDGJets"]        = None
sampleMarkerStyle["GJets"]             = None
sampleMarkerStyle["Diphoton"]          = None
sampleMarkerStyle["tt+X"]              = None
sampleMarkerStyle["DY"]                = None
sampleMarkerStyle["VG"]                = None
sampleMarkerStyle["Diboson"]           = None
sampleMarkerStyle["Single H"]          = None
sampleMarkerStyle["ttH_M125"]          = None
sampleMarkerStyle["HHbbgg"]            = None
sampleMarkerStyle["NMSSM_XToYHTo2G2B"] = None
sampleMarkerStyle["NMSSM_XToYHTo2B2G"] = None

sampleMarkerSize=dict()
sampleMarkerSize["Data"]              = 1.2
sampleMarkerSize["QCD"]               = None
sampleMarkerSize["DDQCDGJets"]        = None
sampleMarkerSize["GJets"]             = None
sampleMarkerSize["Diphoton"]          = None
sampleMarkerSize["tt+X"]              = None
sampleMarkerSize["DY"]                = None
sampleMarkerSize["VG"]                = None
sampleMarkerSize["Diboson"]           = None
sampleMarkerSize["Single H"]          = None
sampleMarkerSize["ttH_M125"]          = None
sampleMarkerSize["HHbbgg"]            = None
sampleMarkerSize["NMSSM_XToYHTo2G2B"] = None
sampleMarkerSize["NMSSM_XToYHTo2B2G"] = None

sampleLegend=dict()
sampleLegend["Data"]              = "Data"
sampleLegend["QCD"]               = "QCD"
sampleLegend["DDQCDGJets"]        = "QCD+GJets"
sampleLegend["GJets"]             = "GJets"
sampleLegend["Diphoton"]          = "Diphoton"
sampleLegend["tt+X"]              = "tZq,tt+X"
sampleLegend["DY"]                = "DY"
sampleLegend["VG"]                = "VG"
sampleLegend["Diboson"]           = "Diboson"
sampleLegend["Single H"]          = "Single H"
sampleLegend["ttH_M125"]          = "ttH"
sampleLegend["HHbbgg"]            = "ggHH"
sampleLegend["NMSSM_XToYHTo2G2B"] = "XYH"
sampleLegend["NMSSM_XToYHTo2B2G"] = "XYH (inv)"

epsilon = 1e-6
# Normalize to event counts (experimental - here for posterity)
#def DDQCDGJetsSF(tree):
#  tree.Draw("weight_central>>h","","goff")
#  h = ROOT.gDirectory.Get("h")
#  w = float(h.GetMean())
#  print("%f",w)
#  return float(1/w)

def BTagSF(tree):
  tree.Draw("weight_beforeBTagSF>>hBefore","","goff")
  hBefore = ROOT.gDirectory.Get("hBefore")
  wBefore = hBefore.GetEntries()*hBefore.GetMean()
  tree.Draw("weight_afterBTagSF>>hAfter","","goff")
  hAfter = ROOT.gDirectory.Get("hAfter")
  wAfter = hAfter.GetEntries()*hAfter.GetMean()
  SF = wBefore / wAfter
  return SF

def get_plots(samples, year, plotname, cut, plotBin, plotXTitle, plotYTitle):

  htempDict=OrderedDict()
  htemp = None
  nBinsX = 0; lowBinX = 0.0; highBinX = 0.0; binsX = None; arrayX = None
  nBinsY = 0; lowBinY = 0.0; highBinY = 0.0; binsY = None; arrayY = None

  if isinstance(plotBin, list): # Fixed binning
    nBinsX = plotBin[0]; lowBinX = plotBin[1]; highBinX = plotBin[2]
    if ":" in plotname: # 2D histo
      nBinsY = plotBin[3]; lowBinY = plotBin[4]; highBinY = plotBin[5]
  else: # Variable binning
    binsX = plotBin.lstrip("{").rstrip("}").split(";")[0]; nBinsX = binsX.count(",")
    arrayX = array("f",[ float(x) for x in binsX.split(",")])
    if ":" in plotname: # 2D histo
      binsY = plotBin.lstrip("{").rstrip("}").split(";")[1]; nBinsY = binsY.count(",")
      arrayY = array("f",[ float(x) for x in binsY.split(",")])

  if year!="all" and year!="2016":
    years=[year]
  elif year=="2016":
    years=["2016nonAPV","2016APV"]
  else:
    years=["2016nonAPV","2016APV","2017","2018"]
  for tyear in years:
    # Place to apply cut on the b tagging score -- 2016APV vs. 2016nonAPV separation not included in the tree branches.
    #if tyear=="2018":
    #  cut = "weight_central*( dijet_lead_btagDeepFlavB>0.2783 && dijet_sublead_btagDeepFlavB>0.2783 )"
    #if tyear=="2017":
    #  cut = "weight_central*( dijet_lead_btagDeepFlavB>0.3040 && dijet_sublead_btagDeepFlavB>0.3040 )"
    #if tyear=="2016nonAPV":
    #  cut = "weight_central*( dijet_lead_btagDeepFlavB>0.2489 && dijet_sublead_btagDeepFlavB>0.2489 )"
    #if tyear=="2016APV":
    #  cut = "weight_central*( dijet_lead_btagDeepFlavB>0.2598 && dijet_sublead_btagDeepFlavB>0.2598 )"
    for i,sample in enumerate(samples):
      if "NMSSM_XToYH" in sample:
        for mass in args.signalMass:
          infile = ROOT.TFile(args.inDir+"output_"+sample+"_"+mass+"_"+tyear+".root")
          tree = infile.Get("tout")
          if tree.GetEntries() == 0:
            print("0 entries for sample %s%s, skipping..."%(sample,mass))
            continue
          if arrayX: # Variable binning
            if arrayY: # 2D histo
              htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX,nBinsY,arrayY);
            else: # 1D histo
              htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX);
          else: # Fixed binning
            if nBinsY > 0: # 2D histo
              htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX,nBinsY,lowBinY,highBinY);
            else: # 1D histo
              htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX);
          tree.Draw(plotname+">>htemp",cut,"goff")
          if (sample+"_"+mass) not in htempDict.keys():
            htempDict[sample+"_"+mass]=[]
          htemp = ROOT.gDirectory.Get("htemp")
          htempDict[sample+"_"+mass].append(copy.deepcopy(htemp))
          htempDict[sample+"_"+mass][-1].Scale(BTagSF(tree))
      elif sample=="GJets":
        for m1,m2 in zip(["40","100","200","400","600"],["100","200","400","600","Inf"]): 
          infile = ROOT.TFile(args.inDir+"output_GJets_HT-"+m1+"To"+m2+"_"+tyear+".root")
          tree = infile.Get("tout")
          if arrayX: # Variable binning
            if arrayY: # 2D histo
              htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX,nBinsY,arrayY);
            else: # 1D histo
              htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX);
          else: # Fixed binning
            if nBinsY > 0: # 2D histo
              htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX,nBinsY,lowBinY,highBinY);
            else: # 1D histo
              htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX);
          tree.Draw(plotname+">>htemp",cut,"goff")
          if sample not in htempDict.keys():
            htempDict[sample]=[]
          htemp = ROOT.gDirectory.Get("htemp")
          htempDict[sample].append(copy.deepcopy(htemp))
          htempDict[sample][-1].Scale(BTagSF(tree))
      else:
        infile = ROOT.TFile(args.inDir+"output_"+sample+"_"+tyear+".root")
        tree = infile.Get("tout")
        if arrayX: # Variable binning
          if arrayY: # 2D histo
            htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX,nBinsY,arrayY);
          else: # 1D histo
            htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,arrayX);
        else: # Fixed binning
          if nBinsY > 0: # 2D histo
            htemp = ROOT.TH2F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX,nBinsY,lowBinY,highBinY);
          else: # 1D histo
            htemp = ROOT.TH1F("htemp",";"+plotXTitle+";"+plotYTitle,nBinsX,lowBinX,highBinX);
        tree.Draw(plotname+">>htemp",cut,"goff")
        if sample not in htempDict.keys():
          htempDict[sample]=[]
        htemp = ROOT.gDirectory.Get("htemp")
        htempDict[sample].append(copy.deepcopy(htemp))
        if not (sample=="DDQCDGJets" or sample=="Data"):
          htempDict[sample][-1].Scale(BTagSF(tree))
        # Normalize to event counts (experimental - here for posterity)
        #if sample=="DDQCDGJets": htempDict[sample][-1].Scale(DDQCDGJetsSF(tree))
        # Process scaling based on 2D template fit
        ## Using data events and separately scaling low and high mass DiPhoton samples
        #if sample=="DDQCDGJets": htempDict[sample][-1].Scale(0.16)
        #if sample=="DiPhotonLow": htempDict[sample][-1].Scale(6.66)
        #if sample=="DiPhoton": htempDict[sample][-1].Scale(0.87)

  plotDict=OrderedDict()
  groupedSamples = OrderedDict()
  tempGroups = OrderedDict()
  tempGroups["QCD"] = ["QCD_Pt-30to40_MGG-80toInf","QCD_Pt-30toInf_MGG-40to80","QCD_Pt-40toInf_MGG-80toInf"]
  tempGroups["Diphoton"] = ["DiPhotonLow","DiPhoton"]
  tempGroups["tt+X"] = ["TTGG","TTGJets","TTJets","tZq","TTW","TTZ"]
  tempGroups["VG"] = ["WG","ZG"]
  tempGroups["Diboson"]   = ["WW","WZ","ZZ"]
  tempGroups["Single H"]   = ["ggHToDiPhoM125","VBFH_M125","VH_M125"]
  for sample in htempDict.keys():
    if sample in tempGroups["QCD"]:
      if "QCD" not in groupedSamples.keys():
        groupedSamples["QCD"]=[]
      groupedSamples["QCD"].append(sample)
    elif sample in tempGroups["Diphoton"]:
      if "Diphoton" not in groupedSamples.keys():
        groupedSamples["Diphoton"]=[]
      groupedSamples["Diphoton"].append(sample)
    elif sample in tempGroups["tt+X"]:
      if "tt+X" not in groupedSamples.keys():
        groupedSamples["tt+X"]=[]
      groupedSamples["tt+X"].append(sample)
    elif sample in tempGroups["VG"]:
      if "VG" not in groupedSamples.keys():
        groupedSamples["VG"]=[]
      groupedSamples["VG"].append(sample)
    elif sample in tempGroups["Diboson"]:
      if "Diboson" not in groupedSamples.keys():
        groupedSamples["Diboson"]=[]
      groupedSamples["Diboson"].append(sample)
    elif sample in tempGroups["Single H"]:
      if "Single H" not in groupedSamples.keys():
        groupedSamples["Single H"]=[]
      groupedSamples["Single H"].append(sample)
    else:
      groupedSamples[sample] = [sample]

  for gsample in groupedSamples.keys():
    tplot=None
    for sample in groupedSamples[gsample]:
      for tsample in htempDict.keys():
        if not tsample==sample:
          continue
        for htemp in htempDict[tsample]:
          if not tplot:
            tplot = copy.deepcopy(htemp)
          else:
            tplot.Add(htemp)

    #for b in range(0, tplot.GetNbinsX()+2):
    #  if tplot.GetBinContent(b)<0.0 or numpy.isnan(tplot.GetBinContent(b)) or not numpy.isfinite(tplot.GetBinContent(b)):
    #    tplot.SetBinContent(b,0.0)
    #    tplot.SetBinError(b,0.0)

    plotDict[gsample] = tplot

  # Process scaling based on 2D template fit
  ## Not using data events
  #plotDict["QCD"].Scale(0.305)
  #plotDict["GJets"].Scale(0.951)
  #plotDict["Diphoton"].Scale(0.001)
  ## Using data events
  #plotDict["DDQCDGJets"].Scale(0.9)
  #plotDict["Diphoton"].Scale(1.47)

  return plotDict

def customize_plot(sample, plot, fillColor, lineColor, lineWidth, markerStyle, markerSize):

  if plot.ClassName() == "TH1F": # Remove overflow/underflow for 2D histos
    error = ROOT.TMath.Sqrt(plot.GetBinError(0)*plot.GetBinError(0)+plot.GetBinError(1)*plot.GetBinError(1))
    plot.SetBinContent(1, plot.GetBinContent(1) + plot.GetBinContent(0))
    plot.SetBinError(1, error)
    plot.SetBinContent(0, 0.0)
    plot.SetBinError(0, 0.0)

    error = ROOT.TMath.Sqrt(plot.GetBinError(plot.GetNbinsX()+1)*plot.GetBinError(plot.GetNbinsX()+1)+plot.GetBinError(plot.GetNbinsX())*plot.GetBinError(plot.GetNbinsX()))
    plot.SetBinContent(plot.GetNbinsX(), plot.GetBinContent(plot.GetNbinsX()+1) + plot.GetBinContent(plot.GetNbinsX()))
    plot.SetBinError(plot.GetNbinsX(), error)
    plot.SetBinContent(plot.GetNbinsX()+1, 0.0)
    plot.SetBinError(plot.GetNbinsX()+1, 0.0)

  if fillColor: 
    plot.SetFillColor(fillColor)
    plot.SetLineColor(fillColor)
    plot.SetMarkerColor(fillColor)
  if lineColor: 
    plot.SetLineColor(lineColor)
    plot.SetMarkerColor(lineColor)
  if lineWidth:
    plot.SetLineWidth(lineWidth)
  if markerStyle:
    plot.SetMarkerStyle(markerStyle)
  if markerSize:
    plot.SetMarkerSize(markerSize)
  #plot.Sumw2()

  ### Remove spikes
  if sample!="Data" and not "met_pt" in plot.GetName():
    for b in range(1, plot.GetNbinsX()+1):
      if plot.GetBinContent(b)>0 and plot.GetBinError(b)/plot.GetBinContent(b)>0.75:
        plot.SetBinContent(b,0.0)
        plot.SetBinError(b,0.0)

  return plot

def get_yields(plotDict, plotname, lumi, year, plotData=False):
  curYields=OrderedDict()
  curErrors=OrderedDict()
  totalSMYield = 0.0
  totalSMError = 0.0

  for i,sample in enumerate(plotDict.keys()):
    # Signal
    if "NMSSM_XToYH" in sample:
      curYields[sample] = plotDict[sample].GetBinContent(1)
      curErrors[sample] = plotDict[sample].GetBinError(1)
      print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
    # Data
    elif sample=="Data": 
      if plotData:
        curYields[sample] = plotDict[sample].GetBinContent(1)
        curErrors[sample] = plotDict[sample].GetBinError(1)
        print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
    # Bkg
    else:
      curYields[sample] = plotDict[sample].GetBinContent(1)
      curErrors[sample] = plotDict[sample].GetBinError(1)
      print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
      if not totalSMYield:
        totalSMYield = curYields[sample]
        totalSMError = curErrors[sample]*curErrors[sample]
      else:
        totalSMYield = totalSMYield + curYields[sample]
        totalSMError = totalSMError + curErrors[sample]
  totalSMError = ROOT.TMath.Sqrt(totalSMError)
  print("Total_SM:\t%.2f +/- %.2f"%(totalSMYield,totalSMError))

def draw_plot(plotDict, plotname, lumi, year, logY=True, logX=False, plotData=False, doRatio=True):
  f = None
  if args.outFormat == "root": f = ROOT.TFile(args.outDir+plotname+args.label+".root","RECREATE")

  # Labels
  latex = ROOT.TLatex()
  latex.SetTextFont(42)
  latex.SetTextAlign(31)
  latex.SetTextSize(0.04)
  latex.SetNDC(True)

  latexCMS = ROOT.TLatex()
  latexCMS.SetTextFont(61)
  latexCMS.SetTextSize(0.04)
  latexCMS.SetNDC(True)

  latexCMSExtra = ROOT.TLatex()
  latexCMSExtra.SetTextFont(52)
  latexCMSExtra.SetTextSize(0.04)
  latexCMSExtra.SetNDC(True)

  legoffset = 0.0
  latexSel = ROOT. TLatex()
  latexSel.SetTextAlign(11)
  latexSel.SetTextFont(42)
  latexSel.SetTextSize(0.03)
  latexSel.SetNDC(True)

  yearenergy=""
  if year!="all" or lumi<100.0:
    if year!="all":
      yearenergy="%.1f fb^{-1} (%s, 13 TeV)"%(lumi,year)
    else:
      yearenergy="%.1f fb^{-1} (2016-2018, 13 TeV)"%(lumi)
  else:
    yearenergy="%.0f fb^{-1} (13 TeV)"%(lumi)
  if plotData:
    cmsExtra="Preliminary"
  else:
    cmsExtra="Simulation"

  curPlots=OrderedDict()

  totalSM = None
  lowToHighBinsCumulative = True
  for i,sample in enumerate(plotDict.keys()):
    # Signal
    if "NMSSM_XToYH" in sample:
      model = sample.split("_M")[0]
      mass = sample.split("B_")[1] if "B_" in sample else sample.split("G_")[1]
      curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[model],sampleLineColor[model]+i%len(args.signalMass),sampleLineWidth[model],sampleMarkerStyle[model],sampleMarkerSize[model]))
      if args.shape and curPlots[sample].Integral(0,-1)>0.0:
        curPlots[sample].Scale(1.0/curPlots[sample].Integral(0,-1))
      if args.cumulative:
        curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
    # Data
    elif sample=="Data": 
      if plotData:
        curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[sample],sampleLineColor[sample],sampleLineWidth[sample],sampleMarkerStyle[sample],sampleMarkerSize[sample]))
        if args.shape and curPlots[sample].Integral(0,-1)>0.0:
          curPlots[sample].Scale(1.0/curPlots[sample].Integral(0,-1))
        if args.cumulative:
          curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
    # Bkg
    else:
      curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[sample],sampleLineColor[sample],sampleLineWidth[sample],sampleMarkerStyle[sample],sampleMarkerSize[sample]))
      if not totalSM:
        totalSM = curPlots[sample].Clone("totalSM")
      else:
        totalSM.Add(curPlots[sample])

  if args.dataOnly:
    totalSM = curPlots["Data"].Clone("totalSM")
  is1Dplot = totalSM.ClassName() == "TH1F"
  totalScale = totalSM.Integral(0,-1) if is1Dplot else totalSM.Integral(0,-1,0,-1)
  if args.cumulative:
    totalSM = plotUtils.GetCumulative(totalSM,lowToHighBinsCumulative)
  if args.shape and totalScale>0.0:
    totalSM.Scale(1.0/totalScale)

  # Build stack
  stack = ROOT.THStack("stack","")
  if not args.dataOnly:
    for i,sample in enumerate(reversed(plotDict.keys())):
      # Bkg
      if not ("NMSSM_XToYH" in sample or sample=="Data"):
        if args.shape and totalScale>0.0:
          curPlots[sample].Scale(1.0/totalScale)
        if args.cumulative:
          curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
        stack.Add(curPlots[sample])


  # Plot legends, ranges
  legendXOffsetNoSelPrint = 0.18
  legendYOffsetNoSelPrint = 0.1
  if args.data:
    legend = ROOT.TLegend(0.5,0.7,0.91,0.91)
    if args.dataOnly:
      legend = ROOT.TLegend(0.5,0.8,0.89,0.89)
  else:
    legend = ROOT.TLegend(0.5,0.7,0.89,0.89)
  if not is1Dplot:
    legend = ROOT.TLegend(0.5,0.7,0.88,0.88)
  legend.SetLineColor(0)
  legend.SetLineWidth(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetNColumns(2)
  
  for sample in curPlots.keys():
    # Signal
    if "NMSSM_XToYH" in sample and not args.dataOnly:
      model = sample.split("_M")[0]
      mass = sample.split("B_")[1] if "B_" in sample else sample.split("G_")[1]
      massX = mass.split("_")[1]
      massY = mass.split("_")[3]
      legend.AddEntry(curPlots[sample],sampleLegend[model]+" ("+massX+"/"+massY+") GeV","L")
    # Data
    elif sample=="Data": 
      if plotData:
        legend.AddEntry(curPlots[sample],sampleLegend[sample],"EPL")
    # Bkg
    else:
      if not args.dataOnly and not args.signalOnly:
        legend.AddEntry(curPlots[sample], sampleLegend[sample],"F")
    

  # Define canvas
  canvas = ROOT.TCanvas("canvas","canvas",800,800)

  MCplot = copy.deepcopy(totalSM)
  g_unc = ROOT.TGraphAsymmErrors()
  g_data = ROOT.TGraphAsymmErrors()
  g_data_clone = ROOT.TGraphAsymmErrors()
  g_ratio = ROOT.TGraphAsymmErrors()
  g_ratio_unc = ROOT.TGraphAsymmErrors()
  g_ratio_signal = ROOT.TMultiGraph()

  h_axis = ROOT.TH1F() if is1Dplot else \
           ROOT.TH2F()
  h_axis_ratio = ROOT.TH1F()
  h_axis = ROOT.TH1F("h_axis","", MCplot.GetNbinsX(), MCplot.GetXaxis().GetBinLowEdge(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX())) if is1Dplot else \
           ROOT.TH2F("h_axis","", MCplot.GetNbinsX(), MCplot.GetXaxis().GetBinLowEdge(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()),
                                  MCplot.GetNbinsY(), MCplot.GetYaxis().GetBinLowEdge(1), MCplot.GetYaxis().GetBinUpEdge(MCplot.GetNbinsY()))
  h_axis_ratio = ROOT.TH1F("h_axis_ratio","", MCplot.GetNbinsX(), MCplot.GetXaxis().GetBinLowEdge(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
  if logX and MCplot.GetXaxis().GetBinLowEdge(1) < epsilon:
    h_axis.GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
    h_axis_ratio.GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))

  doRatio=False
  if is1Dplot:
    if plotData:
      if not args.dataOnly:
        doRatio=True

      #plotUtils.ConvertToPoissonGraph(curPlots["Data"], g_data, drawZeros=True, drawXerr=False)
      plotUtils.ConvertToPoissonGraph(curPlots["Data"], g_data, drawZeros=False, drawXerr=False)
      g_data.SetMarkerStyle(20)
      g_data.SetMarkerSize(1.2)
      g_data.SetLineWidth(1)
      # draw with zero marker size so error bars drawn all the way to x axis in the case of 0 content
      g_data_clone = g_data.Clone()
      g_data_clone.SetMarkerSize(0.0)

      #plotUtils.GetPoissonRatioGraph(MCplot, curPlots["Data"], g_ratio, drawZeros=True, drawXerr=False, useMCErr=False)
      plotUtils.GetPoissonRatioGraph(MCplot, curPlots["Data"], g_ratio, drawZeros=False, drawXerr=False, useMCErr=False)
      g_ratio.SetMarkerStyle(20)
      g_ratio.SetMarkerSize(1.2)
      g_ratio.SetLineWidth(1)

    if not plotData and args.shape and doSignalMCRatio:
      doRatio=True

      for i,sample in enumerate(plotDict.keys()):
        # Signal
        if "NMSSM_XToYH" in sample:
          model = sample.split("_M")[0] 
          mass = sample.split("B_")[1] if "B_" in sample else sample.split("G_")[1]
          g_signal_temp = ROOT.TGraphAsymmErrors()
          plotUtils.ConvertToPoissonGraph(curPlots[sample], g_signal_temp, drawZeros=False, drawXerr=False, drawYerr=False)
          g_signal_temp.SetMarkerStyle(20)
          g_signal_temp.SetMarkerSize(1.2)
          g_signal_temp.SetLineWidth(1)

          # draw with zero marker size so error bars drawn all the way to x axis in the case of 0 content
          g_signal_temp_clone = g_signal_temp.Clone()
          g_signal_temp_clone.SetMarkerSize(0.0)

          g_ratio_signal_temp = ROOT.TGraphAsymmErrors()
          plotUtils.GetPoissonRatioGraph(MCplot, curPlots[sample], g_ratio_signal_temp, drawZeros=False, drawXerr=False, drawYerr=False, useMCErr=False)
          g_ratio_signal_temp.SetMarkerStyle(20)
          g_ratio_signal_temp.SetMarkerSize(1.2)
          g_ratio_signal_temp.SetMarkerColor(sampleLineColor[model]+i%len(args.signalMass))
          g_ratio_signal_temp.SetLineWidth(1)
          g_ratio_signal.Add(copy.deepcopy(g_ratio_signal_temp))

    for b in range(1,MCplot.GetNbinsX()+1):
      thisPoint = g_ratio_unc.GetN()
      yerror = MCplot.GetBinError(b)
      g_unc.SetPoint(thisPoint, MCplot.GetBinCenter(b), MCplot.GetBinContent(b))
      g_unc.SetPointError(thisPoint, 0.5*MCplot.GetBinWidth(b), 0.5*MCplot.GetBinWidth(b), yerror, yerror)
      if MCplot.GetBinContent(b)>0.0:
        yerror = yerror/MCplot.GetBinContent(b)
      else:
        yerror = 0.0
      g_ratio_unc.SetPoint(thisPoint, MCplot.GetBinCenter(b), 1.0)
      g_ratio_unc.SetPointError(thisPoint, 0.5*MCplot.GetBinWidth(b), 0.5*MCplot.GetBinWidth(b), yerror, yerror)
    g_unc.SetFillStyle(3244)
    g_unc.SetFillColor(ROOT.kGray+3)
    g_ratio_unc.SetFillStyle(1001)
    g_ratio_unc.SetFillColor(ROOT.kGray)

  pads = []
  if doRatio==True:
    minR=0.0
    maxR=2.0
    ty = numpy.array([])
    if not doSignalMCRatio:
      tmax=maxR
      if args.data:
        ty = g_ratio.GetY()
      else:
        ty = g_ratio_signal.GetY()
      if len(ty)>0:
        tmax = numpy.amax(ty)
      if tmax>maxR:
        maxR=tmax*1.05
      if maxR>5.0:
        minR=0.1
    h_axis_ratio.GetYaxis().SetRangeUser(minR,maxR)
    h_axis_ratio.SetMinimum(minR)
    h_axis_ratio.SetMaximum(maxR)
    h_axis_ratio.SetTitle(";;Data / MC")
    h_axis_ratio.GetYaxis().SetTitleSize(0.16)
    h_axis_ratio.GetYaxis().SetTitleOffset(0.25)
    if logY:
      h_axis_ratio.GetYaxis().SetTitleOffset(0.3)
    h_axis_ratio.GetYaxis().SetLabelSize(0.12)
    h_axis_ratio.GetYaxis().CenterTitle()
    h_axis_ratio.GetYaxis().SetTickLength(0.02)
    h_axis_ratio.GetXaxis().SetLabelSize(0)
    h_axis_ratio.GetXaxis().SetTitle("")
    h_axis_ratio.GetXaxis().SetTickSize(0.06)

    line = ROOT.TLine(h_axis.GetXaxis().GetBinLowEdge(1), 1.0, h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()), 1.0)

    pads.append(ROOT.TPad("1","1",0.0,0.18,1.0,1.0))
    pads.append(ROOT.TPad("2","2",0.0,0.0,1.0,0.19))
    pads[0].SetTopMargin(0.08)
    pads[0].SetBottomMargin(0.13)
    pads[0].SetRightMargin(0.05)
    pads[0].SetLeftMargin(0.10)
    pads[1].SetRightMargin(0.05)
    pads[1].SetLeftMargin(0.10)
    pads[0].Draw()
    pads[1].Draw()
    pads[1].cd()
    if maxR>5.0:
      pads[1].SetLogy()
    pads[1].SetTickx()
    if logX:
      h_axis_ratio.GetXaxis().SetMoreLogLabels()
      pads[1].SetLogx()
    if plotData:
      h_axis_ratio.Draw("")
      g_ratio_unc.Draw("SAME,2")
      g_ratio.Draw("SAME,P0")
    else:
      g_ratio_signal.Draw("SAME,P0")
      g_ratio_signal.GetXaxis().SetLimits(h_axis.GetXaxis().GetBinLowEdge(1),h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()));
      g_ratio_signal.GetHistogram().GetXaxis().SetRangeUser(h_axis.GetXaxis().GetBinLowEdge(1),h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()));
      g_ratio_signal.GetHistogram().GetYaxis().SetRangeUser(0.,2.0);

      g_ratio_signal.GetHistogram().SetTitle(";;Signal / MC")
      g_ratio_signal.GetHistogram().GetYaxis().SetTitleSize(0.16)
      g_ratio_signal.GetHistogram().GetYaxis().SetTitleOffset(0.25)
      g_ratio_signal.GetHistogram().GetYaxis().SetLabelSize(0.12)
      g_ratio_signal.GetHistogram().GetYaxis().CenterTitle()
      g_ratio_signal.GetHistogram().GetYaxis().SetTickLength(0.02)

      g_ratio_signal.GetHistogram().GetXaxis().SetLabelSize(0)
      g_ratio_signal.GetHistogram().GetXaxis().SetTitle("")
      g_ratio_signal.GetHistogram().GetXaxis().SetTickSize(0.06)
      if logX:
        if MCplot.GetXaxis().GetBinLowEdge(1) < epsilon:
          g_ratio_signal.GetHistogram().GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
        g_ratio_signal.GetHistogram().GetXaxis().SetMoreLogLabels()
        pads[1].SetLogx()

    #
    line.SetLineStyle(2)
    line.SetLineColor(sampleLineColor["Data"])
    line.SetLineWidth(1)
    line.Draw("SAME")
    #
    #pads[1].RedrawAxis()
    pads[1].Modified();
    pads[1].Update();

  else:
    pads.append(ROOT.TPad("1","1",0,0,1,1))
    if not is1Dplot:
      pads[0].SetRightMargin(0.15)
    pads[0].Draw()

  pads[0].cd()
  if logY:
    pads[0].SetLogy()
  if logX:
    h_axis.GetXaxis().SetMoreLogLabels()
    pads[0].SetLogx()


  #plot data, stack, signal, data  
  h_axis.GetYaxis().SetTitleSize(0.04)
  h_axis.GetXaxis().SetTitleSize(0.04)
  h_axis.GetXaxis().SetTitleOffset(1.25)
  h_axis.GetXaxis().SetTitle(MCplot.GetXaxis().GetTitle())
  #
  h_axis.GetYaxis().SetLabelSize(0.03 if is1Dplot else 0.0)
  if args.dataOnly:
    h_axis.GetYaxis().SetTitleOffset(1.35)
  if args.shape:
    h_axis.GetYaxis().SetTitle("A.U.")
  else:
    h_axis.GetYaxis().SetTitle(MCplot.GetYaxis().GetTitle())
  if not args.shape:
    h_axis.GetYaxis().SetMaxDigits(3)
  #
  h_axis.Draw("" if is1Dplot else "COLZ")
  #
  if not args.dataOnly and not args.signalOnly:
    if is1Dplot:
      stack.Draw("HIST,SAME")
      g_unc.Draw("SAME,2")
    else:
      MCplot.Draw("COLZ")
      MCplot.GetYaxis().SetTitleOffset(1.40)
      MCplot.GetZaxis().SetLabelSize(0.03)
  histMax = 0.0
  if plotData and not args.signalOnly:
    if histMax < curPlots["Data"].GetMaximum():
      histMax = curPlots["Data"].GetMaximum()
    if is1Dplot:
      g_data.Draw("P,SAME")
      g_data_clone.Draw("P,SAME")
    else:
      curPlots["Data"].Draw("COLZ")
      curPlots["Data"].GetYaxis().SetTitleOffset(1.45)
      curPlots["Data"].GetZaxis().SetLabelSize(0.03)
  for sample in curPlots.keys():
    if "NMSSM_XToYH" in sample and not args.dataOnly:
      if histMax < curPlots[sample].GetMaximum(): 
        histMax = curPlots[sample].GetMaximum()
      curPlots[sample].Draw("HIST,SAME" if is1Dplot else "COLZ")
      if not is1Dplot:
        curPlots[sample].GetYaxis().SetTitleOffset(1.45)
        curPlots[sample].GetZaxis().SetLabelSize(0.03)

  if histMax < MCplot.GetMaximum(): 
    histMax = MCplot.GetMaximum()
  if logY:
    histMax = histMax*1e3
    h_axis.SetMinimum(1e-3)
  h_axis.SetMaximum(1.1*histMax)

  legend.Draw()
  pads[0].Update()
  pads[0].RedrawAxis()


  # Draw CMS headers
  expoffset=0.0
  if logY or 1.1*histMax<1000.0:
    expoffset=0
  if doRatio:
    latex.DrawLatex(0.95, 0.93+expoffset, yearenergy);
    latexCMS.DrawLatex(0.13,0.88+expoffset,"CMS");
    latexCMSExtra.DrawLatex(0.13,0.835+expoffset, cmsExtra);
  else:
    latex.DrawLatex(0.90, 0.91+expoffset, yearenergy);
    latexCMS.DrawLatex(0.13,0.86+expoffset,"CMS");
    latexCMSExtra.DrawLatex(0.13,0.815+expoffset, cmsExtra);

  #xStart = h_axis.GetXaxis().GetBinLowEdge(1)
  #xEnd = h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX())
  #xDiff = xEnd - xStart
  #yStart = h_axis.GetYaxis().GetBinLowEdge(1)
  #yEnd = h_axis.GetYaxis().GetBinUpEdge(h_axis.GetNbinsY())
  #yDiff = yEnd - yStart
  #pt = ROOT.TPaveText(xStart+0.04*xDiff,yEnd-0.25*yDiff,xStart+0.2*xDiff,yEnd-0.2*yDiff)
  #pt.AddText("r=%.3f"%totalSM.GetCorrelationFactor())
  #pt.Draw()


  # Print and save
  extension = "_"+year
  if plotData:
    if args.dataOnly:
      extension = extension+"_data"
    else:
      extension = extension+"_mcdata"
  elif args.signalOnly:
    extension = extension+"_signal"
  else:
    extension = extension+"_mc"
  if logX:
    extension = extension+"_logX"
  if logY:
    extension = extension+"_logY"
  if args.shape:
    extension = extension+"_areaNormalized"
  if args.cumulative:
    extension = extension+"_cumulative"
  
  canvas.SaveAs(args.outDir + plotname.replace("/","Over") + args.label + extension + "." + ("pdf" if args.outFormat=="pdf" else "png"))

  if args.outFormat == "root": f.Write()

if __name__=="__main__":
  ROOT.gStyle.SetOptStat(0)
  ROOT.gROOT.SetBatch(1)

  user = os.environ.get("USER")
  today = date.today().strftime("%b-%d-%Y")

  parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--inDir", default="./cpp/temp_data/", help="Choose input directory. Default: './cpp/temp_data/'")
  parser.add_argument("--outDir", default="/home/users/"+user+"/public_html/XToYHToggbb/plots_"+today, help="Choose output directory. Default: '/home/users/"+user+"/public_html/XToYHToggbb/plots_"+today+"'")
  parser.add_argument("--data", default=False, action="store_true", help="Plot data")
  parser.add_argument("--dataOnly", default=False, action="store_true", help="Plot only data, no MC bkg")
  parser.add_argument("--signalOnly", default=False, action="store_true", help="Plot only signal, no MC bkg")
  parser.add_argument("--noSignal", default=False, action="store_true", help="Do not plot signals")
  parser.add_argument("--signalMass", default=[], nargs="+", help="Signal mass points to plot. 'all' plots/prints all of the mass points. Default: 'MX_700_MY_100'")
  parser.add_argument("--invertedSignal", default=False, action="store_true", help="Plot Y->bb and H->gg, instead of Y->gg and H->bb.")
  parser.add_argument("--samples", default=[], nargs="+", help="Samples to plot. Default: 'all'")
  parser.add_argument("--shape", default=False, action="store_true", help="Shape normalization")
  parser.add_argument("--cumulative", default=False, action="store_true", help="Cumulative distributions")
  parser.add_argument("--years", default=[], nargs="+", help="List of years to be plotted. Default: all years")
  parser.add_argument("--noLogY", default=False, action="store_true", help="Don't use logY.")
  parser.add_argument("--logX", default=False, action="store_true", help="Use logX.")
  parser.add_argument("--unweighted", default=False, action="store_true", help="Plot unweighted histograms.")
  parser.add_argument("--cut", default="1", help="Selection to apply. Default: No selection on top of the preselection")
  parser.add_argument("--outFormat", default="png", choices=["png","pdf","root"], help="Output format: png, pdf or root. Default: png")
  parser.add_argument("--label", default="", help="Extra label for plots.")
  parser.add_argument("--yields", default=False, action="store_true", help="Print yields instead of plotting")
  parser.add_argument("--DDHistos", default=False, action="store_true", help="Set up plotter for histograms for the DD QCD+GJets estimation")
  parser.add_argument("--plotsToInclude", default=[], nargs="+", help="List of which plots are to be produced. An empty list prints all of the hardcoded ones. Default: Empty list")
  parser.add_argument("--plotsToExclude", default=[], nargs="+", help="List of which plots are NOT to be produced. An empty list excludes only the hardcoded ones. Default: Empty list")
  args = parser.parse_args()

  args.inDir = args.inDir.rstrip("/")+"/"
  args.outDir = args.outDir.rstrip("/")+"/"

  if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)
  os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+args.outDir)

  if args.dataOnly:
    args.data = True

  # Final state: Y->gg & H->bb -- Default
  signalSample = "NMSSM_XToYHTo2G2B"
  signalSampleSplit = "B_"
  # Final state: Y->bb & H->gg -- For cross-checks
  if args.invertedSignal:
    signalSample = "NMSSM_XToYHTo2B2G"
    signalSampleSplit = "G_"

  if len(args.signalMass)==0: 
    args.signalMass = ["MX_700_MY_100"]
  if args.signalMass==["all"]: 
    args.signalMass = []
    fileNames = glob.glob(args.inDir+"*"+signalSample+"*")
    for fileName in fileNames:
      fileName = fileName.split("/")[-1].split(".")[0].split(signalSampleSplit)[1].split("_201")[0]
      args.signalMass.append(fileName)

  if len(args.years)==0:
    args.years = ["all"]

  # Do signal/MC ratio
  doSignalMCRatio = False


  # Samples
  samples=[]
  if args.data:
    samples.append("Data")
  # SM MC
  if not args.dataOnly:
    samples.append("QCD_Pt-30to40_MGG-80toInf")
    samples.append("QCD_Pt-30toInf_MGG-40to80")
    samples.append("QCD_Pt-40toInf_MGG-80toInf")
    samples.append("GJets")
    samples.append("DDQCDGJets")
    samples.append("DiPhotonLow")
    samples.append("DiPhoton")
    samples.append("TTGG")
    samples.append("TTGJets")
    samples.append("TTJets")
    samples.append("tZq")
    samples.append("TTW")
    samples.append("TTZ")
    samples.append("DY")
    samples.append("WG")
    samples.append("ZG")
    samples.append("WW")
    samples.append("WZ")
    samples.append("ZZ")
    samples.append("VBFH_M125")
    samples.append("VH_M125")
    samples.append("ggHToDiPhoM125")
    samples.append("ttH_M125")
    samples.append("HHbbgg")
  if not args.DDHistos:
    samples.remove("QCD_Pt-30to40_MGG-80toInf")
    samples.remove("QCD_Pt-30toInf_MGG-40to80")
    samples.remove("QCD_Pt-40toInf_MGG-80toInf")
    samples.remove("GJets")
  # Signal MC
  if not args.noSignal:
    samples.append(signalSample)
  if len(args.samples) > 0:
    samples = args.samples


  for year in args.years:
    lumi=0.0 #fb^-1
    if year == "2018":
      lumi = 54.5
    elif year == "2017":
      lumi = 41.5
    elif year == "2016APV":
      lumi = 19.5
    elif year == "2016nonAPV":
      lumi = 16.8
    elif year == "2016":
      lumi = 19.5+16.8
    elif year == "all":
      lumi = 59.83 + 41.48 + 19.5 + 16.8


    # Cuts
    weight = "1" if args.unweighted else "weight_central"
    cut = args.cut # Default value if no cut is to be applied

    # Hardcoded selection of cuts commonly applied
    ## Separate BB. EB, EE photons
    #cut = "fabs(LeadPhoton_eta) < 1.442 && fabs(SubleadPhoton_eta) < 1.442"
    #cut = "( fabs(LeadPhoton_eta) < 1.442 && fabs(SubleadPhoton_eta) > 1.566 ) || ( fabs(LeadPhoton_eta) > 1.566 && fabs(SubleadPhoton_eta) < 1.442 )"
    #cut = "fabs(LeadPhoton_eta) > 1.566 && fabs(SubleadPhoton_eta) > 1.566"

    ## Proper DeepFlavor Loose WP cut
    #cut = "( year==2018 && dijet_lead_btagDeepFlavB>0.0490 && dijet_sublead_btagDeepFlavB>0.0490 ) || ( year==2017 && dijet_lead_btagDeepFlavB>0.0532 && dijet_sublead_btagDeepFlavB>0.0532 ) || ( year==2016nonAPV && dijet_lead_btagDeepFlavB>0.0480 && dijet_sublead_btagDeepFlavB>0.0480 ) || ( year==2016APV && dijet_lead_btagDeepFlavB>0.0508 && dijet_sublead_btagDeepFlavB>0.0508 )"
    ## Proper DeepFlavor Medium WP cut
    #cut = "( year==2018 && dijet_lead_btagDeepFlavB>0.2783 && dijet_sublead_btagDeepFlavB>0.2783 ) || ( year==2017 && dijet_lead_btagDeepFlavB>0.3040 && dijet_sublead_btagDeepFlavB>0.3040 ) || ( year==2016nonAPV && dijet_lead_btagDeepFlavB>0.2489 && dijet_sublead_btagDeepFlavB>0.2489 ) || ( year==2016APV && dijet_lead_btagDeepFlavB>0.2598 && dijet_sublead_btagDeepFlavB>0.2598 )"

    ## Diphoton mass cut
    #cut = "Diphoton_mass > 95"

    ## Photon pT / Mgg cuts
    #cut = "LeadPhoton_pt/Diphoton_mass > 0.33"
    #cut = "SubleadPhoton_pt/Diphoton_mass > 0.25"
    #cut = "LeadPhoton_pt/Diphoton_mass > 0.33 && SubleadPhoton_pt/Diphoton_mass > 0.25"

    ## Low MVA ID sideband
    #cut = "( ( LeadPhoton_mvaID >= SubleadPhoton_mvaID )*( SubleadPhoton_mvaID < (fabs(SubleadPhoton_eta)<1.442 ? -0.02 : -0.26) ) + ( SubleadPhoton_mvaID > LeadPhoton_mvaID )*( LeadPhoton_mvaID < (fabs(LeadPhoton_eta)<1.442 ? -0.02 : -0.26) ) )"
    ## Fake photons in low MVA ID sideband
    #cut = "( ( LeadPhoton_mvaID >= SubleadPhoton_mvaID )*( SubleadPhoton_mvaID < (fabs(SubleadPhoton_eta)<1.442 ? -0.02 : -0.26) && SubleadPhoton_genPartFlav==0 ) + ( SubleadPhoton_mvaID > LeadPhoton_mvaID )*( LeadPhoton_mvaID < (fabs(LeadPhoton_eta)<1.442 ? -0.02 : -0.26) && LeadPhoton_genPartFlav==0 ) )"
    ## Fake photons
    #cut = "( LeadPhoton_genPartFlav==0 || SubleadPhoton_genPartFlav==0 )"
    #cut = "( LeadPhoton_mvaID > -0.9 && SubleadPhoton_mvaID > -0.9 ) && LeadPhoton_genPartFlav==0 "
    #cut = "( LeadPhoton_mvaID > -0.9 && SubleadPhoton_mvaID > -0.9 ) && SubleadPhoton_genPartFlav==0 "

    cut = weight + "*(" + cut + ")"


    # List of plots
    plotNames = dict(); plotBins = dict(); plotXTitles = dict(); plotYTitles = dict()

    # 1D histos
    plotNames["xcand_pt"] = ""; plotBins["xcand_pt"] = [50,0,500]; plotXTitles["xcand_pt"] = "p_{T}(X) [GeV]"
    plotNames["xcand_eta"] = ""; plotBins["xcand_eta"] = [50,-3,3]; plotXTitles["xcand_eta"] = "#eta(X)"
    plotNames["xcand_phi"] = ""; plotBins["xcand_phi"] = [50,3.2,3.2]; plotXTitles["xcand_phi"] = "#phi(X)"
    plotNames["xcand_mass"] = ""; plotBins["xcand_mass"] = [50,200,1000]; plotXTitles["xcand_mass"] = "M(X) [GeV]"

    plotNames["LeadPhoton_pt"] = ""; plotBins["LeadPhoton_pt"] = [50,0,500]; plotXTitles["LeadPhoton_pt"] = "p_{T}(#gamma_{1}) [GeV]"
    plotNames["LeadPhoton_eta"] = ""; plotBins["LeadPhoton_eta"] = [50,-3,3]; plotXTitles["LeadPhoton_eta"] = "#eta(#gamma_{1})"
    plotNames["LeadPhoton_phi"] = ""; plotBins["LeadPhoton_phi"] = [50,-3.2,3.2]; plotXTitles["LeadPhoton_phi"] = "#phi(#gamma_{1})"
    plotNames["LeadPhoton_pixelSeed"] = ""; plotBins["LeadPhoton_pixelSeed"] = [2,0,2]; plotXTitles["LeadPhoton_pixelSeed"] = "hasPixelSeed(#gamma_{1})"
    plotNames["LeadPhoton_r9"] = ""; plotBins["LeadPhoton_r9"] = [50,0,2]; plotXTitles["LeadPhoton_r9"] = "R_{9}(#gamma_{1})"
    plotNames["LeadPhoton_sieie"] = ""; plotBins["LeadPhoton_sieie"] = [50,0,0.05]; plotXTitles["LeadPhoton_sieie"] = "#sigma_{ieie}(#gamma_{1})"
    plotNames["LeadPhoton_pfPhoIso03"] = ""; plotBins["LeadPhoton_pfPhoIso03"] = [50,0,20]; plotXTitles["LeadPhoton_pfPhoIso03"] = "PF Iso_{abs}^{#gamma}(#gamma_{1}) [GeV]"
    plotNames["LeadPhoton_chargedHadronIso"] = ""; plotBins["LeadPhoton_chargedHadronIso"] = [50,0,20]; plotXTitles["LeadPhoton_chargedHadronIso"] = "PF Iso_{abs}^{ch}(#gamma_{1}) [GeV]"
    plotNames["LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_mvaID"] = [20,-1,1]; plotXTitles["LeadPhoton_mvaID"] = "MVA ID(#gamma_{1})"
    plotNames["LeadPhoton_pt/Diphoton_mass"] = ""; plotBins["LeadPhoton_pt/Diphoton_mass"] = [50,0,2]; plotXTitles["LeadPhoton_pt/Diphoton_mass"] = "p_{T}(#gamma_{1}) / M(#gamma#gamma)"

    plotNames["SubleadPhoton_pt"] = ""; plotBins["SubleadPhoton_pt"] = [50,0,500]; plotXTitles["SubleadPhoton_pt"] = "p_{T}(#gamma_{2}) [GeV]"
    plotNames["SubleadPhoton_eta"] = ""; plotBins["SubleadPhoton_eta"] = [50,-3,3]; plotXTitles["SubleadPhoton_eta"] = "#eta(#gamma_{2})"
    plotNames["SubleadPhoton_phi"] = ""; plotBins["SubleadPhoton_phi"] = [50,-3.2,3.2]; plotXTitles["SubleadPhoton_phi"] = "#phi(#gamma_{2})"
    plotNames["SubleadPhoton_pixelSeed"] = ""; plotBins["SubleadPhoton_pixelSeed"] = [2,0,2]; plotXTitles["SubleadPhoton_pixelSeed"] = "hasPixelSeed(#gamma_{2})"
    plotNames["SubleadPhoton_r9"] = ""; plotBins["SubleadPhoton_r9"] = [50,0,2]; plotXTitles["SubleadPhoton_r9"] = "R_{9}(#gamma_{2})"
    plotNames["SubleadPhoton_sieie"] = ""; plotBins["SubleadPhoton_sieie"] = [50,0,0.05]; plotXTitles["SubleadPhoton_sieie"] = "#sigma_{ieie}(#gamma_{2})"
    plotNames["SubleadPhoton_pfPhoIso03"] = ""; plotBins["SubleadPhoton_pfPhoIso03"] = [50,0,20]; plotXTitles["SubleadPhoton_pfPhoIso03"] = "PF Iso_{abs}^{#gamma}(#gamma_{2}) [GeV]"
    plotNames["SubleadPhoton_chargedHadronIso"] = ""; plotBins["SubleadPhoton_chargedHadronIso"] = [50,0,20]; plotXTitles["SubleadPhoton_chargedHadronIso"] = "PF Iso_{abs}^{ch}(#gamma_{2}) [GeV]"
    plotNames["SubleadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_mvaID"] = [20,-1,1]; plotXTitles["SubleadPhoton_mvaID"] = "MVA ID(#gamma_{2})"
    plotNames["SubleadPhoton_pt/Diphoton_mass"] = ""; plotBins["SubleadPhoton_pt/Diphoton_mass"] = [50,0,2]; plotXTitles["SubleadPhoton_pt/Diphoton_mass"] = "p_{T}(#gamma_{2}) / M(#gamma#gamma)"

    plotNames["Diphoton_pt"] = ""; plotBins["Diphoton_pt"] = [50,0,500]; plotXTitles["Diphoton_pt"] = "p_{T}(#gamma#gamma) [GeV]"
    plotNames["Diphoton_eta"] = ""; plotBins["Diphoton_eta"] = [50,-3,3]; plotXTitles["Diphoton_eta"] = "#eta(#gamma#gamma)"
    plotNames["Diphoton_phi"] = ""; plotBins["Diphoton_phi"] = [50,-3.2,3.2]; plotXTitles["Diphoton_phi"] = "#phi(#gamma#gamma)"
    plotNames["Diphoton_mass"] = ""; plotBins["Diphoton_mass"] = [50,60,1000]; plotXTitles["Diphoton_mass"] = "M(#gamma#gamma) [GeV]"
    plotNames["Diphoton_mass_lowRange"] = "Diphoton_mass"; plotBins["Diphoton_mass_lowRange"] = [50,0,200]; plotXTitles["Diphoton_mass_lowRange"] = "M(#gamma#gamma) [GeV]"
    plotNames["Diphoton_pt_mgg"] = ""; plotBins["Diphoton_pt_mgg"] = [20,0,3]; plotXTitles["Diphoton_pt_mgg"] = "p_{T}(#gamma#gamma) / M(#gamma#gamma)"
    plotNames["Diphoton_dR"] = ""; plotBins["Diphoton_dR"] = [50,0,6]; plotXTitles["Diphoton_dR"] = "#DeltaR(#gamma#gamma)"

    plotNames["Diphoton_maxMvaID"] = "max(LeadPhoton_mvaID,SubleadPhoton_mvaID)"; plotBins["Diphoton_maxMvaID"] = [200,-1,1]; plotXTitles["Diphoton_maxMvaID"] = "max #gamma MVA ID"
    plotNames["Diphoton_minMvaID"] = "min(LeadPhoton_mvaID,SubleadPhoton_mvaID)"; plotBins["Diphoton_minMvaID"] = [200,-1,1]; plotXTitles["Diphoton_minMvaID"] = "min #gamma MVA ID"

    plotNames["n_jets"] = ""; plotBins["n_jets"] = [5,0,5]; plotXTitles["n_jets"] = "N_{jets}"

    plotNames["dijet_lead_pt"] = ""; plotBins["dijet_lead_pt"] = [50,0,500]; plotXTitles["dijet_lead_pt"] = "p_{T}(j_{1}) [GeV]"
    plotNames["dijet_lead_eta"] = ""; plotBins["dijet_lead_eta"] = [50,-3,3]; plotXTitles["dijet_lead_eta"] = "#eta(j_{1})"
    plotNames["dijet_lead_phi"] = ""; plotBins["dijet_lead_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_lead_phi"] = "#phi(j_{1})"
    plotNames["dijet_lead_mass"] = ""; plotBins["dijet_lead_mass"] = [50,0,100]; plotXTitles["dijet_lead_mass"] = "M(j_{1}) [GeV]"
    plotNames["dijet_lead_btagDeepFlavB"] = ""; plotBins["dijet_lead_btagDeepFlavB"] = [50,0,1]; plotXTitles["dijet_lead_btagDeepFlavB"] = "btagDeepFlavB(j_{1})"

    plotNames["dijet_sublead_pt"] = ""; plotBins["dijet_sublead_pt"] = [50,0,500]; plotXTitles["dijet_sublead_pt"] = "p_{T}(j_{2}) [GeV]"
    plotNames["dijet_sublead_eta"] = ""; plotBins["dijet_sublead_eta"] = [50,-3,3]; plotXTitles["dijet_sublead_eta"] = "#eta(j_{2})"
    plotNames["dijet_sublead_phi"] = ""; plotBins["dijet_sublead_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_sublead_phi"] = "#phi(j_{2})"
    plotNames["dijet_sublead_mass"] = ""; plotBins["dijet_sublead_mass"] = [50,0,100]; plotXTitles["dijet_sublead_mass"] = "M(j_{2}) [GeV]"
    plotNames["dijet_sublead_btagDeepFlavB"] = ""; plotBins["dijet_sublead_btagDeepFlavB"] = [50,0,1]; plotXTitles["dijet_sublead_btagDeepFlavB"] = "btagDeepFlavB(j_{2})"

    plotNames["dijet_pt"] = ""; plotBins["dijet_pt"] = [50,0,500]; plotXTitles["dijet_pt"] = "p_{T}(jj) [GeV]"
    plotNames["dijet_eta"] = ""; plotBins["dijet_eta"] = [50,-3,3]; plotXTitles["dijet_eta"] = "#eta(jj)"
    plotNames["dijet_phi"] = ""; plotBins["dijet_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_phi"] = "#phi(jj)"
    plotNames["dijet_mass"] = ""; plotBins["dijet_mass"] = [50,60,1000]; plotXTitles["dijet_mass"] = "M(jj) [GeV]"
    plotNames["dijet_dR"] = ""; plotBins["dijet_dR"] = [50,0,6]; plotXTitles["dijet_dR"] = "#DeltaR(jj)"

    plotNames["pfmet_pt"] = ""; plotBins["pfmet_pt"] = [50,0,200]; plotXTitles["pfmet_pt"] = "PF MET P_{T} [GeV]"
    plotNames["puppimet_pt"] = ""; plotBins["puppimet_pt"] = [50,0,200]; plotXTitles["puppimet_pt"] = "PUPPI MET P_{T} [GeV]"

    # 2D histos (y:x)
    plotNames["n_jets:Diphoton_maxMvaID"] = "n_jets:max(LeadPhoton_mvaID,SubleadPhoton_mvaID)"; plotBins["n_jets:Diphoton_maxMvaID"] = [200,-1,1,5,0,5]; plotXTitles["n_jets:Diphoton_maxMvaID"] = "max #gamma MVA ID"; plotYTitles["n_jets:Diphoton_maxMvaID"] = "N_{jets}"
    plotNames["n_jets:Diphoton_minMvaID"] = "n_jets:min(LeadPhoton_mvaID,SubleadPhoton_mvaID)"; plotBins["n_jets:Diphoton_minMvaID"] = [200,-1,1,5,0,5]; plotXTitles["n_jets:Diphoton_minMvaID"] = "min #gamma MVA ID"; plotYTitles["n_jets:Diphoton_minMvaID"] = "N_{jets}"

    #plotNames["xcand_pt:LeadPhoton_mvaID"] = ""; plotBins["xcand_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["xcand_pt:LeadPhoton_mvaID"] = "p_{T}(X) [GeV]"
    #plotNames["xcand_eta:LeadPhoton_mvaID"] = ""; plotBins["xcand_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["xcand_eta:LeadPhoton_mvaID"] = "#eta(X)"
    #plotNames["xcand_phi:LeadPhoton_mvaID"] = ""; plotBins["xcand_phi:LeadPhoton_mvaID"] = [50,-1,1,50,3.2,3.2]; plotXTitles["xcand_phi:LeadPhoton_mvaID"] = "#phi(X)"
    #plotNames["xcand_mass:LeadPhoton_mvaID"] = ""; plotBins["xcand_mass:LeadPhoton_mvaID"] = [50,-1,1,50,200,1000]; plotXTitles["xcand_mass:LeadPhoton_mvaID"] = "M(X) [GeV]"

    #plotNames["LeadPhoton_pt:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["LeadPhoton_pt:LeadPhoton_mvaID"] = "p_{T}(#gamma_{1}) [GeV]"
    #plotNames["LeadPhoton_eta:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["LeadPhoton_eta:LeadPhoton_mvaID"] = "#eta(#gamma_{1})"
    #plotNames["LeadPhoton_phi:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["LeadPhoton_phi:LeadPhoton_mvaID"] = "#phi(#gamma_{1})"
    #plotNames["LeadPhoton_pixelSeed:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_pixelSeed:LeadPhoton_mvaID"] = [50,-1,1,2,0,2]; plotXTitles["LeadPhoton_pixelSeed:LeadPhoton_mvaID"] = "hasPixelSeed(#gamma_{1})"
    #plotNames["LeadPhoton_r9:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_r9:LeadPhoton_mvaID"] = [50,-1,1,50,0,2]; plotXTitles["LeadPhoton_r9:LeadPhoton_mvaID"] = "R_{9}(#gamma_{1})"
    #plotNames["LeadPhoton_sieie:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_sieie:LeadPhoton_mvaID"] = [50,-1,1,50,0,0.05]; plotXTitles["LeadPhoton_sieie:LeadPhoton_mvaID"] = "#sigma_{ieie}(#gamma_{1})"
    #plotNames["LeadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = [50,-1,1,50,0,20]; plotXTitles["LeadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = "PF Iso_{abs}^{#gamma}(#gamma_{1}) [GeV]"
    #plotNames["LeadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = [50,-1,1,50,0,20]; plotXTitles["LeadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = "PF Iso_{abs}^{ch}(#gamma_{1}) [GeV]"
    #plotNames["LeadPhoton_mvaID:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_mvaID:LeadPhoton_mvaID"] = [50,-1,1,20,-1,1]; plotXTitles["LeadPhoton_mvaID:LeadPhoton_mvaID"] = "MVA ID(#gamma_{1})"
    #plotNames["LeadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = ""; plotBins["LeadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = [50,-1,1,50,0,2]; plotXTitles["LeadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = "p_{T}(#gamma_{1}) / M(#gamma#gamma)"

    #plotNames["SubleadPhoton_pt:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["SubleadPhoton_pt:LeadPhoton_mvaID"] = "p_{T}(#gamma_{2}) [GeV]"
    #plotNames["SubleadPhoton_eta:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["SubleadPhoton_eta:LeadPhoton_mvaID"] = "#eta(#gamma_{2})"
    #plotNames["SubleadPhoton_phi:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["SubleadPhoton_phi:LeadPhoton_mvaID"] = "#phi(#gamma_{2})"
    #plotNames["SubleadPhoton_pixelSeed:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_pixelSeed:LeadPhoton_mvaID"] = [50,-1,1,2,0,2]; plotXTitles["SubleadPhoton_pixelSeed:LeadPhoton_mvaID"] = "hasPixelSeed(#gamma_{2})"
    #plotNames["SubleadPhoton_r9:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_r9:LeadPhoton_mvaID"] = [50,-1,1,50,0,2]; plotXTitles["SubleadPhoton_r9:LeadPhoton_mvaID"] = "R_{9}(#gamma_{2})"
    #plotNames["SubleadPhoton_sieie:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_sieie:LeadPhoton_mvaID"] = [50,-1,1,50,0,0.05]; plotXTitles["SubleadPhoton_sieie:LeadPhoton_mvaID"] = "#sigma_{ieie}(#gamma_{2})"
    #plotNames["SubleadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = [50,-1,1,50,0,20]; plotXTitles["SubleadPhoton_pfPhoIso03:LeadPhoton_mvaID"] = "PF Iso_{abs}^{#gamma}(#gamma_{2}) [GeV]"
    #plotNames["SubleadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = [50,-1,1,50,0,20]; plotXTitles["SubleadPhoton_chargedHadronIso:LeadPhoton_mvaID"] = "PF Iso_{abs}^{ch}(#gamma_{2}) [GeV]"
    #plotNames["SubleadPhoton_mvaID:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_mvaID:LeadPhoton_mvaID"] = [50,-1,1,20,-1,1]; plotXTitles["SubleadPhoton_mvaID:LeadPhoton_mvaID"] = "MVA ID(#gamma_{2})"
    #plotNames["SubleadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = ""; plotBins["SubleadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = [50,-1,1,50,0,2]; plotXTitles["SubleadPhoton_pt/Diphoton_mass:LeadPhoton_mvaID"] = "p_{T}(#gamma_{2}) / M(#gamma#gamma)"

    #plotNames["Diphoton_pt:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["Diphoton_pt:LeadPhoton_mvaID"] = "p_{T}(#gamma#gamma) [GeV]"
    #plotNames["Diphoton_eta:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["Diphoton_eta:LeadPhoton_mvaID"] = "#eta(#gamma#gamma)"
    #plotNames["Diphoton_phi:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["Diphoton_phi:LeadPhoton_mvaID"] = "#phi(#gamma#gamma)"
    #plotNames["Diphoton_mass:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_mass:LeadPhoton_mvaID"] = [50,-1,1,50,60,1000]; plotXTitles["Diphoton_mass:LeadPhoton_mvaID"] = "M(#gamma#gamma) [GeV]"
    #plotNames["Diphoton_pt_mgg:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_pt_mgg:LeadPhoton_mvaID"] = [50,-1,1,20,0,3]; plotXTitles["Diphoton_pt_mgg:LeadPhoton_mvaID"] = "p_{T}(#gamma#gamma) / M(#gamma#gamma)"
    #plotNames["Diphoton_dR:LeadPhoton_mvaID"] = ""; plotBins["Diphoton_dR:LeadPhoton_mvaID"] = [50,-1,1,50,0,6]; plotXTitles["Diphoton_dR:LeadPhoton_mvaID"] = "#DeltaR(#gamma#gamma)"

    #plotNames["n_jets:LeadPhoton_mvaID"] = ""; plotBins["n_jets:LeadPhoton_mvaID"] = [50,-1,1,5,0,5]; plotXTitles["n_jets:LeadPhoton_mvaID"] = "N_{jets}"

    #plotNames["dijet_lead_pt:LeadPhoton_mvaID"] = ""; plotBins["dijet_lead_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["dijet_lead_pt:LeadPhoton_mvaID"] = "p_{T}(j_{1}) [GeV]"
    #plotNames["dijet_lead_eta:LeadPhoton_mvaID"] = ""; plotBins["dijet_lead_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["dijet_lead_eta:LeadPhoton_mvaID"] = "#eta(j_{1})"
    #plotNames["dijet_lead_phi:LeadPhoton_mvaID"] = ""; plotBins["dijet_lead_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["dijet_lead_phi:LeadPhoton_mvaID"] = "#phi(j_{1})"
    #plotNames["dijet_lead_mass:LeadPhoton_mvaID"] = ""; plotBins["dijet_lead_mass:LeadPhoton_mvaID"] = [50,-1,1,50,0,100]; plotXTitles["dijet_lead_mass:LeadPhoton_mvaID"] = "M(j_{1}) [GeV]"
    #plotNames["dijet_lead_btagDeepFlavB:LeadPhoton_mvaID"] = ""; plotBins["dijet_lead_btagDeepFlavB:LeadPhoton_mvaID"] = [50,-1,1,50,0,1]; plotXTitles["dijet_lead_btagDeepFlavB:LeadPhoton_mvaID"] = "btagDeepFlavB(j_{1})"

    #plotNames["dijet_sublead_pt:LeadPhoton_mvaID"] = ""; plotBins["dijet_sublead_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["dijet_sublead_pt:LeadPhoton_mvaID"] = "p_{T}(j_{2}) [GeV]"
    #plotNames["dijet_sublead_eta:LeadPhoton_mvaID"] = ""; plotBins["dijet_sublead_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["dijet_sublead_eta:LeadPhoton_mvaID"] = "#eta(j_{2})"
    #plotNames["dijet_sublead_phi:LeadPhoton_mvaID"] = ""; plotBins["dijet_sublead_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["dijet_sublead_phi:LeadPhoton_mvaID"] = "#phi(j_{2})"
    #plotNames["dijet_sublead_mass:LeadPhoton_mvaID"] = ""; plotBins["dijet_sublead_mass:LeadPhoton_mvaID"] = [50,-1,1,50,0,100]; plotXTitles["dijet_sublead_mass:LeadPhoton_mvaID"] = "M(j_{2}) [GeV]"
    #plotNames["dijet_sublead_btagDeepFlavB:LeadPhoton_mvaID"] = ""; plotBins["dijet_sublead_btagDeepFlavB:LeadPhoton_mvaID"] = [50,-1,1,50,0,1]; plotXTitles["dijet_sublead_btagDeepFlavB:LeadPhoton_mvaID"] = "btagDeepFlavB(j_{2})"

    #plotNames["dijet_pt:LeadPhoton_mvaID"] = ""; plotBins["dijet_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,500]; plotXTitles["dijet_pt:LeadPhoton_mvaID"] = "p_{T}(jj) [GeV]"
    #plotNames["dijet_eta:LeadPhoton_mvaID"] = ""; plotBins["dijet_eta:LeadPhoton_mvaID"] = [50,-1,1,50,-3,3]; plotXTitles["dijet_eta:LeadPhoton_mvaID"] = "#eta(jj)"
    #plotNames["dijet_phi:LeadPhoton_mvaID"] = ""; plotBins["dijet_phi:LeadPhoton_mvaID"] = [50,-1,1,50,-3.2,3.2]; plotXTitles["dijet_phi:LeadPhoton_mvaID"] = "#phi(jj)"
    #plotNames["dijet_mass:LeadPhoton_mvaID"] = ""; plotBins["dijet_mass:LeadPhoton_mvaID"] = [50,-1,1,50,60,1000]; plotXTitles["dijet_mass:LeadPhoton_mvaID"] = "M(jj) [GeV]"
    #plotNames["dijet_dR:LeadPhoton_mvaID"] = ""; plotBins["dijet_dR:LeadPhoton_mvaID"] = [50,-1,1,50,0,6]; plotXTitles["dijet_dR:LeadPhoton_mvaID"] = "#DeltaR(jj)"

    #plotNames["pfmet_pt:LeadPhoton_mvaID"] = ""; plotBins["pfmet_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,200]; plotXTitles["pfmet_pt:LeadPhoton_mvaID"] = "PF MET P_{T} [GeV]"
    #plotNames["puppimet_pt:LeadPhoton_mvaID"] = ""; plotBins["puppimet_pt:LeadPhoton_mvaID"] = [50,-1,1,50,0,200]; plotXTitles["puppimet_pt:LeadPhoton_mvaID"] = "PUPPI MET P_{T} [GeV]"

    toinclude = []; toexclude = []
    if args.DDHistos:
      toinclude = ["Diphoton_maxMvaID", "Diphoton_minMvaID", "n_jets:Diphoton_maxMvaID", "n_jets:Diphoton_minMvaID"]
    else:
      toexclude = ["n_jets:Diphoton_maxMvaID", "n_jets:Diphoton_minMvaID"]
    if len(args.plotsToInclude) > 0: toinclude = args.plotsToInclude
    if len(args.plotsToExclude) > 0: toexclude = args.plotsToExclude


    if args.yields:
      plotname = "LeadPhoton_pixelSeed"
      plotDict = get_plots(samples, year, plotname, cut, plotBins[plotname], plotXTitles[plotname], plotYTitles[plotname] if plotname in plotYTitles.keys() else "Event")
      get_yields(plotDict, plotname, lumi, year, args.data)
    else:
      for plotname in (plotNames.keys() if len(toinclude)==0 else toinclude):
        if plotname in toexclude:
            continue
        # Open files and get trees and return plots
        plotDict = get_plots(samples, year, plotname if plotNames[plotname]=="" else plotNames[plotname], cut, plotBins[plotname], plotXTitles[plotname], plotYTitles[plotname] if plotname in plotYTitles.keys() else "Event")
        #plotDict = get_plots(samples, year, plotname if plotNames[plotname]=="" else plotNames[plotname], cut, plotBins[plotname], "MVA ID(#gamma_{1})", plotXTitles[plotname])
        draw_plot(plotDict, plotname, lumi, year, not args.noLogY, args.logX, args.data)
