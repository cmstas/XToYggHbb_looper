import sys,os,copy
import math
import ROOT 

### This script returns a cutflow table in tex format.
### Then:
### pdflatex <filename>.tex
### pdfcrop --margins '0 0' <filename>.pdf
### mv <filename-crop>.pdf <filename>.pdf 

usePunziDefinition = True

def print_header(fout,allsamples=[],allsampleLabels=[],doSoverB=False):
    fout.write('\\documentclass{article}\n')  
    fout.write('\\usepackage{adjustbox}\n')
    fout.write('\\usepackage{multirow}\n')
    fout.write('\\thispagestyle{empty}\n')
    fout.write('\\begin{document}\n')
    fout.write('\\begin{table*}[h]\n')
    fout.write('\\footnotesize\n')
    fout.write('\\begin{adjustbox}{width=\\textwidth}\n')
    fout.write('\\begin{tabular}{|l')
    if not doSoverB:
        for i in range(len(allsamples)):
            # Add one cell for absolute yield and one cell for efficiency
            fout.write('|c|c')
        fout.write('|}\n')
        fout.write('\\hline\n')
        fout.write('Selection ')
        for i in range(len(allsamples)):
            fout.write('& \\multicolumn{2}{|c|}{'+allsampleLabels[i].replace("_"," ")+'}')
        fout.write('\\\\\n')
        for i in range(len(allsamples)):
            fout.write('& Yield & Eff. (\%)')
    else:
        for i in range(len(allsamples)):
            # Add one cell for absolute yield and one cell for s/sqrt(b), except for total SM:
            if "Total SM" not in allsamples[i]:
                fout.write('|c|c')
            else:
                fout.write('|c')
        fout.write('|}\n')
        fout.write('\\hline\n')
        fout.write('Selection ')
        for i in range(len(allsamples)):
            if "Total SM" not in allsamples[i]:
                fout.write('& \\multicolumn{2}{|c|}{'+allsampleLabels[i].replace("_"," ")+'}')
            else:
                fout.write('& '+allsampleLabels[i])                
        fout.write('\\\\\n')
        for i in range(len(allsamples)):
            if "Total SM" not in allsamples[i]:
                if not usePunziDefinition:
                    fout.write('& Yield & S/$\\sqrt{\\mathrm{B}}$')
                else:
                    fout.write('& Yield & S/($3/2+\\sqrt{\\mathrm{B}}$)')
            else:
                fout.write('& Yield')
    fout.write('\\\\\n')
    fout.write('\\hline\n')
    fout.write('\\hline\n')

def print_footer(fout):
    fout.write('\\end{tabular}\n')
    fout.write('\\end{adjustbox}\n')
    fout.write('\\end{table*}\n')
    fout.write('\\end{document}\n')

def reformat_label(tlabel):
    tlabel = tlabel.replace("pT","$p_{\\mathrm{T}}$")
    tlabel = tlabel.replace("dR","$\\Delta$R")
    tlabel = tlabel.replace("&","and")
    tlabel = tlabel.replace("|eta|","$|\\eta|$")
    tlabel = tlabel.replace("m_{#mu#mu}","$m_{\\mu\\mu}$")
    tlabel = tlabel.replace("m_{#mu b}","$m_{\\mu\\mathrm{b}}$")
    tlabel = tlabel.replace("#mu","$\\mu$")
    tlabel = tlabel.replace("N_{b-tag}","$N_{\\mathrm{b-tag}}$")
    tlabel = tlabel.replace("p_{T}","$p_{\\mathrm{T}}$")
    tlabel = tlabel.replace("E_{T}^{miss}","$E_{\\mathrm{T}}^{\\mathrm{miss}}$")
    tlabel = tlabel.replace("#geq","$\\geq$")
    tlabel = tlabel.replace("=","$=$")
    tlabel = tlabel.replace(">","$>$")
    tlabel = tlabel.replace("<","$<$")
    tlabel = tlabel.replace("$$","$ $")
    return tlabel

def make_cutflow_table(cutflow="cutflow", samples=[], sampleLabels=[], indir = "./cpp/temp_data/", year = "2018", outdir="tables/", extension="", doBkgTable=True, doSignalTable=False, doSoverB=False, doSignalOnlyTable=False):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if len(sampleLabels)<len(samples):
        sampleLabels = samples

    issignal = []
    bkgsamples = []
    sigsamples = []
    bkgyields  = []
    sigyields  = []
    bkglabels  = []
    siglabels  = []
    cutlabels  = []
    totsm      = []
    totsmeffs  = []
    bkgeffs    = []
    sigeffs    = []
    soverb     = []
    #
    nCuts = 0
    iS = 0
    iB = 0
    for i in range(len(samples)):
        subsamples = []
        if samples[i]=='ZToMuMu':
            subsamples = ['ZToMuMu_50_120',
                          'ZToMuMu_120_200',
                          'ZToMuMu_200_400',
                          'ZToMuMu_400_800',
                          'ZToMuMu_800_1400',
                          'ZToMuMu_1400_2300',
                          'ZToMuMu_2300_3500',
                          'ZToMuMu_3500_4500',
                          'ZToMuMu_4500_6000',
                          'ZToMuMu_6000_Inf']
        elif samples[i]=='VV':
            subsamples = ['WW','WZ','ZZ']
        elif samples[i]=='tW+tZq':
            subsamples = ['tW','tbarW','tZq']
        elif samples[i]=='ttX':
            subsamples = ['TTW','TTZ','TTHToNonbb','TTHTobb']
        else:
            subsamples.append(samples[i])

        if "Y3_M" in samples[i] or "DY3_M" in samples[i] or "DYp3_M" in samples[i] or "B3mL2_M" in samples[i]:
            issignal.append(True)
            sigsamples.append(samples[i])
            siglabels.append(sampleLabels[i])
            sigyields.append([])
            sigeffs  .append([])
            iS = iS+1
        else:
            issignal.append(False)
            bkgsamples.append(samples[i])
            bkglabels.append(sampleLabels[i])
            bkgyields.append([])
            bkgeffs  .append([])
            iB = iB+1

        for ii in range(len(subsamples)):
            tfile = ROOT.TFile(indir+"/output_"+subsamples[ii]+"_"+year+".root") 
            thist = tfile.Get(cutflow).Clone(subsamples[ii])
            
            if i==0 and ii==0:
                nCuts = thist.GetNbinsX()
                for b in range(1, nCuts+1):
                    if thist.GetXaxis().GetBinLabel(b)=="":
                        nCuts = nCuts-1
                        continue
                    else:
                        cutlabels.append(thist.GetXaxis().GetBinLabel(b))

            for b in range(1, nCuts+1):
                if ii==0:
                    if issignal[i]:
                        sigyields[iS-1].append(thist.GetBinContent(b))
                    else:
                        bkgyields[iB-1].append(thist.GetBinContent(b))
                else:
                    if issignal[i]:
                        sigyields[iS-1][b-1] = sigyields[iS-1][b-1]+thist.GetBinContent(b)
                    else:
                        bkgyields[iB-1][b-1] = bkgyields[iB-1][b-1]+thist.GetBinContent(b)

        for b in range(1, nCuts+1):
            if i==0 and not issignal[i]:
                totsm.append(0.0)
            if not issignal[i]:
                totsm[b-1] = totsm[b-1]+bkgyields[iB-1][b-1]
                if bkgyields[iB-1][0]>0.0:
                    bkgeffs[iB-1].append(bkgyields[iB-1][b-1]/bkgyields[iB-1][0])
                else:
                    bkgeffs[iB-1].append(0.0)
            else:
                if sigyields[iS-1][0]>0.0:
                    sigeffs[iS-1].append(sigyields[iS-1][b-1]/sigyields[iS-1][0])
                else:
                    sigeffs[iS-1].append(0.0)

    for b in range(1, nCuts+1):
        if totsm[0]>0.0:
            totsmeffs.append(totsm[b-1]/totsm[0])
        else:
            totsmeffs.append(0.0)

    for i in range(len(sigsamples)):
        soverb.append([])
        for b in range(1, nCuts+1):
            tsoverb = 1e9
            if totsm[b-1]>0.0:
                if not usePunziDefinition:
                    tsoverb = sigyields[i][b-1]/ROOT.TMath.Sqrt(totsm[b-1])
                else:
                    # https://arxiv.org/pdf/physics/0308063.pdf
                    tsoverb = sigyields[i][b-1]/(1.5+ROOT.TMath.Sqrt(totsm[b-1]))
            soverb[i].append(tsoverb)

    canDoSignalTable = False
    canDoBkgTable = False
    if True in issignal:
        canDoSignalTable = True
    if False in issignal:
        canDoBkgTable = True
    if not canDoSignalTable and doSignalTable:
        print("Can not produce a signal table with no signal")
        exit()
    if not canDoBkgTable and doBkgTable:
        print("Can not produce a background table with no background")
        exit()
    if not canDoBkgTable and doSoverB:
        print("Can not produce S/sqrt(B) table with no background")
        exit()
    if not canDoSignalTable and doSoverB:
        print("Can not produce S/sqrt(B) table with no signal")
    if doBkgTable and doSoverB:
        print("Can not append S/sqrt(B) for background, will only show total SM absolute yields")
        doBkgTable = False
    if not canDoBkgTable and canDoSignalTable and doSignalTable:
        doSignalTable = False
        doSignalOnlyTable = True
    if doSoverB:
        doSignalTable = True
        doSignalOnlyTable = False

    skimLabel = []
    skimLabel.append('$\\geq 2$ muons, $p_{\\mathrm{T}}$(leading)$>$50 GeV')
    skimLabel.append('$\\geq 1$ pair with m(ll)$>$100 GeV')

    tabletype = ""
    if doSoverB:
        tabletype=tabletype+"SoverB"
    if doBkgTable:
        tabletype=tabletype+"bkg"
    if doSignalTable:
        tabletype=tabletype+"sig"
    if doSignalOnlyTable:
        tabletype=tabletype+"sigOnly"

    allsamples = []
    allsampleLabels = []
    allyields = []
    alleffs = []
    if not doSignalOnlyTable:
        if doBkgTable:
            allsamples = bkgsamples
            allsampleLabels = bkglabels
            allyields = bkgyields
            alleffs = bkgeffs
        allsamples.append("Total SM")
        allsampleLabels.append("Total SM")
        allyields.append(totsm)
        alleffs.append(totsmeffs)
    if doSignalTable or doSignalOnlyTable:
        allsamples = allsamples+sigsamples
        allsampleLabels = allsampleLabels+siglabels
        allyields = allyields+sigyields
        alleffs = alleffs+sigeffs

    foutname = cutflow+"_"+year+"_"+tabletype
    if extension != "":
        foutname = foutname+"_"+extension
    fout = open(outdir+"/"+foutname+".tex",'w')

    if not doSoverB:
        print_header(fout,allsamples,allsampleLabels,doSoverB)
        for c in range(nCuts):
            tlabel = cutlabels[c]
            tlabel = reformat_label(tlabel)
            if "skim" in tlabel and len(skimLabel)==1:
                tlabel = skimLabel[0]
            if "skim" in tlabel and len(skimLabel)>1:
                for ln,l in enumerate(skimLabel):
                    if ln==0:
                        fout.write(l)
                        for i in range(len(allsamples)):
                            fout.write('& \\multirow{%d}{*}{%.2E} & \\multirow{%d}{*}{%.2E}'%(len(skimLabel),allyields[i][c],len(skimLabel),100.0*alleffs[i][c]))
                        fout.write('\\\\\n')
                    else:
                        fout.write(l)
                        for i in range(len(allsamples)):
                            fout.write('&  & ')
                        fout.write('\\\\\n')
                        fout.write('\\hline\n')
            else:
                fout.write(tlabel)
                for i in range(len(allsamples)):
                    fout.write('& %.2E & %.2E'%(allyields[i][c],100.0*alleffs[i][c]))
                fout.write('\\\\\n')
                fout.write('\\hline\n')
        print_footer(fout)

    else:
        print_header(fout,allsamples,allsampleLabels,doSoverB)
        for c in range(nCuts):
            tlabel = cutlabels[c]
            tlabel = reformat_label(tlabel)
            if "skim" in tlabel and len(skimLabel)==1:
                tlabel = skimLabel[0]
                print(tlabel)
            if "skim" in tlabel and len(skimLabel)>1:
                for ln,l in enumerate(skimLabel):
                    if ln==0:
                        fout.write(l)
                        iS = 0
                        for i in range(len(allsamples)):
                            if "Total SM" not in allsamples[i]:
                                iS = iS+1
                                fout.write('& \\multirow{%d}{*}{%.2E} & \\multirow{%d}{*}{%.2E}'%(len(skimLabel),allyields[i][c],len(skimLabel),soverb[iS-1][c]))
                            else:
                                fout.write('& \\multirow{%d}{*}{%.2E}'%(len(skimLabel),allyields[i][c]))
                        fout.write('\\\\\n')
                    else:
                        fout.write(l)
                        for i in range(len(allsamples)):
                            if "Total SM" not in allsamples[i]:
                                fout.write('&  & ')
                            else:
                                fout.write('& ')
                        fout.write('\\\\\n')
                        fout.write('\\hline\n')
            else:
                fout.write(tlabel)
                iS = 0
                for i in range(len(allsamples)):
                    if "Total SM" not in allsamples[i]:
                        iS = iS+1
                        fout.write('& %.2E & %.2E'%(allyields[i][c],soverb[iS-1][c]))
                    else:
                        fout.write('& %.2E'%(allyields[i][c]))                        
                fout.write('\\\\\n')
                fout.write('\\hline\n')
        print_footer(fout)

    fout.close()


def getListOfCutflows(infile="./cpp/temp_data/output_ttbar_2018.root"):
    listofccutflows = []
    listfile = ROOT.TFile(infile)
    listkeys = listfile.GetListOfKeys()
    size = listkeys.GetSize()
    for i in range(0,size):
        if "cutflow" not in listkeys.At(i).GetName():
            continue
        listofcutflows.append(listkeys.At(i).GetName())
    return listofcutflows

# SM backgrounds
bkgsamples=['ttbar','tW+tZq','ttX','ZToMuMu','VV']
bkgsampleLabels=["t$\\bar{\\mathrm{t}}$","tW+tZq","t$\\bar{\\mathrm{t}}$X","DY($\\mu\\mu$)","VV (V$=$Z,W)"]

# Signal
sigsamples=[]
sigsampleLabels=[]
# For tables of 'standard' signals across mass points, it is enough to change signalname below; accepted values: "Y3", "DY3", "DYp3", "B3mL2"
signalname = "Y3" 
if signalname=="Y3":
    sigsamples=['Y3_M200','Y3_M400','Y3_M700','Y3_M1000','Y3_M1500','Y3_M2000']
elif signalname=="DY3":
    sigsamples=['DY3_M200','DY3_M400','DY3_M700','DY3_M1000','DY3_M1500','DY3_M2000']
elif signalname=="DYp3":
    sigsamples=['DYp3_M200','DYp3_M400','DYp3_M700','DYp3_M1000','DYp3_M1500','DYp3_M2000']
elif signalname=="B3mL2":
    sigsamples=['B3mL2_M200','B3mL2_M400','B3mL2_M700','B3mL2_M1000','B3mL2_M1500','B3mL2_M2000']
else:
    print("Signal is unknown: please, explicitly set your sigsamples and sigsampleLabels lists")
    exit()
sigsampleLabels=["%s (%s GeV)"%(i.split("_")[0],i.split("_")[1].replace("M","")) for i in sigsamples]

# SM + signal
samples=bkgsamples+sigsamples
sampleLabels=bkgsampleLabels+sigsampleLabels

year="2018"
#indir="./cpp/temp_data/"
indir="./cpp/test_17Feb23/"
outdir="tables/"
listofcutflows = []
listofcutflows = getListOfCutflows(indir+"output_ttbar_"+year+".root")
toexclude = []

###make_cutflow_table(cutflow="cutflow", samples=[], sampleLabels=[], indir = "./cpp/temp_data/", year = "2018", outdir="tables/", extension="", doBkgTable=True, doSignalTable=False, doSoverB=False, doSignalOnlyTable=False)
for cutflow in listofcutflows:
    if cutflow in toexclude:
        continue
    make_cutflow_table(cutflow,samples,sampleLabels,year=year,indir=indir,outdir=outdir,doBkgTable=True,doSignalTable=False,doSoverB=False,doSignalOnlyTable=False)
    make_cutflow_table(cutflow,samples,sampleLabels,year=year,indir=indir,outdir=outdir,doBkgTable=False,doSignalTable=True,doSoverB=False,doSignalOnlyTable=False)
    make_cutflow_table(cutflow,samples,sampleLabels,year=year,indir=indir,outdir=outdir,doBkgTable=False,doSignalTable=True,doSoverB=False,doSignalOnlyTable=True,extension=signalname)
    make_cutflow_table(cutflow,samples,sampleLabels,year=year,indir=indir,outdir=outdir,doBkgTable=False,doSignalTable=True,doSoverB=True,doSignalOnlyTable=False,extension=signalname)
