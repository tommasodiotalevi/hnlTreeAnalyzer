import sys
import argparse
import ROOT
from ctypes import *

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser(description="")
parser.add_argument("--data_filename"     , type=str  ,            help="Path to the input data rootfile")
parser.add_argument("--prompt_filename"   , type=str  ,            help="Path to the input prompt rootfile")
parser.add_argument("--nonprompt_filename", type=str  ,            help="Path to the input nonprompt rootfile")
parser.add_argument("--intLumi"           , type=float,            help="Print the integrated lumi")
parser.add_argument("--addTag"            , type=str  ,default="", help="Tag the output file")
args = parser.parse_args()

ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.035)

dataFileName = args.data_filename
mcPromptFileName = args.prompt_filename
mcNonPromptFileName = args.nonprompt_filename
histoName = "h_ds_l_prop"
tag = args.addTag

bpark_lumi = args.intLumi*10e3
prescale_factor = 0.19
BR_BToDs = 0.20509380
BR_Ds = 1.3e-5
prompt_fraction_in_mc = 0.7917
nonprompt_fraction_in_mc = 0.2083

n_era = args.data_filename[args.data_filename.find("ParkingBPH")+len("ParkingBPH")]
era   = args.data_filename[args.data_filename.find("Run2018")+len("Run2018")]
int_lumi = args.intLumi

ww = bpark_lumi*prescale_factor*BR_Ds
print("MC weight: {}".format(ww))

dataFile = ROOT.TFile.Open(dataFileName)
mcPromptFile = ROOT.TFile.Open(mcPromptFileName)
mcNonPromptFile = ROOT.TFile.Open(mcNonPromptFileName)

data = dataFile.Get(histoName)          # data histogram
mc0  = mcPromptFile.Get(histoName)      # first MC histogram
mc1  = mcNonPromptFile.Get(histoName)   # second MC histogram
mc0_wOverflow = ROOT.TH1F("mc0_wOverflow",";D_{s} reconstructed proper decay length [cm];Events",21,0.,0.105)
mc1_wOverflow = ROOT.TH1F("mc1_wOverflow",";D_{s} reconstructed proper decay length [cm];Events",21,0.,0.105)
data_wOverflow = ROOT.TH1F("data_wOverflow",";D_{s} reconstructed proper decay length [cm];Events",21,0.,0.105)
for i in range(1,22):
    mc0_wOverflow.SetBinContent(i,mc0.GetBinContent(i))
    mc1_wOverflow.SetBinContent(i,mc1.GetBinContent(i))
    data_wOverflow.SetBinContent(i,data.GetBinContent(i))
    mc0_wOverflow.SetBinError(i,mc0.GetBinError(i))
    mc1_wOverflow.SetBinError(i,mc1.GetBinError(i))
    data_wOverflow.SetBinError(i,data.GetBinError(i))
mc0_wOverflow.Scale(ww)
mc1_wOverflow.Scale(ww)
mc = ROOT.TObjArray(2)                  # MC histograms are put in this array
mc.Add(mc0_wOverflow)
mc.Add(mc1_wOverflow)
fit = ROOT.TFractionFitter(data_wOverflow, mc) # initialise
fit.Constrain(0,0.,1.)             # constrain fraction 1 to be between 0 and 1
status = int(fit.Fit())              # perform the fit

print("fit status: {}".format(status))
if status == 0:                       # check on fit status 
    print("Chi2/ndof={}".format(fit.GetChisquare()/float(fit.GetNDF())))
    print("plotting result")
    result = fit.GetPlot()
    result.SetName("result")
    f0  = c_double()
    f0e = c_double()
    f1  = c_double()
    f1e = c_double()
    fit.GetResult(0,f0,f0e)
    fit.GetResult(1,f1,f1e)
    #data.Draw("Ep")
    #result.Draw("same")

    #fit
    c_fit = ROOT.TCanvas("c_fit","c_fit",600,600)
    c_fit.SetLeftMargin(0.15)
    c_fit.SetRightMargin(0.05)

    result.SetLineWidth(2)
    result.SetLineColor(ROOT.kBlue)
    data_wOverflow.SetMarkerStyle(20)
    data_wOverflow.SetMarkerColor(ROOT.kBlack)
    data_wOverflow.SetLineColor(ROOT.kBlack)
    ymax = max([data_wOverflow.GetMaximum(),result.GetMaximum()])
    data_wOverflow.SetMaximum(ymax+(ymax/3.))
    result.SetMaximum(ymax+(ymax/3.))
    result.GetYaxis().SetTitleOffset(2.2)
    result.GetXaxis().SetTitle("D_{s} proper decay length [cm]")
    data_wOverflow.GetXaxis().SetTitle("D_{s} proper decay length [cm]")
    result.Draw("hist")
    data_wOverflow.Draw("p same")
    leg_fit = ROOT.TLegend(0.55, 0.65, 0.92, 0.87)
    leg_fit.AddEntry("data_wOverflow","Bkg subtracted data","ep")
    leg_fit.AddEntry("result","Fit","l")
    leg_fit.SetBorderSize(0)
    leg_fit.Draw("same")
    latex = ROOT.TLatex()
    latex.SetTextAlign(12)
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(0.7,0.92,str(int_lumi)+" fb^{-1} (13 TeV)")
    latex.DrawLatexNDC(0.55,0.55,"Direct D_{s} fraction: "+str(round(f0.value,3))+"#pm"+str(round(f0e.value,3)))
    latex.DrawLatexNDC(0.55,0.5,"B #rightarrow D_{s} fraction: "+str(round(f1.value,3))+"#pm"+str(round(f1e.value,3)))
    c_fit.SaveAs("fraction_"+era+n_era+"_"+tag+"_fit_plot.png")
    c_fit.SaveAs("fraction_"+era+n_era+"_"+tag+"_fit_plot.pdf")

    #template
    c_tem = ROOT.TCanvas("c_tem","c_tem",600,600)
    c_tem.SetLeftMargin(0.15)
    c_tem.SetRightMargin(0.05)

    mc0_wOverflow.Scale(1./mc0_wOverflow.Integral())
    mc1_wOverflow.Scale(1./mc1_wOverflow.Integral())
    mc0_wOverflow.SetLineColor(ROOT.kGreen)
    mc1_wOverflow.SetLineColor(ROOT.kRed)
    mc0_wOverflow.SetLineStyle(ROOT.kDashed)
    mc1_wOverflow.SetLineStyle(ROOT.kDashed)
    mc0_wOverflow.SetLineWidth(2)
    mc1_wOverflow.SetLineWidth(2)
    mc0_wOverflow.GetYaxis().SetTitle("Normalized to unit")
    mc1_wOverflow.GetYaxis().SetTitle("Normalized to unit")
    mc0_wOverflow.GetXaxis().SetTitle("D_{s} proper decay length [cm]")
    mc1_wOverflow.GetXaxis().SetTitle("D_{s} proper decay length [cm]")
    mc0_wOverflow.Draw("hist")
    mc1_wOverflow.Draw("hist same")
    leg_tem = ROOT.TLegend(0.55, 0.65, 0.92, 0.87)
    leg_tem.AddEntry("mc0_wOverflow","Direct D_{s}","l")
    leg_tem.AddEntry("mc1_wOverflow","B #rightarrow D_{s}","l")
    leg_tem.SetBorderSize(0)
    leg_tem.Draw("same")
    c_tem.SaveAs("fraction_"+era+n_era+"_"+tag+"_tem_plot.png")
    c_tem.SaveAs("fraction_"+era+n_era+"_"+tag+"_tem_plot.pdf")

