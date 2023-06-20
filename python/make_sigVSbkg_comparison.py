import subprocess
import sys
import os
import json
import argparse

#script input arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("cfg_filename"      ,                          help="Path to the input configuration file")
parser.add_argument("--logy" , action='store_true', default=False, help="Plot with log y scale")
args = parser.parse_args()

configFileName = args.cfg_filename

with open(configFileName, "r") as f:
    config = json.loads(f.read())
import ROOT

ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.02)

inputDirName = str(config["inputDirName"])
outDirName = str(config["outDirName"])

for plotName in config["plotNameList"]:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c","c",800,800)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    
    histoName = str(plotName)
    print("--> {}".format(histoName))
    colorList = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kViolet, ROOT.kPink]
    iColor = 0
    iLineStyle=0
    leg_sigVSbkg = ROOT.TLegend(0.5, 0.65, 0.89, 0.89)
    
    histoStacked = ROOT.THStack("histoStacked","histoStacked")
    
    keyList_bkg = [key for key in config["background"]]
    inputHistoList_bkg = [ROOT.TH1D() for key in config["background"]]
    inputFileList_bkg = [ROOT.TFile() for key in config["background"]]
    inputHistoDic_bkg = dict(zip(keyList_bkg,inputHistoList_bkg))
    inputFileDic_bkg  = dict(zip(keyList_bkg,inputFileList_bkg))
    
    keyList_sig = [key for key in config["signal"]]
    inputHistoList_sig = [ROOT.TH1D() for key in config["signal"]]
    inputFileList_sig = [ROOT.TFile() for key in config["signal"]]
    inputHistoDic_sig = dict(zip(keyList_sig,inputHistoList_sig))
    inputFileDic_sig  = dict(zip(keyList_sig,inputFileList_sig))

    xaxis_label = str()
    integral_tot = float()
    hmax = float()
    hmax_sig = float()
    

    for filename in config["background"]:
        inputFileDic_bkg [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_bkg[filename] = ROOT.TH1D(inputFileDic_bkg[filename].Get(histoName))
        integral_tot += float(inputHistoDic_bkg[filename].Integral())
        #print("--> filename: {}".format(filename))
        #print("--> integral tot: {}".format(integral_tot))
    
    #Stacking background histos
    for filename in config["background"]:
        inputFileDic_bkg [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_bkg[filename] = ROOT.TH1D(inputFileDic_bkg[filename].Get(histoName))
        inputHistoDic_bkg[filename].Scale(1./integral_tot)
        inputHistoDic_bkg[filename].SetLineColor(ROOT.kBlack)
        inputHistoDic_bkg[filename].SetFillColor(ROOT.kGreen+iColor)
        #overflowBin = inputHistoDic_bkg[filename].GetXaxis().GetLast() + 1
        overflowBin = inputHistoDic_bkg[filename].GetXaxis().GetLast() # This is not the overflow bin
        inputHistoDic_bkg[filename].GetXaxis().SetRange(1,overflowBin)
        xaxis_label = str(inputHistoDic_bkg[filename].GetXaxis().GetTitle())
        histoStacked.Add(inputHistoDic_bkg[filename])
        iColor += 1
        bkgLabel = config["background"][filename]["label"]
        leg_sigVSbkg.AddEntry(inputHistoDic_bkg[filename],bkgLabel)

    c.cd()
    if args.logy:
        c.SetLogy()
    histoStacked.Draw() #otherwise next commands before Draw will fail
    histoStacked.SetTitle("")
    histoStacked.GetYaxis().SetTitle("Normalized to unit")
    #histoStacked.SetMinimum(0.001)
    histoStacked.GetXaxis().SetTitle(xaxis_label)
    if xaxis_label.find("IPS")>0:
        xaxis_label = xaxis_label.replace("[cm] ","")
        histoStacked.GetXaxis().SetTitle(xaxis_label)
    hmax = float(histoStacked.GetMaximum())

    #Superimposing signal
    for filename in config["signal"]:
        inputFileDic_sig [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_sig[filename] = ROOT.TH1D(inputFileDic_sig[filename].Get(histoName))
        integral = float(inputHistoDic_sig[filename].Integral())
        #print("----> signal filename: {}".format(filename))
        #print("----> signal integral: {}".format(integral))
        inputHistoDic_sig[filename].Scale(1./integral)
        hmax_sig = float(inputHistoDic_sig[filename].GetMaximum())
        if hmax_sig>hmax:
            hmax = hmax_sig
    yedge = hmax+(hmax/2)

    histoStacked.SetMaximum(yedge)
    if plotName.find("2DDist")>0:
        histoStacked.GetXaxis().SetRange(1,15)
    histoStacked.Draw("hist")

    #overflowBin = histoStacked.GetXaxis().GetLast() + 1
    overflowBin = histoStacked.GetXaxis().GetLast() # This is not the overflow bin
    histoStacked.GetXaxis().SetRange(1,overflowBin)
    c.Update()

    #Superimposing signal
    for filename in config["signal"]:
        c.cd()
        integral = float(inputHistoDic_sig[filename].Integral())
        inputHistoDic_sig[filename].Scale(1./integral)
        if filename.find("mN1p5")>0:
            inputHistoDic_sig[filename].SetLineColor(ROOT.kBlack)
            if filename.find("ctau10p0")>0:
                inputHistoDic_sig[filename].SetLineStyle(1)
            elif filename.find("ctau100p0")>0:
                inputHistoDic_sig[filename].SetLineStyle(2)
            else:
                inputHistoDic_sig[filename].SetLineStyle(3)
        if filename.find("mN1p0")>0:
            inputHistoDic_sig[filename].SetLineColor(ROOT.kBlue)
            if filename.find("ctau10p0")>0:
                inputHistoDic_sig[filename].SetLineStyle(1)
            elif filename.find("ctau100p0")>0:
                inputHistoDic_sig[filename].SetLineStyle(2)
            else:
                inputHistoDic_sig[filename].SetLineStyle(3)
        inputHistoDic_sig[filename].SetLineWidth(2)
        #overflowBin = inputHistoDic_sig[filename].GetXaxis().GetLast() + 1
        overflowBin = inputHistoDic_sig[filename].GetXaxis().GetLast() #This is not the overflow bin
        inputHistoDic_sig[filename].GetXaxis().SetRange(1,overflowBin)
        sigLabel = config["signal"][filename]["label"]
        leg_sigVSbkg.AddEntry(inputHistoDic_sig[filename],sigLabel)
        if plotName.find("2DDist")>0:
            inputHistoDic_sig[filename].GetXaxis().SetRange(1,15)
        inputHistoDic_sig[filename].Draw("hist same")
        inputHistoDic_sig[filename].SetMaximum(yedge)
        iLineStyle+=1

    if plotName.find("2DDist")>0:
        l=ROOT.TLine()
        l.SetLineColor(ROOT.kRed)
        l.SetLineWidth(3)
        l.SetLineStyle(ROOT.kDashed)
        l.DrawLine(1.,0.,1.,yedge)
        l.DrawLine(5.,0.,5.,yedge)
    leg_sigVSbkg.Draw("same")
    c.Update()
  
    #latexpv = ROOT.TLatex()
    #latexpv.SetTextAlign(12)
    #latexpv.SetTextSize(0.04)
    #latexpv.DrawLatexNDC(0.1,0.91,"CMS Private Work")

    outputFullPath = os.path.join(outDirName,"sigVSbkg")
    subprocess.call(["mkdir","-p",outputFullPath])
    out_file_name = os.path.join(outputFullPath,histoName +"_sigVSbkg")
    if args.logy:
        out_file_name += "_logy"
    c.SaveAs(out_file_name + ".png")
    c.SaveAs(out_file_name + ".pdf")
    c.SaveAs(out_file_name + ".root")

    del c
