import subprocess
import sys
import os
import json
script, configFileName = sys.argv

with open(configFileName, "r") as f:
    config = json.loads(f.read())
import ROOT

inputDirName = str(config["inputDirName"])
outDirName = str(config["outDirName"])

for plotName in config["plotNameList"]:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c","c",800,800)
    
    histoName = str(plotName)
    colorList = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kViolet, ROOT.kPink]
    iColor = 0
    iLineStyle=0
    leg_dataVSmc = ROOT.TLegend(0.65, 0.75, 0.87, 0.87)
    leg_sigVSbkg = ROOT.TLegend(0.65, 0.75, 0.87, 0.87)
    
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
    
    #Stacking background histos
    for filename in config["background"]:
        inputFileDic_bkg [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_bkg[filename] = ROOT.TH1D(inputFileDic_bkg[filename].Get(histoName))
        inputHistoDic_bkg[filename].Scale(1./integral_tot)
        inputHistoDic_bkg[filename].SetLineColor(ROOT.kBlack)
        inputHistoDic_bkg[filename].SetFillColor(colorList[iColor])
        #overflowBin = inputHistoDic_bkg[filename].GetXaxis().GetLast() + 1
        overflowBin = inputHistoDic_bkg[filename].GetXaxis().GetLast() # This is not the overflow bin
        inputHistoDic_bkg[filename].GetXaxis().SetRange(1,overflowBin)
        xaxis_label = str(inputHistoDic_bkg[filename].GetXaxis().GetTitle())
        histoStacked.Add(inputHistoDic_bkg[filename])
        iColor += 1
        bkgLabel = config["background"][filename]["label"]
        leg_dataVSmc.AddEntry(inputHistoDic_bkg[filename],bkgLabel)
        leg_sigVSbkg.AddEntry(inputHistoDic_bkg[filename],bkgLabel)

    c.cd()
    histoStacked.Draw() #otherwise next commands befora Draw will fail
    histoStacked.SetTitle("")
    histoStacked.GetYaxis().SetTitle("Normalized to unit")
    histoStacked.SetMinimum(0.)
    histoStacked.GetXaxis().SetTitle(xaxis_label)
    hmax = float(histoStacked.GetMaximum())

    #Superimposing signal
    for filename in config["signal"]:
        inputFileDic_sig [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_sig[filename] = ROOT.TH1D(inputFileDic_sig[filename].Get(histoName))
        c.cd()
        integral = float(inputHistoDic_sig[filename].Integral())
        inputHistoDic_sig[filename].Scale(1./integral)
        inputHistoDic_sig[filename].SetLineColor(ROOT.kBlack)
        inputHistoDic_sig[filename].SetLineWidth(2)
        inputHistoDic_sig[filename].SetLineStyle(1+iLineStyle)
        #overflowBin = inputHistoDic_sig[filename].GetXaxis().GetLast() + 1
        overflowBin = inputHistoDic_sig[filename].GetXaxis().GetLast() #This is not the overflow bin
        inputHistoDic_sig[filename].GetXaxis().SetRange(1,overflowBin)
        sigLabel = config["signal"][filename]["label"]
        hmax_sig = float(inputHistoDic_sig[filename].GetMaximum())
        leg_sigVSbkg.AddEntry(inputHistoDic_sig[filename],sigLabel)
        inputHistoDic_sig[filename].Draw("hist")

        if hmax_sig>hmax:
            hmax = hmax_sig
        yedge = hmax+(hmax/2)
        inputHistoDic_sig[filename].SetMaximum(yedge)
        iLineStyle+=1


    #overflowBin = histoStacked.GetXaxis().GetLast() + 1
    overflowBin = histoStacked.GetXaxis().GetLast() # This is not the overflow bin
    histoStacked.GetXaxis().SetRange(1,overflowBin)
    histoStacked.Draw("hist same")
    c.Update()



    leg_sigVSbkg.Draw("same")
    c.Update()

    outputFullPath = os.path.join(outDirName,"sigVSbkg")
    subprocess.call(["mkdir","-p",outputFullPath])
    c.SaveAs(os.path.join(outputFullPath,histoName +"_sigVSbkg.png"))
    c.SaveAs(os.path.join(outputFullPath,histoName +"_sigVSbkg.root"))

    del c
