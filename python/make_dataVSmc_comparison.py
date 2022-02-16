import subprocess
import sys
import os
import json
script, configFileName = sys.argv

with open(configFileName, "r") as f:
    config = json.loads(f.read())
import ROOT

inputDirName = str(config["inputDirName"])
outDirName = str(os.path.join(config["outDirName"],inputDirName.split("/")[-1])) 

inputMCFileName = inputDirName + "/" + "hadd_bkg.root"
subprocess.call(["hadd","-f",inputMCFileName] + [str(inputDirName + "/" + bkgFileName) for bkgFileName in config["background"].keys()])

inputDataFileName = inputDirName + "/" + "hadd_data.root"
subprocess.call(["hadd","-f",inputDataFileName] + [str(inputDirName + "/" + dataFileName) for dataFileName in config["data"].keys()])

for plotName in config["plotNameList"]:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)

    c = ROOT.TCanvas("c","c",900,800)

    padUpper = ROOT.TPad('padUpper', 'padUpper', 0, 0.3, 1, 1.0)
    padUpper.SetBottomMargin(0.01)
    padUpper.SetTopMargin(0.12)
    padUpper.Draw()
    
    padLower = ROOT.TPad('padLower', 'padLower', 0, 0.0, 1, 0.3)
    padLower.SetBottomMargin(0.35)
    padLower.SetTopMargin(0.12)
    padLower.SetGridy()
    padLower.Draw()

    
    histoName = str(plotName)
    colorList = [ROOT.kRed, ROOT.kGreen, ROOT.kBlue, ROOT.kYellow, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kSpring, ROOT.kTeal, ROOT.kAzure, ROOT.kViolet, ROOT.kPink]
    iColor = 0
    iLineStyle=0
    leg_dataVSmc = ROOT.TLegend(0.65, 0.75, 0.87, 0.87)
    
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
    data_filename = list(config["data"].keys())[0]
    lumi_data = float(config["data"][data_filename]["integrated_lumi"]) #/pb
    
    #Stacking background histos
    for filename in config["background"]:
        inputFileDic_bkg [filename] = ROOT.TFile.Open(str(os.path.join(inputDirName,filename)))
        inputHistoDic_bkg[filename] = ROOT.TH1D(inputFileDic_bkg[filename].Get(histoName))
        inputHistoDic_bkg[filename].Scale(lumi_data)
        inputHistoDic_bkg[filename].SetLineColor(ROOT.kBlack)
        inputHistoDic_bkg[filename].SetFillColor(colorList[iColor])
        overflowBin = inputHistoDic_bkg[filename].GetXaxis().GetLast() + 1
        inputHistoDic_bkg[filename].GetXaxis().SetRange(1,overflowBin)
        xaxis_label = str(inputHistoDic_bkg[filename].GetXaxis().GetTitle())
        histoStacked.Add(inputHistoDic_bkg[filename])
        iColor += 1
        bkgLabel = config["background"][filename]["label"]
        leg_dataVSmc.AddEntry(inputHistoDic_bkg[filename],bkgLabel)

    #Drawing stacked histos
    padUpper.cd()
    histoStacked.Draw("hist")
    histoStacked.SetTitle("")
    histoStacked.GetYaxis().SetTitle("Events")
    histoStacked.GetYaxis().SetLabelSize(0.05)
    histoStacked.GetYaxis().SetTitleSize(0.06)
    histoStacked.GetYaxis().SetTitleOffset(0.8)
    histoStacked.GetXaxis().SetLabelSize(0)
    histoStacked.SetMinimum(0.)
    hmax = float(histoStacked.GetMaximum())
    histoStacked.SetMaximum(hmax+(hmax/3))
    overflowBin = histoStacked.GetXaxis().GetLast() + 1
    histoStacked.GetXaxis().SetRange(1,overflowBin)
    padUpper.Update()
    
    #Superimposing data
    inputDataFile = ROOT.TFile.Open(str(inputDataFileName))
    inputDataHisto = ROOT.TH1D(inputDataFile.Get(histoName))
    padUpper.cd()
    inputDataHisto.SetLineColor(ROOT.kBlack)
    inputDataHisto.SetMarkerStyle(20)
    inputDataHisto.GetXaxis().SetRange(1,overflowBin)
    leg_dataVSmc.AddEntry(inputDataHisto,"DATA")
    inputDataHisto.Draw("ex0p same")
    leg_dataVSmc.Draw("same")

    inputMCFile = ROOT.TFile.Open(str(inputMCFileName))
    inputMCHisto = ROOT.TH1D(inputMCFile.Get(histoName))
    inputMCHisto.Scale(lumi_data)
    
    hRatio = inputDataHisto.Clone()
    hRatio.Divide(inputMCHisto)
    padLower.cd()
    hRatio.SetTitle("")
    hRatio.GetYaxis().SetTitle("")
    hRatio.SetLineColor(ROOT.kBlack)
    hRatio.SetLineWidth(1)
    hRatio.SetMarkerStyle(20)
    hRatio.GetYaxis().SetRangeUser(0.5,1.5)
    hRatio.GetYaxis().SetTitle("DATA/MC")
    hRatio.GetXaxis().SetTitle(inputMCHisto.GetXaxis().GetTitle())
    hRatio.GetXaxis().SetTitleSize(0.12)
    hRatio.GetYaxis().SetTitleSize(0.12)
    hRatio.GetYaxis().SetLabelSize(0.12)
    hRatio.GetXaxis().SetLabelSize(0.12)
    hRatio.GetXaxis().SetRange(-1,overflowBin)
    hRatio.GetYaxis().SetTitleOffset(0.4)
    hRatio.GetYaxis().SetNdivisions(5)
    ROOT.gStyle.SetOptStat(0)
    hRatio.Draw()
    
    subprocess.call(["mkdir","-p",outDirName])
    c.SaveAs(outDirName + "/" + histoName +"_dataVSmc.png")
    c.SaveAs(outDirName + "/" + histoName +"_dataVSmc.root")
    del c
