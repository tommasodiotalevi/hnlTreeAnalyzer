import ROOT
import sys
import os
import subprocess
import json
script, configFileName = sys.argv

with open(configFileName, "r") as f:
    config = json.loads(f.read())

inputDirName = str(config["inputDirName"])
outDirName = str(config["outDirName"])

for varName in config["variables"]:
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("c","c",900,800)

    #initialize TGraphs here otherwise it will save empty canvas
    gr_dict = {"signal" : [ROOT.TGraph() for fileName in config["signal"]],
               "background" : [ROOT.TGraph() for fileName in config["background"]]}

    leg = ROOT.TLegend(0.60, 0.55, 0.87, 0.67)

    k = int(0)
    for cat in ["signal","background"]:
        for i in range(0,len(config[cat])):
            fileNameList = list(config[cat].keys())
            fileName = fileNameList[i]
            t = ROOT.TTree()
            t.ReadFile(os.path.join(inputDirName,fileName))
            df = ROOT.RDataFrame(t)
            gr_dict[cat][i] = df.Graph(varName+"_cut",varName+"_eff")
            gr_dict[cat][i].GetXaxis().SetTitle(config["variables"][varName]["x_label"])
            gr_dict[cat][i].GetYaxis().SetTitle("")
            gr_dict[cat][i].GetYaxis().SetRangeUser(-0.1,1.1)
            gr_dict[cat][i].SetTitle("")
            gr_dict[cat][i].SetMarkerStyle(config[cat][fileName]["marker"])
            gr_dict[cat][i].SetMarkerSize(1.5)
            opt = "ap"
            if k>0:
                opt ="p same"

            gr_dict[cat][i].Draw(opt)

            leg_label = str(config[cat][fileName]["label"])
            leg.AddEntry(gr_dict[cat][i].GetValue(),leg_label,"p")

            k+=1

    leg.Draw("same")
    
    subprocess.call(["mkdir","-p",outDirName])
    c.SaveAs(outDirName + "/" + varName +"_selEff.png")
    c.SaveAs(outDirName + "/" + varName +"_selEff.root")
    del c
