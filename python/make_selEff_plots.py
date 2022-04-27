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

            sel_eff = varName+"_pass/("+varName+"_pass+"+varName+"_fail)"
            rej_eff = "1.0-("+varName+"_pass/("+varName+"_pass+"+varName+"_fail))"
            df = df.Define("sel_eff",sel_eff)
            df = df.Define("rej_eff",rej_eff)
            var_to_plot = "sel_eff"
            # plot 1-selection_eff=bkg_rejection in case of background
            if cat == "background":
                var_to_plot = "rej_eff"

            gr_dict[cat][i] = df.Graph(varName+"_cut",var_to_plot)
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
