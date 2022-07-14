import ROOT
import sys
import os
import subprocess
import json
import argparse
import math
from array import array

parser = argparse.ArgumentParser(description="")
parser.add_argument("cfg_filename", help="Path to the input configuration file")
parser.add_argument("--nS"        , type=int           , default=2435  , help="Number of expected signal events")
parser.add_argument("--nB"        , type=int           , default=21542 , help="Number of expected background events")
args = parser.parse_args()

configFileName = args.cfg_filename
nB = args.nB
nS = args.nS

with open(configFileName, "r") as f:
    config = json.loads(f.read())

inputDirName = str(config["inputDirName"])
outDirName = str(config["outDirName"])

def significance(s,b):
    sig = 0
    if b!=0:
        sig = math.sqrt(2*( (s+b)*math.log(1+(s/b)) - s))
    return sig


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    print("Input file for background: {}".format(os.path.join(config["inputDirName"],config["background"]["filename"])))
    print("Input file for signal: {}".format(os.path.join(config["inputDirName"],config["signal"]["filename"])))

    st = ROOT.TTree()
    bt = ROOT.TTree()
    st.ReadFile(os.path.join(config["inputDirName"],config["background"]["filename"]))
    bt.ReadFile(os.path.join(config["inputDirName"],config["signal"]["filename"]))
    sdf = ROOT.RDataFrame(st)
    bdf = ROOT.RDataFrame(bt)
    for var in config["variables"]:
        c = ROOT.TCanvas("c","c",800,800)
        x = var["name"]
        sdf = sdf.Define(x+"_eff",x+"_pass/"+"("+x+"_pass+"+x+"_fail)")
        bdf = bdf.Define(x+"_eff",x+"_pass/"+"("+x+"_pass+"+x+"_fail)")
        eff_s = [a for a in sdf.Take['float'](x+'_eff').GetValue()]
        eff_b = [a for a in bdf.Take['float'](x+'_eff').GetValue()]
        sign = array('f')
        for (s,b) in zip(eff_s,eff_b):
            sign.append(significance(nS*s,nB*b))
        cuts = array('f')
        for a in bdf.Take['float'](x+'_cut').GetValue():
            cuts.append(a)
        n = len(cuts)
        gr = ROOT.TGraph(n,cuts,sign)
        gr.Draw()
        gr.SetMarkerStyle(24)
        gr.GetXaxis().SetTitle(var["x_label"])
        gr.GetYaxis().SetTitle("Significance")
        c.SaveAs(os.path.join(config["outDirName"],x+".png"))
        del c



