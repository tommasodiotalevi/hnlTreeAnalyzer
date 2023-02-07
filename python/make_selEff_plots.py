import ROOT
import sys
import os
import subprocess
import json
import argparse
import math
import importlib.util
from array import array

#import hnl tools from combine dir
spec = importlib.util.spec_from_file_location("hnl_tools", "/afs/cern.ch/work/l/llunerti/private/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/python/hnl_tools.py")
hnl_tools = importlib.util.module_from_spec(spec)
sys.modules["hnl_tools"] = hnl_tools
spec.loader.exec_module(hnl_tools)

parser = argparse.ArgumentParser(description="")
parser.add_argument("cfg_filename", help="Path to the input configuration file")
args = parser.parse_args()

#input parameters
configFileName = args.cfg_filename
m  = float(configFileName[configFileName.find("mN")+2:configFileName.find("mN")+5].replace("p","."))
ct = ctau = float(configFileName[configFileName.find("ctau")+4:configFileName.rfind("p0")+2].replace("p","."))

#open configuration files
with open(configFileName, "r") as f:
    config = json.loads(f.read())

with open(config["dsToHnlMu_ntuples_cfg_file"], "r") as f:
    dsToHnlMu_ntuples = json.loads(f.read())

with open(config["dsToPhiPi_ntuples_cfg_file"], "r") as f:
    dsToPhiPi_ntuples = json.loads(f.read())

inputDirName = str(config["inputDirName"])
outDirName = str(config["outDirName"])
    
sig_short_name = "DsToNMu_NToMuPi_mN{mass}_ctau{ctau}mm".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"))

#sometimes values are set to nan in csv and then they get skipped
#to temporarily prevent this I set all 'nan' values to '0.0' so that
#they don't get skipped
for fn in [config["dsToPhiPi_final_csv"],config["dsToHnlMu_final_csv"]]:
    command = "sed -i 's/nan/0.0/g' {}".format(fn)
    print("Running {}".format(command))
    subprocess.call(command,shell=True)

#get number of generated events and selected events to 
#get the Ds->PhiPi and Ds->HnlMu selection efficiency 
n_DsToPhiPi_gen = float(dsToPhiPi_ntuples["DsToPhiPi_PhiToMuMu"]["processed_events"])
n_DsToPhiPi_sel = float(hnl_tools.get_yield_from_csv(config["dsToPhiPi_final_csv"]))
n_DsToHnlMu_gen = float(dsToHnlMu_ntuples[sig_short_name]["processed_events"])
n_DsToHnlMu_sel = float(hnl_tools.get_yield_from_csv(config["dsToHnlMu_final_csv"]))

print("-------> n_DsToPhiPi_gen: {}".format(n_DsToPhiPi_gen)) 
print("-------> n_DsToPhiPi_sel: {}".format(n_DsToPhiPi_sel)) 
print("-------> n_DsToHnlMu_gen: {}".format(n_DsToHnlMu_gen)) 
print("-------> n_DsToHnlMu_sel: {}".format(n_DsToHnlMu_sel)) 

eff_Ds   = n_DsToPhiPi_sel/n_DsToPhiPi_gen
eff_Hnl  = n_DsToHnlMu_sel/n_DsToHnlMu_gen
f_prompt = config["fraction_prompt_ds"]
n_Ds     = 2110

#get signal and background yield before selection optimization
nS = hnl_tools.get_expected_signal_yield(m,ct,n_Ds,f_prompt,eff_Hnl,eff_Ds)
nB = hnl_tools.get_yield_from_workspace(config["background"]["workspace_filename"],config["background"]["workspace"],config["background"]["yield_var"])

def significance(s,b):
    sig = 0
    if b!=0:
        sig = math.sqrt(2*( (s+b)*math.log(1+(s/b)) - s))
    return sig

def significance_error(s,ps,fs,pb,fb,nS,nB):
    e = s*math.sqrt(protected_division(1.,ps)**2 +
                    protected_division(1.,fs)**2 +
                    protected_division(1.,pb)**2 +
                    protected_division(1.,fb)**2 +
                    protected_division(1.,nS)**2 +
                    protected_division(1.,nB)**2)
    return e

def protected_division(n,d):
    output = 0.
    if d!= 0:
        output = n/d
    return output


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    bkg_fileName = os.path.join(config["inputDirName"],config["background"]["filename"])
    sig_fileName = os.path.join(config["inputDirName"],config["signal"]["filename"])
    print("Input file for background: {}".format(bkg_fileName))
    print("Input file for signal: {}".format(sig_fileName))
    cat = config["background"]["filename"].split(".")[0].split("_")[-1]
    
    print("*********CATEGORY*************")
    print("{}".format(cat))
    print("*******INPUT*PARAMETERS*******")
    print("input mass: {} [GeV]".format(m))
    print("input ctau: {} [mm]".format(ct))
    print("n_Ds = {}".format(n_Ds))
    print("f_prompt = {}".format(f_prompt))
    print("eff_Ds = {}".format(eff_Ds))
    print("eff_Hnl = {}".format(eff_Hnl))
    print("signal yield before optimization = {}".format(nS))
    print("background yield before optimization = {}".format(nB))
    print("******************************")
    print('\n')

    st = ROOT.TTree()
    bt = ROOT.TTree()
    bt.ReadFile(bkg_fileName)
    st.ReadFile(sig_fileName)
    sdf = ROOT.RDataFrame(st)
    bdf = ROOT.RDataFrame(bt)
    for var in config["variables"]:
        c = ROOT.TCanvas("c","c",800,800)
        x = var["name"]
        sdf = sdf.Define(x+"_eff",str(x+"_pass/"+"("+x+"_pass+"+x+"_fail)"))
        bdf = bdf.Define(x+"_eff",str(x+"_pass/"+"("+x+"_pass+"+x+"_fail)"))
        eff_s = [a for a in sdf.Take['float'](x+'_eff').GetValue()]
        eff_b = [a for a in bdf.Take['float'](x+'_eff').GetValue()]
        n_pass_s = [a for a in sdf.Take['float'](x+'_pass').GetValue()]
        n_fail_s = [a for a in sdf.Take['float'](x+'_fail').GetValue()]
        n_pass_b = [a for a in bdf.Take['float'](x+'_pass').GetValue()]
        n_fail_b = [a for a in bdf.Take['float'](x+'_fail').GetValue()]
        sign = array('f')
        for (s,b) in zip(eff_s,eff_b):
            sign.append(significance(nS*s,nB*b))
        cuts = array('f')
        for a in bdf.Take['float'](x+'_cut').GetValue():
            cuts.append(a)
        sign_err = array('f')
        for (ps,fs,pb,fb,s) in zip(n_pass_s,n_fail_s,n_pass_b,n_fail_b,sign):
            e = significance_error(s,ps,fs,pb,fb,nS,nB)
            sign_err.append(e)
        cuts_err = array('f')
        for i in range(0,len(cuts)):
            cuts_err.append(0.)
        n = len(cuts)
        #gr = ROOT.TGraph(n,cuts,sign)
        gr = ROOT.TGraphErrors(n,cuts,sign,cuts_err,sign_err)
        gr.SetTitle("")
        gr.Draw()
        gr.SetMarkerStyle(24)
        gr.GetXaxis().SetTitle(var["x_label"])
        gr.GetYaxis().SetTitle("Significance")

        lcat = ROOT.TLatex()
        lcat.SetTextSize(0.03)
        lcat.SetTextAlign(11)
        lcat.SetTextFont(42)
        lcat.DrawLatexNDC(.2,.925,"#bf{#bf{"+cat+" category (m="+str(m)+" GeV,"+"c#tau="+str(ct)+" mm)}}")

        c.SaveAs(os.path.join(config["outDirName"],x+"_mN{mass}_ctau{ctau}_{cat}.png".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"),cat=cat)))
        c.SaveAs(os.path.join(config["outDirName"],x+"_mN{mass}_ctau{ctau}_{cat}.pdf".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"),cat=cat)))
        gr.SaveAs(os.path.join(config["outDirName"],x+"_mN{mass}_ctau{ctau}_{cat}.root".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"),cat=cat)))
        del c



