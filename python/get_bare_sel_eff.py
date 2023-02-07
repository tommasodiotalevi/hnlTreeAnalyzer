import sys
import ROOT
import os
import subprocess
import hnl_tools
import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser(description="")
parser.add_argument("hnl_mass", help="Mass of the HNL [GeV]")
parser.add_argument("hnl_ctau", help="ctau of the HNL [mm]")
parser.add_argument("--reweightToCtau",type=float, default=-1.0, help="Reweight to new ctau, make sure that the variable exists already in the tree")
parser.add_argument("--tag"           ,type=str, default="", help="Take as input files with the specified tag")
parser.add_argument("--nansToZeros",action='store_true', default=False, help="Replace nans with 0.0s")
args = parser.parse_args()

tag = str(args.tag)

#input parameters
m      = float(args.hnl_mass)
ct     = float(args.hnl_ctau)
new_ct = float(args.reweightToCtau)

#categories = ["inclusive"]
categories = ["OSlxy0to1","OSlxy1to5","OSlxy5toInf","SSlxy0to1","SSlxy1to5","SSlxy5toInf"]

out_hnl_mass = list()
out_hnl_ctau = list()
out_hnl_ntot = list()
out_hnl_nsel = list()
out_hnl_eff  = list()
out_hnl_cat  = list()
out_hnl_tag  = list() 

out_ds_ntot = list()
out_ds_nsel = list()
out_ds_eff  = list()
out_ds_cat  = list()
out_ds_tag  = list() 


dsToHnlMu_ntuples_cfg_path = "/home/CMS-T3/lunerti/hnlTreeAnalyzer/cfg/DsToHnlMu_HnlToMuPi_prompt_UL_tree_input_fromCrab.json"
dsToPhiPi_ntuples_cfg_path = "/home/CMS-T3/lunerti/hnlTreeAnalyzer/cfg/DsToPhiPi_PhiToMuMu_prompt_UL_tree_input_fromCrab.json"

with open(dsToHnlMu_ntuples_cfg_path, "r") as f:
    dsToHnlMu_ntuples = json.loads(f.read())

with open(dsToPhiPi_ntuples_cfg_path, "r") as f:
    dsToPhiPi_ntuples = json.loads(f.read())


for cat in categories:
    
    sig_short_name = "DsToNMu_NToMuPi_mN{mass}_ctau{ctau}mm".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"))
    #sig_short_name = "DsToNMu_NToMuPi_mN{mass}_ctau{ctau}mm_incl".format(mass=str(m).replace(".","p"),ctau=str(ct).replace(".","p"))
   

    i_norm_cat = -1 # always the inclusive category
    i_sig_cat  = -1
    final_file_name_list = "final_file_name_list"
    if tag != "":
        final_file_name_list = "final_file_name_list_{}".format(tag)

    for i in range(0,len(dsToPhiPi_ntuples["DsToPhiPi_ToMuMu"][final_file_name_list])):
        s = dsToPhiPi_ntuples["DsToPhiPi_ToMuMu"][final_file_name_list][i]
        if "inclusive" in s:
            i_norm_cat = i
    if i_norm_cat<0:
        print("CATEGORY {} NOT FOUND IN {} CFG FILE".format(cat,dsToPhiPi_ntuples_cfg_path))
        sys.exit(0)

    for i in range(0,len(dsToHnlMu_ntuples[sig_short_name][final_file_name_list])):
        s = dsToHnlMu_ntuples[sig_short_name][final_file_name_list][i]
        if cat in s:
            #print("------> {}".format(s))
            i_sig_cat = i
    if i_sig_cat<0:
        print("CATEGORY {} NOT FOUND IN {} CFG FILE".format(cat,dsToHnlMu_ntuples_cfg_path))
        sys.exit(0)

    #sometimes values are set to nan in csv and then they get skipped
    #to temporarily prevent this I set all 'nan' values to '0.0' so that
    #they don't get skipped
    if args.nansToZeros:
        for fn in [dsToPhiPi_ntuples["DsToPhiPi_ToMuMu"][final_file_name_list][i_norm_cat].replace("root","csv"),dsToHnlMu_ntuples[sig_short_name][final_file_name_list][i_sig_cat].replace("root","csv")]:
            command = "sed -i 's/nan/0.0/g' {}".format(fn)
            print("Running {}".format(command))
            subprocess.call(command,shell=True)

    input_DsPhiPi_ntuple_path = dsToPhiPi_ntuples["DsToPhiPi_ToMuMu"][final_file_name_list][i_norm_cat].replace("root","csv")
    input_DsHnlMu_ntuple_path = dsToHnlMu_ntuples[sig_short_name][final_file_name_list][i_sig_cat].replace("root","csv")
    print("Ds->PhiPi input: {}".format(input_DsPhiPi_ntuple_path))
    print("Ds->HNL input: {}".format(input_DsHnlMu_ntuple_path))

    n_DsToPhiPi_gen = float(dsToPhiPi_ntuples["DsToPhiPi_ToMuMu"]["processed_events"])
    n_DsToPhiPi_sel = float(hnl_tools.get_weighted_yield_from_csv(input_DsPhiPi_ntuple_path,"tot_weight"))
    n_DsToHnlMu_gen = float(dsToHnlMu_ntuples[sig_short_name]["processed_events"])
    n_DsToHnlMu_sel = float(hnl_tools.get_weighted_yield_from_csv(input_DsHnlMu_ntuple_path,"tot_weight"))
    if new_ct>0.:
        ctauwvar = "ctau_weight_{old_ctau}TO{new_ctau}".format(old_ctau=str(ct).replace(".","p"),new_ctau=str(new_ct).replace(".","p"))
        n_DsToHnlMu_sel = float(hnl_tools.get_ctauweighted_yield_from_csv(dsToHnlMu_ntuples[sig_short_name][final_file_name_list][i_sig_cat].replace("root","csv"),ctauwvar,"tot_weight"))
    
    eff_Ds   = n_DsToPhiPi_sel/n_DsToPhiPi_gen
    eff_Hnl  = n_DsToHnlMu_sel/n_DsToHnlMu_gen

    out_hnl_mass.append(float(args.hnl_mass))
    out_hnl_ctau.append(float(args.hnl_ctau))
    out_hnl_ntot.append(n_DsToHnlMu_gen)
    out_hnl_nsel.append(n_DsToHnlMu_sel)
    out_hnl_eff.append(eff_Hnl)
    out_hnl_cat.append(cat)
    out_hnl_tag.append(tag) 

    out_ds_ntot.append(n_DsToPhiPi_gen)
    out_ds_nsel.append(n_DsToPhiPi_sel)
    out_ds_eff.append(eff_Ds)
    out_ds_tag.append(tag) 
    
    
    print("*********CATEGORY*************")
    print("{}".format(cat))
    print("*******INPUT*PARAMETERS*******")
    print("input mass: {} [GeV]".format(m))
    print("input ctau: {} [mm]".format(ct))
    if new_ct>0.:
        print("ctau reweighted to: {} [mm]".format(new_ct))
    print("eff_Ds: {}/{}={}".format(n_DsToPhiPi_sel,n_DsToPhiPi_gen,eff_Ds))
    print("eff_Hnl: {}/{}={}".format(n_DsToHnlMu_sel,n_DsToHnlMu_gen,eff_Hnl))
    print("******************************")
    print('\n')

out_hnl_dict = dict()
out_hnl_dict["m"] = out_hnl_mass
out_hnl_dict["ctau"] = out_hnl_ctau
out_hnl_dict["ntot"] = out_hnl_ntot
out_hnl_dict["nsel"] = out_hnl_nsel
out_hnl_dict["eff"] = out_hnl_eff 
out_hnl_dict["cat"] = out_hnl_cat 
out_hnl_dict["tag"] = out_hnl_tag 

out_ds_dict = dict()
out_ds_dict["ntot"] = out_ds_ntot
out_ds_dict["nsel"] = out_ds_nsel
out_ds_dict["eff"] = out_ds_eff 
out_ds_dict["tag"] = out_ds_tag 

df_hnl = pd.DataFrame(out_hnl_dict)
df_ds  = pd.DataFrame(out_ds_dict)

out_tag = tag
if tag == "":
    out_tag = "noCorr"

df_hnl.to_csv('hnl_mN{}_ctau{}_{}.csv'.format(m,ct,out_tag),index=False) 
df_ds.to_csv('ds_{}.csv'.format(out_tag),index=False) 
