import ROOT
import sys
import math
import os
import json
import time
import subprocess
import numpy
start_time = time.time()

# let private function to be avaliable in RDataFrame evn
header_path = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/interface/df_tools.h"
ROOT.gInterpreter.Declare('#include "{}"'.format(header_path))

configFileName = sys.argv[1]
dataset_to_process = sys.argv[2]

#open analyzer configuration file
with open(configFileName, "r") as f:
    config = json.loads(f.read())

#get histogram configuration
with open(config["histogram_cfg_file_full_path"], "r") as f:
    histos = json.loads(f.read())

#get selection and categorization cuts
with open(config["selection_cfg_file_full_path"], "r") as f:
    selection_cfg = json.loads(f.read())

#get input files
inputFileName_list = config[dataset_to_process]["file_name_list"]
#inputFileName_list = config[dataset_to_process]["file_name_list"][:10]
dataset_category = config[dataset_to_process]["dataset_category"]

#get pu weights histogram
ROOT.gInterpreter.Declare("""
auto pu_weights_file = TFile::Open("/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/pu_weights.root");
auto h_pu_weights = pu_weights_file->Get<TH1D>("pu_weights");
""")

#define unit weights for data
mc_weight = 1.
pu_weight = 1.

if dataset_category != "data":
    cross_section      = float(config[dataset_to_process]["cross_section"])
    filter_efficiency  = float(config[dataset_to_process]["filter_efficiency"])
    total_events       = float(config[dataset_to_process]["processed_events"])
    mc_weight = cross_section*filter_efficiency/total_events

print("mc_weight: {}".format(mc_weight))

#chain = ROOT.TChain('slimmed_tree')
chain = ROOT.TChain('wztree')
for inputFileName in inputFileName_list:
    print("Adding {} to the chain...".format(inputFileName))
    success = chain.Add(inputFileName)
    print("exit: {}".format(success))

print("\n")
print("{} total entries ...".format(chain.GetEntries()))

input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.find("_tree")]

for cat in selection_cfg["categories"]:
    
    outputFileName = "out_tree_analyzer_"+dataset_name_label+"_"+cat["label"]+".root"
    outputDirName = os.path.join(config["output_dir_name"],cat["label"])
    
    subprocess.call(['mkdir','-p',outputDirName])
    outputFile = ROOT.TFile(os.path.join(outputDirName,outputFileName),"RECREATE")

    ROOT.EnableImplicitMT()
    df = ROOT.RDataFrame(chain)

    # define a index selecting best candidate in the event
    df = df.Define(selection_cfg["best_cand_var"]["name"],selection_cfg["best_cand_var"]["definition"])

    #build new variables and get bet candidate idx
    for var in selection_cfg["new_variables"]:
        df = df.Define(var["name"],var["definition"])

    # redefine variables so that only best candidate is saved
    for c in df.GetColumnNames():
        col_name = str(c)
        col_type = df.GetColumnType(col_name)
        # choose candidate branches (beginning with 'C_')
        if col_name.find("C_")<0: 
        #if col_name.find("Bs_")<0 and col_name.find("Ds_")<0 and col_name.find("Phi_")<0 and col_name.find("mu1_Phi_")<0 and col_name.find("mu2_Phi_")<0 and col_name.find("mu1_B_")<0 and col_name.find("pi_")<0 and col_name.find("C_")<0:
            continue
        # save best candidate only
        if not col_type.find("ROOT::VecOps")<0:
            idx = str(selection_cfg["best_cand_var"]["name"])
            df = df.Redefine(col_name,col_name+"["+idx+"]")
            continue
    
    #define cross section normalization mc_weight
    df = df.Define("mc_weight",str(mc_weight))    

    ## define pu weight for MC only
    #if dataset_category != "data":
    #    pu_weight = "h_pu_weights->GetBinContent(h_pu_weights->FindBin(nPU_trueInt))"
    df = df.Define("pu_weight",str(pu_weight)) 

    #define total event weight
    df = df.Define("tot_weight","mc_weight*pu_weight")


    #apply categorization
    df = df.Filter(cat["cut"] ,cat["printout"])

    #apply selection
    for sel in selection_cfg["selection_cuts"]:
        df = df.Filter(sel["cut"],sel["printout"])
    
    #get mc truth in case of signal sample
    if dataset_category == "signal":
        for sel in selection_cfg["gen_matching_cuts"]:
            df = df.Filter(sel["cut"],sel["printout"])

    df.Snapshot("test_tree","test_tree.root",df.GetColumnNames())
    

    #book histograms
    histo_dict = {}
    for histo_name in histos:
        title = str(histos[histo_name]["title"])
        nbins = int(histos[histo_name]["nbins"])
        xlow  = float(histos[histo_name]["xlow"])
        xhigh = float(histos[histo_name]["xhigh"])
        histo_model = (histo_name,title,nbins,xlow,xhigh)

        var_name = str(histos[histo_name]["var"])
        #print(var_name)

        histo_dict[histo_name]= df.Histo1D(histo_model,var_name,"tot_weight")

    #save histograms
    for histo_name in histo_dict:
        histo_dict[histo_name].Write()
    
    outputFile.Close()
    
    r = df.Report()
    print("-- REPORT --")
    r.Print()
    print("Output saved in {}".format(os.path.join(outputDirName,outputFileName)))

print("--- Analysis completed in {} seconds ---".format(time.time() - start_time))
