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

#get selection and categorization cuts
with open(config["selection_cfg_file_full_path"], "r") as f:
    selection_cfg = json.loads(f.read())

#get input files
#inputFileName_list = config[dataset_to_process]["file_name_list"]
inputFileName_list = config[dataset_to_process]["file_name_list"][:10]
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
 
n_cuts = 20

input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.find("_tree")]

for cat in selection_cfg["categories"]:
    
    outputDirName = os.path.join(config["output_dir_name"],cat["label"])
    subprocess.call(['mkdir','-p',outputDirName])

    ROOT.EnableImplicitMT()
    df = ROOT.RDataFrame(chain)

    # define a index selecting best candidate in the event
    df = df.Define(selection_cfg["best_cand_var"]["name"],selection_cfg["best_cand_var"]["definition"])

    #build new variables and get bet candidate idx
    for var in selection_cfg["new_variables"]:
        print(var["name"])
        df = df.Define(var["name"],var["definition"])

    # redefine variables so that only best candidate is saved
    for c in df.GetColumnNames():
        col_name = str(c)
        col_type = df.GetColumnType(col_name)
        # choose candidate branches (beginning with 'C_')
        if col_name.find("C_")<0: 
            continue
        # save best candidate only
        if not col_type.find("ROOT::VecOps")<0:
            idx = str(selection_cfg["best_cand_var"]["name"])
            df = df.Redefine(col_name,col_name+"["+idx+"]")
            continue
    
    #define cross section normalization mc_weight
    df = df.Define("mc_weight",str(mc_weight))    

    # define pu weight for MC only
    if dataset_category != "data":
        pu_weight = "h_pu_weights->GetBinContent(h_pu_weights->FindBin(nPU_trueInt))"
    df = df.Define("pu_weight",str(pu_weight)) 

    #define total event weight
    df = df.Define("tot_weight","mc_weight*pu_weight")

    df.Snapshot("test_tree","test_tree.root",df.GetColumnNames())

    #apply categorization
    df = df.Filter(cat["cut"] ,cat["printout"])

    #get mc truth in case of signal sample
    if dataset_category == "signal":
        for sel in selection_cfg["gen_matching_cuts"]:
            df = df.Filter(sel["cut"],sel["printout"])

    tot_events = df.Count().GetValue()

    #build signal selection efficiency csv file 
    #(easier and faster with RDataFrame but it looks like it's not avaliable on lxplus)
    if cat["label"] == "neutralHnl":
        var_header = "entry"
        fmt = "%d"
        df_array = numpy.arange(n_cuts)
        for var in selection_cfg["selection_eff_scan"]:
            a_cuts   = numpy.array(numpy.linspace(var["low_edge"],var["up_edge"],n_cuts))
            pass_events = [df.Filter(str(var["name"])+var["logic"]+str(cut)).Count().GetValue() for cut in a_cuts]
            fail_events = [(tot_events-df.Filter(str(var["name"])+var["logic"]+str(cut)).Count().GetValue()) for cut in a_cuts]
            a_pass = numpy.array(pass_events)
            a_fail = numpy.array(fail_events)
            var_header += ","+var["name"]+"_cut,"+var["name"]+"_pass,"+var["name"]+"_fail"
            fmt += ", %.2f, %.1f,  %.1f"
            df_array = numpy.column_stack((df_array,a_cuts,a_pass,a_fail))
        output_path = os.path.join(outputDirName,"sel_eff_tree_"+dataset_name_label+".csv")
        numpy.savetxt(output_path, df_array, delimiter=',', header=var_header, comments='',fmt=fmt)
        print("Output saved in {}".format(output_path))

    r = df.Report()
    print("-- REPORT --")
    r.Print()

print("--- Analysis completed in {} seconds ---".format(time.time() - start_time))
