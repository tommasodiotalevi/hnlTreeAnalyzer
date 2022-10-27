import ROOT
import sys
import math
import os
import json
import time
import subprocess
import numpy
import argparse
start_time = time.time()

#script input arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("cfg_filename"         ,                                     help="Path to the input configuration file")
parser.add_argument("--datasetsToProcess"  , nargs='*'          , default=[]   , help="List of short names of the input sample to process")
parser.add_argument("--nThreads"           , type=int           , default=1    , help="Number of threads")
parser.add_argument("--nCuts"              , type=int           , default=20   , help="Number of cuts for the efficiency scan")
args = parser.parse_args()

ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(args.nThreads)

configFileName      = args.cfg_filename
datasets_to_process = args.datasetsToProcess

# let private function to be avaliable in RDataFrame evn
header_path = "/home/CMS-T3/lunerti/hnlTreeAnalyzer/interface/df_tools.h"
ROOT.gInterpreter.Declare('#include "{}"'.format(header_path))

#open analyzer configuration file
with open(configFileName, "r") as f:
    config = json.loads(f.read())

#get ntuples configuration
with open(config["ntuples_cfg_file_full_path"], "r") as f:
    ntuples = json.loads(f.read())

#get selection and categorization cuts
with open(config["selection_cfg_file_full_path"], "r") as f:
    selection = json.loads(f.read())

input_file_list = "final_file_name_list"
tree_name = "final_tree"

#get input files
dataset_category = ntuples[datasets_to_process[0]]["dataset_category"]
inputFileName_list = []
for dataset_to_process in datasets_to_process:
    inputFileName_list += ntuples[dataset_to_process][input_file_list]
    if ntuples[dataset_to_process]["dataset_category"] != dataset_category:
        print("Please provide input samples of the same category (e.g. all background samples)")
        exit(0)

print("----> {}".format(inputFileName_list))
 
n_cuts = args.nCuts
reports = {}


#build a dictionary that groups all the files
#of the same catogory togheter

input_file_catOrdered = {}
for cat in selection["categories"]:
    file_list = []
    for inputFileName in inputFileName_list:
        file_cat_label = inputFileName.split("/")[-1].split(".")[0].split("_")[-1]
        if file_cat_label == cat["label"]:
            file_list.append(inputFileName)
    input_file_catOrdered[cat["label"]] = file_list

print(input_file_catOrdered)


for cat in selection["categories"]:

    chain = ROOT.TChain(tree_name)
    
    for inputFileName in input_file_catOrdered[cat["label"]]:
        print("Adding {} to the chain...".format(inputFileName))
        success = chain.Add(inputFileName)

    print("\n")
    print("{} total entries ...".format(chain.GetEntries()))
    
    outputDirName = os.path.join(config["output_dir_name"])
    subprocess.call(['mkdir','-p',outputDirName])

    df = ROOT.RDataFrame(chain)

    #get mc truth in case of signal sample
    if dataset_category == "signal":
        for sel in selection["gen_matching_cuts"]:
            df = df.Filter(sel["cut"],sel["printout"])

    tot_events = df.Count().GetValue()

    #build signal selection efficiency csv file 
    var_header = "entry"
    fmt = "%d"
    df_array = numpy.arange(n_cuts)
    for var in selection["selection_eff_scan"]:
        a_cuts   = numpy.array(numpy.linspace(var["low_edge"],var["up_edge"],n_cuts))
        pass_events = [df.Filter(str(var["name"])+var["logic"]+str(cut)).Count().GetValue() for cut in a_cuts]
        fail_events = [(tot_events-df.Filter(str(var["name"])+var["logic"]+str(cut)).Count().GetValue()) for cut in a_cuts]
        a_pass = numpy.array(pass_events)
        a_fail = numpy.array(fail_events)
        var_header += ","+var["name"]+"_cut,"+var["name"]+"_pass,"+var["name"]+"_fail"
        fmt += ", %.2f, %.1f, %.1f"
        df_array = numpy.column_stack((df_array,a_cuts,a_pass,a_fail))
    output_path = os.path.join(outputDirName,"sel_eff_tree_"+dataset_category+"_"+cat["label"]+".csv")
    numpy.savetxt(output_path, df_array, delimiter=',', header=var_header, comments='',fmt=fmt)
    print("Output saved in {}".format(output_path))

    #reports[cat["label"]] = df.Report()

#print("+++ FINAL REPORT +++")
#for c in reports:
#    print("--> {} category".format(c))
#    reports[c].Print()

print("--- Analysis completed in {} seconds ---".format(time.time() - start_time))

