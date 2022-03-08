import ROOT
import sys
import math
import os
import json
import time
import subprocess
import numpy
start_time = time.time()

configFileName = sys.argv[1]
dataset_to_process = sys.argv[2]

with open(configFileName, "r") as f:
    config = json.loads(f.read())

histoConfigFileName = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/cfg/histos.json"
with open(histoConfigFileName, "r") as f:
    histos = json.loads(f.read())

#get input files
inputFileName_list = config[dataset_to_process]["slimmed_file_name_list"]
dataset_category = config[dataset_to_process]["dataset_category"]

event_weight = 1.

if dataset_category != "data":
    cross_section      = float(config[dataset_to_process]["cross_section"])
    filter_efficiency  = float(config[dataset_to_process]["filter_efficiency"])
    total_events       = float(config[dataset_to_process]["processed_events"])
    event_weight = cross_section*filter_efficiency/total_events

print("event_weight: {}".format(event_weight))

chain = ROOT.TChain('slimmed_tree')
for inputFileName in inputFileName_list:
    print("Adding {} to the chain...".format(inputFileName))
    chain.Add(inputFileName)

cat_dict = {"chargedHnl":"C_mu1_charge_ptBest+C_pi_charge_ptBest != 0",
            "neutralHnl":"C_mu1_charge_ptBest+C_pi_charge_ptBest == 0"}

n_cuts = 20
sel_cuts = {"C_Hnl_lxy_ptBest":numpy.linspace(-1,30,n_cuts),
            "C_Hnl_pt_ptBest": numpy.linspace(-1,30,n_cuts)}

n_threads = 8

for cat in cat_dict:
    input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
    dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.find("_tree")]
    
    outputFileName = "out_hnl_tree_analyzer_"+dataset_name_label+"_"+cat+".root"
    outputDirName = os.path.join(config["output_dir_name"],cat)
    
    subprocess.call(['mkdir','-p',outputDirName])
    outputFile = ROOT.TFile(os.path.join(outputDirName,outputFileName),"RECREATE")
    
    ROOT.EnableImplicitMT(n_threads)
    df = ROOT.RDataFrame(chain)
    
    df = df.Define("weight",str(event_weight))    
    df = df.Filter(cat_dict[cat] ,cat+" category")
    df = df.Filter("mu9_ip6_matched_ptBest>0" ,"mu9_ip6 matched")
    tot_events = df.Count().GetValue()

    #build signal selection efficiency csv file 
    #(easier and faster with RDataFrame but it looks like it's not avaliable on lxplus)
    if cat == "neutralHnl":
        var_header = "entry"
        fmt = "%d"
        df_array = numpy.arange(n_cuts)
        for var in sel_cuts:
            a_cuts   = numpy.array(sel_cuts[var])
            events = [df.Filter(str(var)+">"+str(cut)).Count().GetValue() for cut in sel_cuts[var]]
            eff = []
            if dataset_category == "signal":
                eff = [x/tot_events for x in events] #signal efficiency
            else:
                eff = [1-(x/tot_events) for x in events] #background rejection
            a_eff = numpy.array(eff)
            var_header += ","+var+"_cut,"+var+"_eff"
            fmt += ", %.2f, %.10f"
            df_array = numpy.column_stack((df_array,a_cuts,eff))

        numpy.savetxt(os.path.join(outputDirName,"sel_eff_tree_"+dataset_name_label+".csv"), df_array, delimiter=',', header=var_header, comments='',fmt=fmt)
    
    #get mc truth in case of signal sample
    if dataset_category == "signal":
        df = df.Filter("C_mu1_isHnlDaughter_ptBest>0 && C_pi_isHnlDaughter_ptBest>0 && C_mu2_isHnlBrother_ptBest>0")

    #book histograms
    histo_dict = {}
    for histo_name in histos:
        title = str(histos[histo_name]["title"])
        nbins = int(histos[histo_name]["nbins"])
        xlow  = float(histos[histo_name]["xlow"])
        xhigh = float(histos[histo_name]["xhigh"])
        histo_model = (histo_name,title,nbins,xlow,xhigh)

        var_name = str(histos[histo_name]["var"])

        histo_dict[histo_name]= df.Histo1D(histo_model,var_name,"weight")

    #save histograms
    for histo_name in histo_dict:
        histo_dict[histo_name].Write()
    
    outputFile.Close()
    
    r = df.Report()
    print("-- REPORT --")
    r.Print()

print("Output saved in {}".format(os.path.join(outputDirName,outputFileName)))
print("---{} threads, {} seconds ---".format(n_threads,time.time() - start_time))
