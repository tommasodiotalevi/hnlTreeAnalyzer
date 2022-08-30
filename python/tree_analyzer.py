import ROOT
import sys
import math
import os
import json
import time
import subprocess
import argparse
import numpy
start_time = time.time()

#script input arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("cfg_filename"      ,                                     help="Path to the input configuration file")
parser.add_argument("dataset_short_name",                                     help="Short name of the input sample to process")
parser.add_argument("--saveOutputTree"  , action='store_true', default=False, help="Save an output tree (after all the cuts)")
parser.add_argument("--saveSlimmedTree" , action='store_true', default=False, help="Save an output tree (before the cuts and after best candidate selection)")
parser.add_argument("--noHistograms"    , action='store_true', default=False, help="Do not save histogram as output")
parser.add_argument("--addSPlotWeight"  , action='store_true', default=False, help="Add splot weight for data")
parser.add_argument("--skipSlimCuts"    , action='store_true', default=False, help="Skip best candidate selection (N.B. a tree where the best candidate has been already been selected has to be provided in the input cfg file)")
parser.add_argument("--skipSelCuts"     , action='store_true', default=False, help="Do not apply selection cuts")
parser.add_argument("--skipPUrw"        , action='store_true', default=False, help="Do not apply PU reweighting")
parser.add_argument("--nThreads"        , type=int           , default=1    , help="Number of threads")
parser.add_argument("--ctauReweighting" , action='store_true', default=False, help="Include ctau reweighting")
args = parser.parse_args()

ROOT.EnableThreadSafety()
ROOT.EnableImplicitMT(args.nThreads)

configFileName     = args.cfg_filename
dataset_to_process = args.dataset_short_name

#open analyzer configuration file
with open(configFileName, "r") as f:
    config = json.loads(f.read())

#get ntuples configuration
with open(config["ntuples_cfg_file_full_path"], "r") as f:
    ntuples = json.loads(f.read())

#get histogram configuration
with open(config["histogram_cfg_file_full_path"], "r") as f:
    histos = json.loads(f.read())

#get selection and categorization cuts
with open(config["selection_cfg_file_full_path"], "r") as f:
    selection = json.loads(f.read())

# let private function to be avaliable in RDataFrame evn
header_path = config["user_defined_function_path"]
ROOT.gInterpreter.Declare('#include "{}"'.format(header_path))

input_file_list = "file_name_list"
tree_name = "wztree"

if args.skipSlimCuts:
    input_file_list = "slimmed_file_name_list"
    tree_name = "slimmed_tree"

if args.skipSelCuts and args.skipSlimCuts:
    input_file_list = "final_file_name_list"
    tree_name = "final_tree"

#get input files
inputFileName_list = ntuples[dataset_to_process][input_file_list]
dataset_category = ntuples[dataset_to_process]["dataset_category"]

#get pu weights histogram
if dataset_category != "data" and not args.skipPUrw:
    ROOT.gInterpreter.Declare("""
    auto pu_weights_file = TFile::Open("{}");
    auto h_pu_weights = pu_weights_file->Get<TH1D>("{}");
    """.format(config["pu_weight_input_file"],config["pu_weight_histo_name"])
    )

#define unit weights
mc_weight = 1.
pu_weight = 1.

# get generator weight for MC
if dataset_category != "data":
    cross_section      = float(ntuples[dataset_to_process]["cross_section"])
    filter_efficiency  = float(ntuples[dataset_to_process]["filter_efficiency"])
    total_events       = float(ntuples[dataset_to_process]["processed_events"])
    mc_weight = cross_section*filter_efficiency/total_events

print("mc_weight: {}".format(mc_weight))

#initialize chain
chain = ROOT.TChain(tree_name)
tot_file = len(inputFileName_list)
files_in_the_chain = int(0)

#add file to the chain
for inputFileName in inputFileName_list:
    print("Adding {} to the chain...".format(inputFileName))
    chain.Add(inputFileName)
    files_in_the_chain += 1
    print("[{}/{}] files added to the chain".format(files_in_the_chain,tot_file))
    print("{} entries now in the chain".format(chain.GetEntries()))

print("\n")
print("{} total entries ...".format(chain.GetEntries()))

input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.rfind("_")]

# save slimmed tree only once without applying categorization
slimmed_tree_has_been_saved = False
if args.saveOutputTree:
    ntuples[dataset_to_process]["final_file_name_list"] = []

reports = {}
weighted_events_reports = {}
for cat in selection["categories"]:

    #initialize data frame
    df = ROOT.RDataFrame(chain)
    
    if not args.skipSlimCuts:
        # define a index selecting best candidate in the event
        df = df.Define(selection["best_cand_var"]["name"],selection["best_cand_var"]["definition"])
          
        #build new variables and get bet candidate idx
        for var in selection["new_variables"]:
            df = df.Define(var["name"],var["definition"])

    
        # redefine variables so that only best candidate is saved
        for c in df.GetColumnNames():
            col_name = str(c)
            candidate_columns = []
            col_type = df.GetColumnType(col_name)
            # choose candidate branches (beginning with 'C_')
            if col_name.find("C_")<0: 
                continue
            # save best candidate only
            if not col_type.find("ROOT::VecOps")<0:
                idx = str(selection["best_cand_var"]["name"])
                df = df.Redefine(col_name,col_name+"["+idx+"]")
                continue
        

        #define cross section normalization mc_weight
        df = df.Define("mc_weight",str(mc_weight))    
    
        # define pu weight for MC only
        if dataset_category != "data" and not args.skipPUrw:
            pu_weight = "h_pu_weights->GetBinContent(h_pu_weights->FindBin(nPU_trueInt))"
        df = df.Define("pu_weight",str(pu_weight)) 
    
        df = df.Define("tot_weight","mc_weight*pu_weight")

        # save slimmed tree: only the best candidate is saved for each event
        if args.saveSlimmedTree and not slimmed_tree_has_been_saved:
            slim_outputFileName = "slimmed_"+dataset_name_label+".root"
            slim_outputDirName = os.path.join(config["slimmed_tree_output_dir_name"],dataset_name_label)
            subprocess.call(['mkdir','-p',slim_outputDirName])
            slim_outFullPath = os.path.join(slim_outputDirName,slim_outputFileName)
            df.Snapshot(config["slimmed_tree_output_name"],slim_outFullPath,df.GetColumnNames())
            slimmed_tree_has_been_saved = True
            print("Slimmed tree saved in {}".format(slim_outFullPath))
            ntuples[dataset_to_process]["slimmed_file_name_list"] = [str(slim_outFullPath)]
            with open(config["ntuples_cfg_file_full_path"], "w") as f:
                json.dump(ntuples,f, indent=4, sort_keys=True)
            print("{} updated".format(config["ntuples_cfg_file_full_path"])) 

    #apply categorization 
    df = df.Filter(cat["cut"] ,cat["printout"])
    
    if not args.skipSelCuts:
        #apply selection
        for sel in cat["selection_cuts"]:
            df = df.Filter(sel["cut"],sel["printout"])

        
        #get mc truth in case of signal sample
        if dataset_category == "signal":
            for sel in selection["gen_matching_cuts"]:
                df = df.Filter(sel["cut"],sel["printout"])
            if args.ctauReweighting and dataset_category == "signal":
                old_ctau_label = dataset_name_label[dataset_name_label.find("ctau")+4:dataset_name_label.find("mm")]
                hnl_mass_label = dataset_name_label[dataset_name_label.find("mN")+2:dataset_name_label.find("mN")+5]
                for new_ctau in selection["mN"+hnl_mass_label+"_ctau"+old_ctau_label+"mm_rw_points"]:
                  old_ctau = float(old_ctau_label.replace("p","."))
                  w_expr   = "("+str(old_ctau)+"/"+str(new_ctau)+")*"+"exp(C_Hnl_gen_l_prop*("+str(1./old_ctau)+"-"+str(1./new_ctau)+"))"
                  #print("---> weight = {}".format(w_expr))
                  df = df.Define("ctau_weight_"+old_ctau_label+"TO"+str(new_ctau).replace(".","p"),w_expr)

        if args.saveOutputTree and cat["save"]=="yes":
            finalTree_outputFileName = "tree_"+dataset_name_label+"_"+cat["label"]+".root"
            finalCSV_outputFileName  = "tree_"+dataset_name_label+"_"+cat["label"]+".csv"
            finalTree_outputDirName = os.path.join(config["tree_output_dir_name"],dataset_name_label)
            subprocess.call(['mkdir','-p',finalTree_outputDirName])
            finalTree_outFullPath = os.path.join(finalTree_outputDirName,finalTree_outputFileName)
            finalCSV_outFullPath = os.path.join(finalTree_outputDirName,finalCSV_outputFileName)
            #save output tree
            df.Snapshot(config["tree_output_name"],finalTree_outFullPath,df.GetColumnNames())
            #save output csv
            a = df.AsNumpy([x for x in df.GetColumnNames() if x.find("C_")==0 or x.find("ctau_weight")==0])
            arr = numpy.array([x for x in a.values()]).transpose()
            #arr = numpy.array([numpy.nan_to_num(x) for x in a.values()]).transpose()
            numpy.savetxt(finalCSV_outFullPath, arr, delimiter=',', header=",".join([str(x) for x in a.keys()]), comments='')
            #output_tree_has_been_saved = True
            print("Output tree saved in {}".format(finalTree_outFullPath))
            print("Output csv saved in {}".format(finalCSV_outFullPath))
            ntuples[dataset_to_process]["final_file_name_list"] += [str(finalTree_outFullPath)]
            with open(config["ntuples_cfg_file_full_path"], "w") as f:
                json.dump(ntuples,f, indent=4, sort_keys=True)
            print("{} updated".format(config["ntuples_cfg_file_full_path"]))

    if dataset_category=="data" and args.addSPlotWeight:
        sdf = ROOT.RDataFrame(config["splot_weight_tree_name"],config["splot_weight_input_file"]) # get tree containing splot weights
        asw = sdf.AsNumpy([config["splot_weight_variable"]])[config["splot_weight_variable"]] # get column of splot weights
        print("entries in splot tree: {}".format(len(asw)))
        print("entries in analyzed tree: {}".format(df.Count().GetValue()))
        if len(asw) != df.Count().GetValue():
            print("!!! splot tree and analyzed tree do not have the same number of entries !!!")
            print("!!! please check that you are using the correct input tree !!!")
            exit()
        df = df.Define("splot_weight",'auto to_eval = std::string("asw[") + std::to_string(rdfentry_) + "]"; return float(TPython::Eval(to_eval.c_str()));') 
        df = df.Redefine("tot_weight","tot_weight*splot_weight")

    #save histograms
    if not args.noHistograms:
        histo_outputFileName = "histograms_"+dataset_name_label+"_"+cat["label"]+".root"
        #histo_outputDirName = os.path.join(config["output_dir_name"],cat["label"])        
        histo_outputDirName = config["output_dir_name"]        
        subprocess.call(['mkdir','-p',histo_outputDirName])
        histo_outFullPath = os.path.join(histo_outputDirName,histo_outputFileName)
        histo_outputFile = ROOT.TFile.Open(histo_outFullPath,"RECREATE")

        #book histograms
        histo_dict = {}
        for histo_name in histos:
            title = str(histos[histo_name]["title"])
            nbins = int(histos[histo_name]["nbins"])
            xlow  = float(histos[histo_name]["xlow"])
            xhigh = float(histos[histo_name]["xhigh"])
            histo_model = (histo_name,title,nbins,xlow,xhigh)

            var_name = str(histos[histo_name]["var"])
            #print("---> histo {}, var {}".format(histo_name,histos[histo_name]["var"]))

            histo_dict[histo_name]= df.Histo1D(histo_model,var_name,"tot_weight")

        #Add ctau reweighted histograms
        if args.ctauReweighting and dataset_category=="signal":
            for w in [str(x) for x in df.GetColumnNames() if not str(x).find("ctau_weight_")<0]:
                weight_label = w.split("_")[-1]
                df = df.Define("tot_weight_"+weight_label,"tot_weight*"+w)
                for histo_name in histos:
                    weighted_histo_name = histo_name+"_"+weight_label
                    title = str(histos[histo_name]["title"])+"_"+weight_label
                    nbins = int(histos[histo_name]["nbins"])
                    xlow  = float(histos[histo_name]["xlow"])
                    xhigh = float(histos[histo_name]["xhigh"])
                    histo_model = (weighted_histo_name,title,nbins,xlow,xhigh)

                    var_name = str(histos[histo_name]["var"])
                    #print("---> {}".format(var_name))

                    histo_dict[weighted_histo_name]= df.Histo1D(histo_model,var_name,"tot_weight_"+weight_label)

        # write histograms on file
        for histo_name in histo_dict:
            histo_outputFile.cd()
            histo_dict[histo_name].Write()
        
        histo_outputFile.Close()
        print("Output histograms saved in {}".format(os.path.join(histo_outputDirName,histo_outputFileName)))
    
    # Fill in output reports per category
    reports[cat["label"]] = df.Report()
    if args.ctauReweighting and dataset_category=="signal":
        old_ctau_label = dataset_name_label[dataset_name_label.find("ctau")+4:dataset_name_label.find("mm")]
        hnl_mass_label = dataset_name_label[dataset_name_label.find("mN")+2:dataset_name_label.find("mN")+5]
        events = {}
        for new_ctau in selection["mN"+hnl_mass_label+"_ctau"+old_ctau_label+"mm_rw_points"]:
            weight_var = "ctau_weight_"+old_ctau_label+"TO"+str(new_ctau).replace(".","p")
            weighted_events = df.Sum(weight_var).GetValue()
            events[weight_var+"_weighted_events"] = weighted_events
            weighted_events_reports[cat["label"]] = events

print("+++ FINAL REPORT +++")
for c in reports:
    print("--> {} category".format(c))
    reports[c].Print()
    if args.ctauReweighting and dataset_category=="signal":
        for w in weighted_events_reports[c]:
            print("{} weighted events: {}".format(w,weighted_events_reports[c][w]))

print("--- Analysis completed in {} seconds ---".format(time.time() - start_time))
