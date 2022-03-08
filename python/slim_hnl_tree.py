import ROOT
import sys
import os
import json
import subprocess
import time
start_time = time.time()

configFileName = sys.argv[1]
dataset_to_slim = sys.argv[2]

with open(configFileName, "r") as f:
    config = json.loads(f.read())

ROOT.EnableImplicitMT(8)

#get input files
inputFileName_list = config[dataset_to_slim]["file_name_list"]
dataset_category = config[dataset_to_slim]["dataset_category"]

chain = ROOT.TChain('wztree')
for inputFileName in inputFileName_list:
    chain.Add(inputFileName)

input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.find("_tree")]


df = ROOT.RDataFrame(chain)
df = df.Define("C_Hnl_pt" ,"sqrt(C_Hnl_px*C_Hnl_px + C_Hnl_py*C_Hnl_py)")\
       .Define("C_Hnl_lxy","sqrt(C_Hnl_vertex_x*C_Hnl_vertex_x + C_Hnl_vertex_y*C_Hnl_vertex_y)")\
       .Define("C_mu1_pt" ,"sqrt(C_mu1_px*C_mu1_px + C_mu1_py*C_mu1_py)")\
       .Define("C_mu2_pt" ,"sqrt(C_mu2_px*C_mu2_px + C_mu2_py*C_mu2_py)")\
       .Define("C_pi_pt"  ,"sqrt(C_pi_px*C_pi_px + C_pi_py*C_pi_py)")\
       .Define("mask","ArgMax(C_Hnl_pt)")
       #.Define("mask","ArgMax(C_Hnl_vertex_cos2D)")

slimmed_col_names = []

#keep best hnl pt combo from each event
for c in df.GetColumnNames():
    col_name = str(c)
    col_type = df.GetColumnType(col_name)
    if col_type.find("ROOT::VecOps")<0:
        slimmed_col_names.append(col_name)
        continue
    #df = df.Define(col_name+"_cos2DBest",col_name+"[mask]")
    #slimmed_col_names.append(str(col_name+"_cos2DBest"))
    df = df.Define(col_name+"_ptBest",col_name+"[mask]")
    slimmed_col_names.append(str(col_name+"_ptBest"))

outputFileName = "slimmed_"+dataset_name_label+"_tree.root"
outputDirName = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/slimmed_tree/test"
subprocess.call(['mkdir','-p',outputDirName])
outputPath = os.path.join(outputDirName,outputFileName)

print("Saving {}".format(outputPath))

df.Snapshot("slimmed_tree",outputPath,slimmed_col_names)

#update cfg file for hnl tree analyzer with location of slimmed hnl tree
out_cfg = {}
with open(configFileName, "r") as f:
    out_cfg = json.loads(f.read())

dataset_dic = out_cfg[dataset_to_slim]
dataset_dic["slimmed_file_name_list"] = [str(outputPath)]

out_cfg[dataset_to_slim] = dataset_dic

with open(configFileName, "w") as f:
    json.dump(out_cfg,f, indent=4, sort_keys=True)

print("--- %s seconds ---" % (time.time() - start_time))
