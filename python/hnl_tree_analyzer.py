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
    chain.Add(inputFileName)

cat_dict = {"chargedHnl":"C_mu1_charge_ptBest+C_pi_charge_ptBest != 0",
            "neutralHnl":"C_mu1_charge_ptBest+C_pi_charge_ptBest == 0"}

n_cuts = 20
sel_cuts = {"C_Hnl_lxy_ptBest":numpy.linspace(-1,30,n_cuts),
            "C_Hnl_pt_ptBest": numpy.linspace(-1,30,n_cuts)}

for cat in cat_dict:
    input_file_name = inputFileName_list[0].split("/")[-1].split(".")[0]
    dataset_name_label = input_file_name[input_file_name.find("_")+1:input_file_name.find("_tree")]
    
    outputFileName = "out_hnl_tree_analyzer_"+dataset_name_label+"_"+cat+".root"
    outputDirName = os.path.join(config["output_dir_name"],cat)
    
    subprocess.call(['mkdir','-p',outputDirName])
    outputFile = ROOT.TFile(os.path.join(outputDirName,outputFileName),"RECREATE")
    
    m_hnl_pt           = ("h_hnl_pt"          ,";HNL p_{T} [GeV];Events", 34, 6., 40.)
    m_hnl_preFit_mass  = ("h_hnl_preFit_mass" ,";HNL (pre-fit) mass [GeV];Events", 40, 0.25, 6.)
    m_hnl_postFit_mass = ("h_hnl_postFit_mass",";HNL (post-fit) mass [GeV];Events", 40, 0.25, 6.)
    m_hnl_lxy          = ("h_hnl_lxy"         ,";HNL L_{xy} [cm];Events", 20, 0., 40.)
    m_hnl_vtx_prob     = ("h_hnl_vtx_prob"    ,";SV probability;Events", 50, 0., 1.)
    m_hnl_vtx_cos2D    = ("h_hnl_vtx_cos2D"   ,";SV cos2D;Events", 100, -1., 1.)
    m_hnl_vtx_dispSig  = ("h_hnl_vtx_dispSig" ,";SV displacement significance;Events", 100, 0., 100.)
    
    m_mu1_pt           = ("h_mu1_pt"    ,";#mu_{1} p_{T} [GeV];Events", 24, 6., 30.)
    m_mu1_eta          = ("h_mu1_eta"   ,";#mu_{1} #eta;Events", 20, -2.5, 2.5)
    m_mu1_ip_z         = ("h_mu1_ip_z"  ,";#mu_{1} z IP [cm];Events", 20, 0., 5.)
    m_mu1_ip_xy        = ("h_mu1_ip_xy" ,";#mu_{1} xy IP [cm];Events", 20, 0., 5.)
    m_mu1_ips_z        = ("h_mu1_ips_z" ,";#mu_{1} z IPS;Events", 20, 0., 200.)
    m_mu1_ips_xy       = ("h_mu1_ips_xy",";#mu_{1} xy IPS;Events", 20, 0., 200.)
    
    m_mu2_pt           = ("h_mu2_pt"  ,";#mu_{2} p_{T} [GeV];Events", 24, 6., 30.)
    m_mu2_eta          = ("h_mu2_eta" ,";#mu_{2} #eta;Events", 20, -2.5, 2.5)
    
    m_pi_pt            = ("h_pi_pt"    ,";#pi p_{T} [GeV];Events", 40, 0., 20.)
    m_pi_eta           = ("h_pi_eta"   ,";#pi #eta;Events", 20, -2.5, 2.5)
    m_pi_ip_z          = ("h_pi_ip_z"  ,";#pi z IP [cm];Events", 20, 0., 5.)
    m_pi_ip_xy         = ("h_pi_ip_xy" ,";#pi xy IP [cm];Events", 20, 0., 5.)
    m_pi_ips_z         = ("h_pi_ips_z" ,";#pi z IPS;Events", 20, 0., 200.)
    m_pi_ips_xy        = ("h_pi_ips_xy",";#pi xy IPS;Events", 20, 0., 200.)
    
    m_pv               = ("h_pv",";Number of primary vertices;Events", 60, 0., 60.)
    m_pu_trueInt       = ("h_pu_trueInt",";Number of true PU interactions;Events", 60, 0., 60.)
    
    ROOT.EnableImplicitMT(8)
    
    df = ROOT.RDataFrame(chain)
    
    df = df.Define("weight",str(event_weight))    
    df = df.Filter(cat_dict[cat] ,cat+" category")
    df = df.Filter("mu9_ip6_matched_ptBest>0 || mu12_ip6_matched_ptBest>0" ,"mu9_ip6 or mu12_ip6 matched")
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
    
    h_hnl_pt           = df.Histo1D(m_hnl_pt          ,"C_Hnl_pt_ptBest","weight")
    h_hnl_preFit_mass  = df.Histo1D(m_hnl_preFit_mass ,"C_Hnl_preFit_mass_ptBest","weight")
    h_hnl_postFit_mass = df.Histo1D(m_hnl_postFit_mass,"C_Hnl_mass_ptBest","weight")
    h_hnl_lxy          = df.Histo1D(m_hnl_lxy         ,"C_Hnl_lxy_ptBest","weight")
    h_hnl_vtx_prob     = df.Histo1D(m_hnl_vtx_prob    ,"C_Hnl_vertex_prob_ptBest","weight")
    h_hnl_vtx_cos2D    = df.Histo1D(m_hnl_vtx_cos2D   ,"C_Hnl_vertex_cos2D_ptBest","weight")
    h_hnl_vtx_dispSig  = df.Histo1D(m_hnl_vtx_dispSig ,"C_Hnl_vertex_sig_ptBest","weight")
                                                      
    h_mu2_pt           = df.Histo1D(m_mu2_pt          ,"C_mu2_pt_ptBest","weight")
    h_mu2_eta          = df.Histo1D(m_mu2_eta         ,"C_mu2_eta_ptBest","weight")
                                          
    h_mu1_pt           = df.Histo1D(m_mu1_pt          ,"C_mu1_pt_ptBest","weight")
    h_mu1_eta          = df.Histo1D(m_mu1_eta         ,"C_mu1_eta_ptBest","weight")
    h_mu1_ip_z         = df.Histo1D(m_mu1_ip_z        ,"C_mu1_ip_z_ptBest","weight")
    h_mu1_ip_xy        = df.Histo1D(m_mu1_ip_xy       ,"C_mu1_ip_xy_ptBest","weight")
    h_mu1_ips_z        = df.Histo1D(m_mu1_ips_z       ,"C_mu1_ips_z_ptBest","weight")
    h_mu1_ips_xy       = df.Histo1D(m_mu1_ips_xy      ,"C_mu1_ips_xy_ptBest","weight")
                                                      
    h_pi_pt            = df.Histo1D(m_pi_pt           ,"C_pi_pt_ptBest","weight")
    h_pi_eta           = df.Histo1D(m_pi_eta          ,"C_pi_eta_ptBest","weight")
    h_pi_ip_z          = df.Histo1D(m_pi_ip_z         ,"C_pi_ip_z_ptBest","weight")
    h_pi_ip_xy         = df.Histo1D(m_pi_ip_xy        ,"C_pi_ip_xy_ptBest","weight")
    h_pi_ips_z         = df.Histo1D(m_pi_ips_z        ,"C_pi_ips_z_ptBest","weight")
    h_pi_ips_xy        = df.Histo1D(m_pi_ips_xy       ,"C_pi_ips_xy_ptBest","weight")
    
    h_pv               = df.Histo1D(m_pv              ,"nPV","weight")
    h_pu_trueInt       = df.Histo1D(m_pu_trueInt      ,"nPU_trueInt","weight")
    
    h_hnl_pt          .Write()
    h_hnl_preFit_mass .Write()
    h_hnl_postFit_mass.Write()
    h_hnl_lxy         .Write()
    h_hnl_vtx_prob    .Write()
    h_hnl_vtx_cos2D   .Write()
    h_hnl_vtx_dispSig .Write()
    h_mu2_pt          .Write()
    h_mu2_eta         .Write()
    h_mu1_pt          .Write()
    h_mu1_eta         .Write()
    h_mu1_ip_z        .Write()
    h_mu1_ip_xy       .Write()
    h_mu1_ips_z       .Write()
    h_mu1_ips_xy      .Write()
    h_pi_pt           .Write()
    h_pi_eta          .Write()
    h_pi_ip_z         .Write()
    h_pi_ip_xy        .Write()
    h_pi_ips_z        .Write()
    h_pi_ips_xy       .Write()
    h_pv              .Write()
    h_pu_trueInt      .Write()
    
    outputFile.Close()
    
    r = df.Report()
    print("-- REPORT --")
    r.Print()

#print("--- {} seconds ---".format(time.time() - start_time))
