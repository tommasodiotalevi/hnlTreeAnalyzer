import argparse
import sys
import os
import subprocess

do_plot = False

output_dir = "splot_test"
subprocess.call("mkdir -p "+output_dir,shell=True)

parser = argparse.ArgumentParser(description="")
parser.add_argument("filename",help="Path to the input rootfile")
parser.add_argument("treeName",help="Name of the input tree")
args = parser.parse_args()

import ROOT

#Import TH1F into a RooDataHist	
inputFile  = ROOT.TFile(args.filename)
inputTree = inputFile.Get(args.treeName)

w = ROOT.RooWorkspace("w","workspace")

c = ROOT.TCanvas("c","c",800,600)

C_Ds_mass = ROOT.RooRealVar("C_Ds_mass","#mu#mu#pi invariant mass",1.8,2.1)

inputDataSet = ROOT.RooDataSet("inputDataSet","inputDataSet",ROOT.RooArgSet(C_Ds_mass),ROOT.RooFit.Import(inputTree))

#Defining variables
mean_D    = ROOT.RooRealVar("mean_D","mean_D",1.8695)    
sigma_D   = ROOT.RooRealVar("sigma_D","sigma_D",0.01,0.1)
mean_Ds   = ROOT.RooRealVar("mean_Ds","mean_Ds",1.96835)   
sigma_Ds  = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.01,0.1)
c_exp     = ROOT.RooRealVar("c_exp","c_exp",-1.5,-1.1)
coeff     = ROOT.RooRealVar("coeff","coeff",0.,1.)
n_D       = ROOT.RooRealVar("n_D","n_D",0,1000)
n_Ds      = ROOT.RooRealVar("n_Ds","n_Ds",0,1000)
n_bkg     = ROOT.RooRealVar("n_bkg","n_bkg",0,1000)

#Defining signal and background pdfs
gaussian_D = ROOT.RooGaussian("gaussian_D","gaussian_D",C_Ds_mass, mean_D, sigma_D)
gaussian_Ds = ROOT.RooGaussian("gaussian_Ds","gaussian_Ds",C_Ds_mass, mean_Ds, sigma_Ds)
expo       = ROOT.RooExponential("expo","expo",C_Ds_mass,c_exp)
model      = ROOT.RooAddPdf("model","D_{(s)}#rightarrow #phi#pi decay model", ROOT.RooArgList(gaussian_D,gaussian_Ds, expo), ROOT.RooArgList(n_D,n_Ds,n_bkg))

frame = C_Ds_mass.frame()

#Fitting the data with the model
#result = model.fitTo(inputDataSet,ROOT.RooFit.SumW2Error(ROOT.kTRUE),ROOT.RooFit.Save())
result = model.fitTo(inputDataSet,ROOT.RooFit.Save(),ROOT.RooFit.Extended())

print("before ----->")
inputDataSet.Print()

sData = ROOT.RooStats.SPlot("sData", "An SPlot", inputDataSet, model, ROOT.RooArgList(n_D,n_Ds,n_bkg))

print("after ----->")
inputDataSet.Print()


outputFile = ROOT.TFile.Open(os.path.join(output_dir,"fitted_"+args.treeName+"_parameters.root"),"RECREATE")
inputDataSet.convertToTreeStore()

if do_plot:
    #Plot and fit a RooDataHist
    inputDataSet.plotOn(frame)
    model.plotOn(frame)
    
    #Draw fitted parameters
    chiSquare = frame.chiSquare()
    t1 = ROOT.TPaveText(0.67,0.7,0.9,0.88,"brNDC")
    t1.AddText("#chi^{{2}}/dof = {}".format(round(chiSquare,2)))
    t1.AddText("N_{{D_{{s}}}} = {} #pm {}".format(round(n_Ds.getValV(),0),round(n_Ds.getError(),0)))
    t1.AddText("N_{{D}} = {} #pm {}".format(round(n_D.getValV(),0),round(n_D.getError(),0)))
    t1.AddText("N_{{bkg}} = {} #pm {}".format(round(n_bkg.getValV(),0),round(n_bkg.getError(),0)))
    t1.AddText("M_{{D_{{s}}}} = {} #pm {} MeV".format(round(mean_Ds.getValV()*1e3,1),round(mean_Ds.getError()*1e3,1)))
    t1.AddText("M_{{D}} = {} #pm {} MeV".format(round(mean_D.getValV()*1e3,1),round(mean_D.getError()*1e3,1)))
    t1.SetTextColor(ROOT.kBlack)
    
    #Draw fitted functions
    print("chi-square/ndof = {} \n".format(frame.chiSquare()))
    model.plotOn(frame,ROOT.RooFit.Name("gaussian_D"),ROOT.RooFit.Components("gaussian_D"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
    model.plotOn(frame,ROOT.RooFit.Name("gaussian_Ds"),ROOT.RooFit.Components("gaussian_Ds"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kViolet))
    model.plotOn(frame,ROOT.RooFit.Name("background"),ROOT.RooFit.Components("expo")      ,ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
    
    #Build legend
    leg = ROOT.TLegend(0.12,0.75,0.25,0.88)
    leg.AddEntry(frame.findObject("gaussian_Ds"),"D_{s}","l")
    leg.AddEntry(frame.findObject("gaussian_D"),"D","l")
    leg.AddEntry(frame.findObject("background"),"bkg","l")
    leg.SetBorderSize(0)
    
    c.cd()
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("m_{#mu#mu#pi} [GeV]")
    frame.Draw()
    t1.Draw("same")
    leg.Draw("same")
    #w.writeToFile(os.path.join(output_dir,"model_"+args.treeName+".root"))
    c.SaveAs(os.path.join(output_dir,"fitted_"+args.treeName+".png"))
    c.SaveAs(os.path.join(output_dir,"fitted_"+args.treeName+".root"))
    inputFile.Close()
    c.Write()
    

getattr(w,'import')(model)
w.Write()
sData.Write()
outputFile.Close()
    
