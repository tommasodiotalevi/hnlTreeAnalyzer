import argparse
import sys
import os
import subprocess

output_dir = "splot_output"
subprocess.call("mkdir -p "+output_dir,shell=True)

parser = argparse.ArgumentParser(description="")
parser.add_argument("filename",help="Path to the input rootfile")
parser.add_argument("treeName",help="Name of the input tree")
parser.add_argument("--doPlot",action='store_true',default=False,help="Save fitted plot")
parser.add_argument("--gausSig",action='store_true',default=False,help="Fit the signal with a gaussian")
parser.add_argument("--voigSig",action='store_true',default=False,help="Fit the signal with a voigtian")
parser.add_argument("--intLumi",default="41.6",help="Print the integrated lumi")
parser.add_argument("--addTag",type=str,default="",help="Tag the output file")
parser.add_argument("--showFitResults",action='store_true',default=False,help="Show fit results on the plot")
args = parser.parse_args()

if args.gausSig and args.voigSig:
    print("You can choose either --gausSig or --voigSig fitting option, exiting...")
    sys.exit(1)

import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)

n_era = args.filename[args.filename.find("ParkingBPH")+len("ParkingBPH")]
era   = args.filename[args.filename.find("Run2018")+len("Run2018")]
int_lumi = args.intLumi
tag = args.addTag

#Import TH1F into a RooDataHist	
inputFile  = ROOT.TFile(args.filename)
inputTree = inputFile.Get(args.treeName)

w = ROOT.RooWorkspace("w","workspace")

c = ROOT.TCanvas("c","c",600,600)

C_Ds_mass = ROOT.RooRealVar("C_Ds_mass","#mu#mu#pi invariant mass",1.75,2.15)
event = ROOT.RooRealVar("event","event",-1)

inputDataSet = ROOT.RooDataSet("inputDataSet","inputDataSet",ROOT.RooArgSet(C_Ds_mass,event),ROOT.RooFit.Import(inputTree))

output_full_path = os.path.join(output_dir,"splot_{}{}_output.root".format(era,n_era))
if tag != "":
    output_full_path = os.path.join(output_dir,"splot_{}{}_{}_output.root".format(era,n_era,tag))
outputFile = ROOT.TFile.Open(output_full_path,"RECREATE")

#Defining variables

if not args.gausSig and not args.voigSig:

    m0_D      = ROOT.RooRealVar("m0_D"   ,"m0_D",1.8695)    
    sigma_D   = ROOT.RooRealVar("sigma_D","sigma_D",0.00005,0.1)
    alpha_D   = ROOT.RooRealVar("alpha_D","alpha_D",0.,10.)
    N_D       = ROOT.RooRealVar("N_D"    ,"N_D",0.,50.)
    
    m0_Ds     = ROOT.RooRealVar("m0_Ds"   ,"m0_Ds",1.96835)    
    sigma_Ds  = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.00005,0.1)
    alpha_Ds  = ROOT.RooRealVar("alpha_Ds","alpha_Ds",0.,10.)
    N_Ds      = ROOT.RooRealVar("N_Ds"    ,"N_Ds",0.,50.)
    
    Dpeak    = ROOT.RooCBShape("Dpeak" ,"Dpeak" ,C_Ds_mass, m0_D , sigma_D ,alpha_D ,N_D)
    Dspeak   = ROOT.RooCBShape("Dspeak","Dspeak",C_Ds_mass, m0_Ds, sigma_Ds,alpha_Ds,N_Ds)

elif args.gausSig:
    mean_D    = ROOT.RooRealVar("mean_D","mean_D",1.8695)    
    sigma_D   = ROOT.RooRealVar("sigma_D","sigma_D",0.001,0.1)
    mean_Ds   = ROOT.RooRealVar("mean_Ds","mean_Ds",1.96835)   
    sigma_Ds  = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.001,0.1)

    Dpeak = ROOT.RooGaussian("Dpeak","Dpeak",C_Ds_mass, mean_D, sigma_D)
    Dspeak = ROOT.RooGaussian("Dspeak","Dspeak",C_Ds_mass, mean_Ds, sigma_Ds)

elif args.voigSig:
    mean_D   = ROOT.RooRealVar("mean_D","mean_D",1.8695)    
    sigma_D  = ROOT.RooRealVar("sigma_D","sigma_D",0.001,0.01)
    width_D  = ROOT.RooRealVar("width_D","width_D",0.001,0.01)
    mean_Ds  = ROOT.RooRealVar("mean_Ds","mean_Ds",1.96835)   
    sigma_Ds = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.001,0.1)
    width_Ds = ROOT.RooRealVar("width_Ds","width_Ds",0.001,0.1)

    Dpeak = ROOT.RooVoigtian("Dpeak","Dpeak",C_Ds_mass, mean_D, width_D, sigma_D)
    Dspeak = ROOT.RooVoigtian("Dspeak","Dspeak",C_Ds_mass, mean_Ds, width_Ds, sigma_Ds)


c_exp     = ROOT.RooRealVar("c_exp","c_exp",-5,1.5)
bkg      = ROOT.RooExponential("bkg","bkg",C_Ds_mass,c_exp)

n_D       = ROOT.RooRealVar("n_D","n_D",0,100000)
n_Ds      = ROOT.RooRealVar("n_Ds","n_Ds",0,100000)
n_bkg     = ROOT.RooRealVar("n_bkg","n_bkg",0,100000000)


#Defining signal and background pdfs
model      = ROOT.RooAddPdf("model","model", ROOT.RooArgList(Dpeak,Dspeak, bkg), ROOT.RooArgList(n_D,n_Ds,n_bkg))

#Fitting the data with the model
result = model.fitTo(inputDataSet,ROOT.RooFit.Save(),ROOT.RooFit.Extended(),ROOT.RooFit.SumW2Error(ROOT.kTRUE))

sData = ROOT.RooStats.SPlot("sData", "An SPlot", inputDataSet, model, ROOT.RooArgList(n_D,n_Ds,n_bkg))

inputDataSet.convertToTreeStore()


if args.doPlot:
    
    #Plot and fit a RooDataHist
    frame = C_Ds_mass.frame()
    inputDataSet.plotOn(frame)
    model.plotOn(frame,ROOT.RooFit.Name("inputDataSet"),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.LineColor(ROOT.kBlack),ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineWidth(1))
    
    #Draw fitted parameters
    chiSquare = frame.chiSquare()
    t1 = ROOT.TPaveText(0.65,0.7,0.88,0.88,"brNDC")
    t1.AddText("#chi^{{2}}/dof = {}".format(round(chiSquare,2)))
    t1.AddText("N_{{D_{{s}}}} = {} #pm {}".format(round(n_Ds.getValV(),0),round(n_Ds.getError(),0)))
    t1.AddText("N_{{D}} = {} #pm {}".format(round(n_D.getValV(),0),round(n_D.getError(),0)))
    t1.AddText("N_{{bkg}} = {} #pm {}".format(round(n_bkg.getValV(),0),round(n_bkg.getError(),0)))
    #t1.AddText("M_{{D_{{s}}}} = {} #pm {} MeV".format(round(m0_Ds.getValV()*1e3,1),round(m0_Ds.getError()*1e3,1)))
    #t1.AddText("M_{{D}} = {} #pm {} MeV".format(round(m0_D.getValV()*1e3,1),round(m0_D.getError()*1e3,1)))
    t1.SetTextColor(ROOT.kBlack)
    
    #Draw fitted functions
    print("chi-square/ndof = {} \n".format(frame.chiSquare()))
    #model.plotOn(frame,ROOT.RooFit.Name("gaussian_D"),ROOT.RooFit.Components("gaussian_D"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
    #model.plotOn(frame,ROOT.RooFit.Name("gaussian_Ds"),ROOT.RooFit.Components("gaussian_Ds"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kViolet))
    #model.plotOn(frame,ROOT.RooFit.Name("bkg"),ROOT.RooFit.Components("bkg"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))
    #model.plotOn(frame,ROOT.RooFit.Name("Dpeak")      ,ROOT.RooFit.Components("Dpeak") ,ROOT.RooFit.LineWidth(0),ROOT.RooFit.VLines(),ROOT.RooFit.FillStyle(3481),ROOT.RooFit.FillColor(ROOT.kViolet),ROOT.RooFit.DrawOption("F"),ROOT.RooFit.MoveToBack())
    model.plotOn(frame,ROOT.RooFit.Name("Dpeak")      ,ROOT.RooFit.Components("Dpeak") ,ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.DrawOption("l"))
    model.plotOn(frame,ROOT.RooFit.Name("Dspeak")     ,ROOT.RooFit.Components("Dspeak"),ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.DrawOption("l"))
    model.plotOn(frame,ROOT.RooFit.Name("background"),ROOT.RooFit.Components("bkg") ,ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kBlue))
    model.plotOn(frame,ROOT.RooFit.Name("model"),ROOT.RooFit.LineColor(ROOT.kBlue))

    #Build legend
    leg = ROOT.TLegend(0.12,0.7,0.3,0.88)
    #leg.AddEntry(frame.findObject("gaussian_Ds"),"D_{s}","l")
    #leg.AddEntry(frame.findObject("gaussian_D"),"D","l")
    leg.AddEntry(frame.findObject("Dspeak"),"D and D_{s}","l")
    #leg.AddEntry(frame.findObject("Dpeak"),"D","l")
    leg.AddEntry(frame.findObject("background"),"Background","l")
    leg.AddEntry(frame.findObject("model"),"Full fit","L")
    leg.AddEntry(frame.findObject("inputDataSet"),"Data","EP")
    leg.SetBorderSize(0)
    
    
    c.cd()
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("m(#mu#mu#pi) [GeV]")
    binwidth = frame.getFitRangeBinW()
    frame.GetYaxis().SetTitle("Events/{} GeV".format(f'{binwidth:.1g}'))
    frame.GetYaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetMaxDigits(2)
    frame.Draw()
    if args.showFitResults:
       t1.Draw("same")
    leg.Draw("same")
    latex = ROOT.TLatex()
    latex.SetTextAlign(12)
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(0.67,0.92,str(int_lumi)+" fb^{-1} (13 TeV)")
    #latex.DrawLatex(2.042,212.5,str(int_lumi)+" fb^{-1} (13 TeV)")
    outfilename = "fitted_{}{}_mass".format(era,n_era)
    if tag != "":
        outfilename = "fitted_{}{}_{}_mass".format(era,n_era,tag)
    c.SaveAs(os.path.join(output_dir,outfilename+".png"))
    c.SaveAs(os.path.join(output_dir,outfilename+".pdf"))
    c.SaveAs(os.path.join(output_dir,outfilename+".root"))
    inputFile.Close()
    c.Write()
    

getattr(w,'import')(model)
getattr(w,'import')(inputDataSet)
w.Write()
sData.Write()
outputFile.Close()

print("Output file saved in {} ".format(output_full_path))
    
