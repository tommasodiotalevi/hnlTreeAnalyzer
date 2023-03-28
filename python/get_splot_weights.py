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
parser.add_argument("--johnSig",action='store_true',default=False,help="Fit the signal with a johnson")
parser.add_argument("--intLumi",default="41.6",help="Print the integrated lumi")
parser.add_argument("--addTag",type=str,default="",help="Tag the output file")
parser.add_argument("--showFitResults",action='store_true',default=False,help="Show fit results on the plot")
parser.add_argument("--fixMass",action='store_true',default=False,help="Fix D and Ds masses to PDG value")
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

padUpper = ROOT.TPad('padUpper', 'padUpper', 0, 0.3, 1, 1.0)
padUpper.SetBottomMargin(0.01)
padUpper.SetTopMargin(0.12)
padUpper.Draw()
    
padLower = ROOT.TPad('padLower', 'padLower', 0, 0.0, 1, 0.3)
padLower.SetBottomMargin(0.35)
padLower.SetTopMargin(0.05)
padLower.SetGridy()
padLower.Draw()

C_Ds_mass = ROOT.RooRealVar("C_Ds_mass","#mu#mu#pi invariant mass",1.75,2.15)
event = ROOT.RooRealVar("event","event",-1)

inputDataSet = ROOT.RooDataSet("inputDataSet","inputDataSet",ROOT.RooArgSet(C_Ds_mass,event),ROOT.RooFit.Import(inputTree))

output_full_path = os.path.join(output_dir,"splot_{}{}_output.root".format(era,n_era))
if tag != "":
    output_full_path = os.path.join(output_dir,"splot_{}{}_{}_output.root".format(era,n_era,tag))
outputFile = ROOT.TFile.Open(output_full_path,"RECREATE")

var_to_display = list()
sigpdf = str()

if not args.gausSig and not args.voigSig and not args.johnSig:

    print("Using Double crystal ball function to fit the signal peaks")

    sigpdf="dcb"

    m0_D      = ROOT.RooRealVar("m0_D"   ,"m0_D",1.8,1.9)    
    if args.fixMass:
        m0_D.setVal(1.8695)
        m0_D.setConstant(ROOT.kTRUE)

    sigma_D   = ROOT.RooRealVar("sigma_D","sigma_D",0.00005,0.1)
    alpha_D   = ROOT.RooRealVar("alpha_D","alpha_D",0.,10.)
    N_D       = ROOT.RooRealVar("N_D"    ,"N_D",0.,50.)
    
    m0_Ds     = ROOT.RooRealVar("m0_Ds"   ,"m0_Ds",1.9,2.0)
    if args.fixMass:
        m0_Ds.setVal(1.96835)
        m0_Ds.setConstant(ROOT.kTRUE)
    sigma_Ds  = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.00005,0.1)
    alpha_Ds  = ROOT.RooRealVar("alpha_Ds","alpha_Ds",0.,10.)
    N_Ds      = ROOT.RooRealVar("N_Ds"    ,"N_Ds",0.,50.)
    
    Dpeak    = ROOT.RooCBShape("Dpeak" ,"Dpeak" ,C_Ds_mass, m0_D , sigma_D ,alpha_D ,N_D)
    Dspeak   = ROOT.RooCBShape("Dspeak","Dspeak",C_Ds_mass, m0_Ds, sigma_Ds,alpha_Ds,N_Ds)

    if not args.fixMass:
        var_to_display.append(m0_D)
        var_to_display.append(m0_Ds)

    var_to_display.append(sigma_D)
    var_to_display.append(alpha_D)
    var_to_display.append(N_D)
    var_to_display.append(sigma_Ds)
    var_to_display.append(alpha_Ds)
    var_to_display.append(N_Ds)

elif args.gausSig:

    print("Using Gaussian function to fit the signal peaks")

    sigpdf="gau"

    mean_D    = ROOT.RooRealVar("mean_D","mean_D",1.8,1.9)
    if args.fixMass:
        mean_D.setVal(1.8695)
        mean_D.setConstant(ROOT.kTRUE)
    sigma_D   = ROOT.RooRealVar("sigma_D","sigma_D",0.001,0.1)
    mean_Ds   = ROOT.RooRealVar("mean_Ds","mean_Ds",1.9,2.0) 
    if args.fixMass:
        mean_Ds.setVal(1.96835)
        mean_Ds.setConstant(ROOT.kTRUE)
    sigma_Ds  = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.001,0.1)

    Dpeak = ROOT.RooGaussian("Dpeak","Dpeak",C_Ds_mass, mean_D, sigma_D)
    Dspeak = ROOT.RooGaussian("Dspeak","Dspeak",C_Ds_mass, mean_Ds, sigma_Ds)
    
    if not args.fixMass:
        var_to_display.append(mean_D)
        var_to_display.append(mean_Ds)
    var_to_display.append(sigma_D)
    var_to_display.append(sigma_Ds)

elif args.voigSig:

    print("Using Voigtian function to fit the signal peaks")

    sigpdf="voi"

    mean_D   = ROOT.RooRealVar("mean_D","mean_D",1.8,1.9)   
    if args.fixMass:
        mean_D.setVal(1.8695)
        mean_D.setConstant(ROOT.kTRUE)
    sigma_D  = ROOT.RooRealVar("sigma_D","sigma_D",0.001,0.01)
    width_D  = ROOT.RooRealVar("width_D","width_D",0.001,0.01)

    mean_Ds  = ROOT.RooRealVar("mean_Ds","mean_Ds",1.9,2.0)   
    if args.fixMass:
        mean_Ds.setVal(1.96835)
        mean_Ds.setConstant(ROOT.kTRUE)
    sigma_Ds = ROOT.RooRealVar("sigma_Ds","sigma_Ds",0.001,0.1)
    width_Ds = ROOT.RooRealVar("width_Ds","width_Ds",0.001,0.1)

    if not args.fixMass:
        var_to_display.append(mean_D)
        var_to_display.append(mean_Ds)
    var_to_display.append(sigma_D)
    var_to_display.append(width_D)
    var_to_display.append(sigma_Ds)
    var_to_display.append(width_Ds)

    Dpeak = ROOT.RooVoigtian("Dpeak","Dpeak",C_Ds_mass, mean_D, width_D, sigma_D)
    Dspeak = ROOT.RooVoigtian("Dspeak","Dspeak",C_Ds_mass, mean_Ds, width_Ds, sigma_Ds)

elif args.johnSig:

    print("Using Johnson function to fit the signal peaks")
    sigpdf="joh"

    mu_D   = ROOT.RooRealVar("mu_D","mu_D",1.82,1.88)   
    if args.fixMass:
        mu_D.setVal(1.8695)
        mu_D.setConstant(ROOT.kTRUE)    
    lambda_D  = ROOT.RooRealVar("lambda_D","lambda_D",0.01,3.)
    gamma_D  = ROOT.RooRealVar("gamma_D","gamma_D",-10.,10.)
    delta_D  = ROOT.RooRealVar("delta_D","delta_D",1.,3.)

    mu_Ds  = ROOT.RooRealVar("mu_Ds","mu_Ds",1.9,2.0)   
    if args.fixMass:
        mu_Ds.setVal(1.96835)
        mu_Ds.setConstant(ROOT.kTRUE)
    lambda_Ds  = ROOT.RooRealVar("lambda_Ds","lambda_Ds",0.01,3.)
    gamma_Ds  = ROOT.RooRealVar("gamma_Ds","gamma_Ds",-10.,10.)
    delta_Ds  = ROOT.RooRealVar("delta_Ds","delta_Ds",1.,3.)

    if not args.fixMass:
        var_to_display.append(mu_D)
        var_to_display.append(mu_Ds)
    var_to_display.append(lambda_D)
    var_to_display.append(gamma_D)
    var_to_display.append(delta_D)
    var_to_display.append(lambda_Ds)
    var_to_display.append(gamma_Ds)
    var_to_display.append(delta_Ds)

    Dpeak = ROOT.RooJohnson("Dpeak","Dpeak",C_Ds_mass,mu_D,lambda_D,gamma_D,delta_D)
    Dspeak = ROOT.RooJohnson("Dspeak","Dspeak",C_Ds_mass, mu_Ds,lambda_Ds,gamma_Ds,delta_Ds)


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
    t1 = ROOT.TPaveText(0.65,0.55,0.88,0.86,"brNDC")
    t1.AddText("#chi^{{2}}/dof = {}".format(round(chiSquare,2)))
    t1.AddText("N_{{D_{{s}}}} = {} #pm {}".format(round(n_Ds.getValV(),0),round(n_Ds.getError(),0)))
    t1.AddText("N_{{D}} = {} #pm {}".format(round(n_D.getValV(),0),round(n_D.getError(),0)))
    for x in var_to_display:
        t1.AddText("{} = {} #pm {}".format(x.GetName(),round(x.getValV(),4),round(x.getError(),4)))
    t1.SetTextColor(ROOT.kBlack)
    
    #Draw fitted functions
    print("chi-square/ndof = {} \n".format(frame.chiSquare()))
    model.plotOn(frame,ROOT.RooFit.Name("Dpeak")      ,ROOT.RooFit.Components("Dpeak") ,ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.DrawOption("l"))
    model.plotOn(frame,ROOT.RooFit.Name("Dspeak")     ,ROOT.RooFit.Components("Dspeak"),ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.DrawOption("l"))
    model.plotOn(frame,ROOT.RooFit.Name("background"),ROOT.RooFit.Components("bkg") ,ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kBlue))
    model.plotOn(frame,ROOT.RooFit.Name("model"),ROOT.RooFit.LineColor(ROOT.kBlue))

    #Build legend
    leg = ROOT.TLegend(0.12,0.68,0.3,0.86)
    leg.AddEntry(frame.findObject("Dspeak"),"D and D_{s}","l")
    leg.AddEntry(frame.findObject("background"),"Background","l")
    leg.AddEntry(frame.findObject("model"),"Full fit","L")
    leg.AddEntry(frame.findObject("inputDataSet"),"Data","EP")
    leg.SetBorderSize(0)
    
    
    padUpper.cd()
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("")
    binwidth = frame.getFitRangeBinW()
    frame.GetYaxis().SetTitle("Events/{} GeV".format(f'{binwidth:.1g}'))
    frame.GetYaxis().SetTitleSize(0.035)
    frame.GetYaxis().SetTitleOffset(0.7)
    frame.GetYaxis().SetMaxDigits(2)
    frame.Draw()
    if args.showFitResults:
       t1.Draw("same")
    leg.Draw("same")
    latex = ROOT.TLatex()
    latex.SetTextAlign(12)
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(0.70,0.9,str(int_lumi)+" fb^{-1} (13 TeV)")
    
    #Draw pulls
    hpull = frame.pullHist()
    padLower.cd()
    pullframe = C_Ds_mass.frame()
    pullframe.SetTitle("")
    pullframe.GetXaxis().SetTitle("m(#mu#mu#pi) [GeV]")
    pullframe.GetYaxis().SetTitle("Pull")
    pullframe.GetYaxis().SetTitleOffset(0.5)
    pullframe.GetYaxis().CenterTitle(True)  
    pullframe.GetYaxis().SetTitleSize(0.1)
    pullframe.GetXaxis().SetTitleSize(0.1)
    pullframe.GetYaxis().SetLabelSize(0.08)
    pullframe.GetXaxis().SetLabelSize(0.09)
    pullframe.addPlotable(hpull, "P")
    pullframe.Draw() 
    
    outfilename = "fitted_{}{}_{}_mass".format(era,n_era,sigpdf)
    if tag != "":
        outfilename = "fitted_{}{}_{}_{}_mass".format(era,n_era,sigpdf,tag)
    if args.fixMass:
        outfilename = outfilename.replace("mass","fixM_mass")
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
    
