import argparse
import sys
parser = argparse.ArgumentParser(description="")
parser.add_argument("filename",help="Path to the input rootfile")
parser.add_argument("histoName",help="Name of the input histogram")
args = parser.parse_args()

import ROOT

#Import TH1F into a RooDataHist	
inputFile  = ROOT.TFile(args.filename)
inputHisto = ROOT.TH1F(inputFile.Get(args.histoName))
histoTitle = inputHisto.GetTitle()
w = ROOT.RooWorkspace("w","workspace")

c = ROOT.TCanvas("c","c",800,600)

invMass = ROOT.RooRealVar("invMass","#mu#mu#pi invariant mass",1.8,2.1)
inputDataHist = ROOT.RooDataHist("inputDataHist","inputDataHist",ROOT.RooArgList(invMass),inputHisto)

#Defining variables
mean_1    = ROOT.RooRealVar("mean_1","mean_1",1.96,1.98) 
sigma_1   = ROOT.RooRealVar("sigma_1","sigma_1",0.001,0.05)
mean_2    = ROOT.RooRealVar("mean_2","mean_2",1.96,1.98)
sigma_2   = ROOT.RooRealVar("sigma_2","sigma_2",0.001,0.05)
c_exp     = ROOT.RooRealVar("c_exp","c_exp",-1.5,-1.1)
coeff_1   = ROOT.RooRealVar("coeff_1","coeff_1",0.,1.)
n_sig     = ROOT.RooRealVar("n_sig","n_sig",10,100000)
n_bkg     = ROOT.RooRealVar("n_bkg","n_bkg",10,100000)

#Defining signal and background pdfs
gaussian_1 = ROOT.RooGaussian("gaussian_1","gaussian_1",invMass, mean_1, sigma_1)
gaussian_2 = ROOT.RooGaussian("gaussian_2","gaussian_2",invMass, mean_2, sigma_2)
expo       = ROOT.RooExponential("expo","expo",invMass,c_exp)
sig        = ROOT.RooAddPdf("sig","sig",ROOT.RooArgList(gaussian_1,gaussian_2),ROOT.RooArgList(coeff_1),ROOT.kTRUE)
model      = ROOT.RooAddPdf("model","D_{s}#rightarrow #phi#pi decay model", ROOT.RooArgList(sig, expo), ROOT.RooArgList(n_sig,n_bkg))

frame = invMass.frame()

#Fitting the data with the model
result = model.fitTo(inputDataHist,ROOT.RooFit.SumW2Error(ROOT.kTRUE),ROOT.RooFit.Save())

#Plot and fit a RooDataHist
inputDataHist.plotOn(frame)
model.plotOn(frame)
getattr(w,'import')(model)

#Draw fitted parameters
chiSquare = frame.chiSquare()
t1 = ROOT.TPaveText(0.67,0.7,0.9,0.88,"brNDC")
t1.AddText("#chi^{{2}}/dof = {}".format(round(chiSquare,2)))
t1.AddText("N_{{sig}} = {} #pm {}".format(round(n_sig.getValV(),0),round(n_sig.getError(),0)))
#t1.AddText("M_{{D_{{s}}}} = {} #pm {} MeV".format(round(mean_2.getValV()*1e3,1),round(mean_2.getError()*1e3,1)))
#t1.AddText("M_{{D}} = {} #pm {} MeV".format(round(mean_1.getValV()*1e3,1),round(mean_1.getError()*1e3,1)))
t1.SetTextColor(ROOT.kBlack)

#Draw fitted functions
print("chi-square/ndof = {} \n".format(frame.chiSquare()))
model.plotOn(frame,ROOT.RooFit.Name("gaussian_1"),ROOT.RooFit.Components("gaussian_1"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame,ROOT.RooFit.Name("gaussian_2"),ROOT.RooFit.Components("gaussian_2"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed))
model.plotOn(frame,ROOT.RooFit.Name("background"),ROOT.RooFit.Components("expo")      ,ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen))

#Build legend
leg = ROOT.TLegend(0.12,0.75,0.25,0.88)
leg.AddEntry(frame.findObject("gaussian_2"),"D_{s}","l")
leg.AddEntry(frame.findObject("background"),"bkg","l")
leg.SetBorderSize(0)

c.cd()
frame.SetTitle("")
frame.GetXaxis().SetTitle("m_{#mu#mu#pi} [GeV]")
frame.Draw()
t1.Draw("same")
leg.Draw("same")
w.writeToFile("Out_fitter/model"+args.histoName+".root")
c.SaveAs("Out_fitter/fitted_"+args.histoName+"_mass.png")
c.SaveAs("Out_fitter/fitted_"+args.histoName+"_mass.root")
inputFile.Close()
outputFile = ROOT.TFile.Open("Out_fitter/fitted_"+args.histoName+"_parameters.root","RECREATE")
c.Write()
result.Write()
outputFile.Close()

