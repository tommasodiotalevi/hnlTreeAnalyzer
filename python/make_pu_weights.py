import ROOT
import sys

inputFileNameData = sys.argv[1]
inputFileNameMC = sys.argv[2]

# open files 
inputFileData = ROOT.TFile.Open(inputFileNameData)
inputFileMC = ROOT.TFile.Open(inputFileNameMC)

# get histograms
inputHistoData = inputFileData.Get("pileup")
inputHistoMC = inputFileMC.Get("h_pu_trueInt")

# normalize to unity
inputHistoData.Scale(1./inputHistoData.Integral())
inputHistoMC.Scale(1./inputHistoMC.Integral())

weightHisto = inputHistoData.Clone()
weightHisto.SetName("pu_weights")
weightHisto.SetTitle("pu_weights")

# make data/MC ratio
weightHisto.Divide(inputHistoMC)

# save
weightHisto.SaveAs("pu_weights.root")
