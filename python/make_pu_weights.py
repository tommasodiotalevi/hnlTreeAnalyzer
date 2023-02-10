import ROOT
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("puFileData", help="Path to the input file containing data pu profile", type=str)
parser.add_argument("puFileMC"  , help="Path to the input file containing MC pu profile", type=str)

args = parser.parse_args()

inputFileNameData = args.puFileData
inputFileNameMC = args.puFileMC

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
