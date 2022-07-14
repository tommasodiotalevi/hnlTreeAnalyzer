import sys
import ROOT


dataFileName = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/output/ds_prompt/crab/inclusive/histograms_PhiToMuMu_prompt_ParkingBPH_Run2018A_full_inclusive.root"
mcPromptFileName = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/output/ds_prompt/crab/inclusive/histograms_PhiToMuMu_prompt_DsToPhiPi_PhiToMuMu_BParkMuonFilter_TuneCP5_13TeV-pythia8-evtgen_inclusive.root"
mcNonPromptFileName = "/afs/cern.ch/work/l/llunerti/private/hnlTreeAnalyzer/output/ds_prompt/crab/inclusive/histograms_PhiToMuMu_prompt_BsToDsNuMu_DsToPhiPi_PhiToMuMu_SoftQCD_TuneCP5_13TeV-pythia8-evtgen_inclusive.root"
histoName = "h_ds_l_prop"

dataFile = ROOT.TFile.Open(dataFileName)
mcPromptFile = ROOT.TFile.Open(mcPromptFileName)
mcNonPromptFile = ROOT.TFile.Open(mcNonPromptFileName)

data = dataFile.Get(histoName)          # data histogram
mc0  = mcPromptFile.Get(histoName)      # first MC histogram
mc1  = mcNonPromptFile.Get(histoName)   # second MC histogram
mc = ROOT.TObjArray(2)                  # MC histograms are put in this array
mc.Add(mc0)
mc.Add(mc1)
fit = ROOT.TFractionFitter(data, mc) # initialise
fit.Constrain(1,0.0,1.0)             # constrain fraction 1 to be between 0 and 1
fit.SetRangeX(1,15)                  # use only the first 15 bins in the fit
status = int(fit.Fit())              # perform the fit

print("fit status: {}".format(status))
if status == 0:                       # check on fit status 
    result = fit.GetPlot()
    #data.Draw("Ep")
    #result.Draw("same")
    result.Draw()

