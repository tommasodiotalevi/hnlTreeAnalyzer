import math
import sys
import ROOT
import pandas as pd

#dummy class that contains all the physics constants that I need
class Constants:

    Ds_mass = 1.96835 # [GeV] from pdg
    mu_mass = 0.1056583745 # [GeV] from pdg
    G_f = 1.1663787e-5 # [1/GeV^2] from pdg https://pdg.lbl.gov/2019/reviews/rpp2019-rev-phys-constants.pdf
    f_Ds = 0.2498 # [GeV] table IV https://arxiv.org/pdf/1509.02220.pdf
    f_pi = 0.1302 # [GeV] table I https://arxiv.org/pdf/1509.02220.pdf
    V_cs = 0.987 # from pdg
    V_ud = 0.97370 # from pdg
    Ds_meanLife = 504e-15 # [s] from pdg
    hbar = 6.582119569e-25 # [Gev s] 
    light_speed = 299792458e+3 # [mm/s]
    BR_DsToPhiPi = 0.045 # from pdg
    BR_PhiToMuMu = 0.000286 # from pdg

    def __init__(self):
        pass

# from eq A.11 https://arxiv.org/pdf/1805.08567.pdf
def compute_BR_DsToN (m_HNL,v2):
    constants = Constants()
    square_brackets = (float(m_HNL)/constants.Ds_mass)**2+(constants.mu_mass/constants.Ds_mass)**2-((float(m_HNL)/constants.Ds_mass)**2-(constants.mu_mass/constants.Ds_mass)**2)**2
    square_root     = math.sqrt(1 + (float(m_HNL)/constants.Ds_mass)**4 + (constants.mu_mass/constants.Ds_mass)**4 - 2*((float(m_HNL)/constants.Ds_mass)**2) - 2*((constants.mu_mass/constants.Ds_mass)**2) -2*((float(m_HNL)/constants.Ds_mass)**2)*((constants.mu_mass/constants.Ds_mass)**2))
    gamma           = (((constants.G_f)**2)*(constants.f_Ds**2)*(constants.Ds_mass**3)*(constants.V_cs**2)*v2*square_brackets*square_root)/(8*math.pi)
    gamma_Ds        = (constants.hbar)/(constants.Ds_meanLife)
    print("*** Computing BR(Ds-->NMu)")
    print("HNL mass: {}".format(m_HNL))
    print("|V|^2: {}".format(v2))
    print("BR(Ds-->NMu)  = {}".format(gamma/gamma_Ds))
    print('\n')
    return gamma/gamma_Ds

# ratio eq.2/eq.13 https://journals-aps-org.ezproxy.unibo.it/prd/pdf/10.1103/PhysRevD.94.113007
def compute_BR_NToPiMu (m_HNL):
    constants = Constants()
    br = (96.0*(math.pi**2)*(constants.f_pi**2)*(constants.V_ud**2))/(16.0*10.95*(float(m_HNL)**2))
    print("*** Computing BR(N-->MuPi)")
    print("HNL mass: {}".format(m_HNL))
    print("BR(N-->MuPi) = {}".format(br))
    print('\n')
    return br

def compute_v2_from_ctau(m_HNL,ctau):
    constants = Constants()
    num = constants.light_speed*constants.hbar
    den = ((constants.G_f**2)*(m_HNL**5)*ctau)/(96*(math.pi**3))
    return num/den

def get_expected_signal_yield(m_HNL,ctau,nDs,fPrompt,effHnl,effDs):
    constants = Constants()
    BR_DsToNMu = compute_BR_DsToN(m_HNL,compute_v2_from_ctau(m_HNL,ctau))
    BR_NToPiMu = compute_BR_NToPiMu(m_HNL)
    BR_factor = (BR_DsToNMu*BR_NToPiMu)/(constants.BR_DsToPhiPi*constants.BR_PhiToMuMu)
    print("*** Computing (BR(Ds->Hnl Mu)*BR(Hnl->Mu Pi))/(BR(Ds->Phi Pi)*BR(Phi->Mu Mu))")
    print("BR factor = (BR_DsToNMu*BR_NToPiMu)/(BR_DsToPhiPi*BR_PhiToMuMu) = {}".format(BR_factor))
    print('\n')
    nHnl = nDs*fPrompt*BR_factor*(effHnl/effDs)
    return nHnl

def get_yield_from_workspace(fileName,wsName,varName):
    infile = ROOT.TFile.Open(fileName)
    inws   = infile.Get(str(wsName))
    var    = inws.var(varName)
    val    = var.getVal()
    return val

def get_yield_from_tree(fileName,treeName):
    infile = ROOT.TFile.Open(fileName)
    intree = infile.Get(treeName)
    val    = intree.GetEntries()
    return val
    
def get_yield_from_csv(fileName):
    t = ROOT.TTree()
    t.ReadFile(fileName)
    val = t.GetEntries()
    return val

def get_weighted_yield_from_csv(fileName,varName):
    print("filename: {}".format(fileName))
    df = pd.read_csv(fileName)
    #print("----> {}".format(list(df.columns)))
    df[varName]=df[varName]/df["mc_weight"]
    weighted_yield = df[varName].sum()
    return weighted_yield

def get_ctauweighted_yield_from_csv(fileName,ctauVarName, varName):
    print("filename: {}".format(fileName))
    df = pd.read_csv(fileName)
    df[varName]=df[varName]/df["mc_weight"]
    df[varName] = df[varName]*df[ctauVarName]
    #print("----> {}".format(list(df.columns)))
    weighted_yield = df[varName].sum()
    return weighted_yield

def is_good_cand_var(varname):
    good = False
    if varname.find("C_")==0 and varname!="C_mu_Hnl_BS_ips" and varname!="C_mu_Hnl_BS_ip" and varname!="C_mu_Ds_BS_ips" and varname!="C_mu_Ds_BS_ip" and varname!="C_mu_Hnl_PV_ips" and varname!="C_mu_Hnl_PV_ip" and varname!="C_mu_Ds_PV_ips" and varname!="C_mu_Ds_PV_ip":
        good = True
    return good
