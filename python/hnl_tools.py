import argparse
import math
import sys


parser = argparse.ArgumentParser(description="")
parser.add_argument("hnl_mass"                 , help="Mass of the HNL [GeV]")
parser.add_argument("hnl_ctau"                 , help="ctau of the HNL [mm]")
parser.add_argument("--nDs"    , default=1000  , help="Number of reconstructed Ds")
parser.add_argument("--effDs"  , default=0.034 , help="Ds reconstruction efficiency")
parser.add_argument("--effHnl" , default=0.0142, help="Hnl reconstruction efficiency")
parser.add_argument("--fPrompt", default=0.552 , help="Fraction of reconstructed prompt Ds")
args = parser.parse_args()

m  = float(args.hnl_mass)
ct = float(args.hnl_ctau)

print("input mass: {} [GeV]".format(m))
print("input ctau: {} [mm]".format(ct))
print('\n')

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

#input parameters
eff_Ds   = float(args.effDs) # very preliminary value check MC matching
eff_Hnl  = float(args.effHnl) # very preliminary value computed with Bs->Ds->Hnl sample (m=1.5,ctau=10mm)
f_prompt = float(args.fPrompt) # very preliminary value
n_Ds     = float(args.nDs) # Ds found in RunA

# from eq A.11 https://arxiv.org/pdf/1805.08567.pdf
def compute_BR_DsToN (m_HNL,v2):
    square_brackets = (float(m_HNL)/Ds_mass)**2+(mu_mass/Ds_mass)**2-((float(m_HNL)/Ds_mass)**2-(mu_mass/Ds_mass)**2)**2
    square_root     = math.sqrt(1 + (float(m_HNL)/Ds_mass)**4 + (mu_mass/Ds_mass)**4 - 2*((float(m_HNL)/Ds_mass)**2) - 2*((mu_mass/Ds_mass)**2) -2*((float(m_HNL)/Ds_mass)**2)*((mu_mass/Ds_mass)**2))
    gamma           = (((G_f)**2)*(f_Ds**2)*(Ds_mass**3)*(V_cs**2)*v2*square_brackets*square_root)/(8*math.pi)
    gamma_Ds        = (hbar)/(Ds_meanLife)
    print("*** Computing BR(Ds-->NMu)")
    print("HNL mass: {}".format(m_HNL))
    print("|V|^2: {}".format(v2))
    print("BR(Ds-->NMu)  = {}".format(gamma/gamma_Ds))
    print('\n')
    return gamma/gamma_Ds

# ratio eq.2/eq.13 https://journals-aps-org.ezproxy.unibo.it/prd/pdf/10.1103/PhysRevD.94.113007
def compute_BR_NToPiMu (m_HNL):
    br = (96.0*(math.pi**2)*(f_pi**2)*(V_ud**2))/(16.0*10.95*(float(m_HNL)**2))
    print("*** Computing BR(N-->MuPi)")
    print("HNL mass: {}".format(m_HNL))
    print("BR(N-->MuPi) = {}".format(br))
    print('\n')
    return br

def compute_v2_from_ctau(m_HNL,ctau):
    num = light_speed*hbar
    den = ((G_f**2)*(m_HNL**5)*ctau)/(96*(math.pi**3))
    return num/den

if __name__ == "__main__":
    BR_DsToNMu = compute_BR_DsToN(m,compute_v2_from_ctau(m,ct))
    BR_NToPiMu = compute_BR_NToPiMu(m)
    BR_factor = (BR_DsToNMu*BR_NToPiMu)/(BR_DsToPhiPi*BR_PhiToMuMu)
    n_Hnl = n_Ds*f_prompt*BR_factor*(eff_Hnl/eff_Ds)
    print("*** Computing expected number of HNL")
    print("n_Ds = {}".format(n_Ds))
    print("f_prompt = {}".format(f_prompt))
    print("BR factor = (BR_DsToNMu*BR_NToPiMu)/(BR_DsToPhiPi*BR_PhiToMuMu) = {}".format(BR_factor))
    print("eff_Ds = {}".format(eff_Ds))
    print("eff_Hnl = {}".format(eff_Hnl))
    print("n_Hnl = n_Ds*f_prompt*BR_factor*(eff_Hnl/eff_Ds) = {}".format(n_Hnl))
