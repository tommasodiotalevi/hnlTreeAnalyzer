#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLorentzVector.h"

using namespace ROOT;
using namespace ROOT::VecOps;

//select the best candidate based on cand_var,
size_t get_best_cand_idx(RVec<double> cand_var)
{
  float var_best = -1;
  size_t idx_best = 9999;
  for (unsigned i=0; i<cand_var.size(); ++i)
  {
    if(cand_var.at(i)>var_best)
    {
      var_best = cand_var.at(i);
      idx_best = i;
    }
  }

  return idx_best;
}

//select the best candidate based on cand_var,
//if more than one mu1mu2pi permutation is present
//then then select the best combination based on mu1mu2_var
size_t get_best_cand_idx(RVec<double> cand_var, RVec<double> mu1mu2_var, RVec<unsigned> mu1_idx, RVec<unsigned> mu2_idx, RVec<unsigned> pi_idx)
{
  float var_best = -9999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var.size(); ++i)
  {
    if(cand_var.at(i)>var_best)
    {
      var_best = cand_var.at(i);
      idx_best = i;
    }
  }

  //if more than one permutation exists, check mu1m2_var
  for (unsigned j=0; j<cand_var.size(); ++j)
  {
    if (mu1_idx.at(idx_best) == mu2_idx.at(j) && mu2_idx.at(idx_best) == mu1_idx.at(j) && pi_idx.at(idx_best) == pi_idx.at(j))
    {
      if(mu1mu2_var.at(j) > mu1mu2_var.at(idx_best))
      {
        idx_best = j;
      }
    }
  }

  return idx_best;
}

RVec<short> get_mu_trigger_matching(RVec<short> trigger_match, RVec<unsigned> trig_mu_idx, RVec<unsigned> mu_idx)
{
  RVec<short> mu_has_matched_trigger(mu_idx.size());

  for(unsigned i_trigMatch=0; i_trigMatch<trigger_match.size(); ++i_trigMatch){

    if (!(trigger_match.at(i_trigMatch)>0)) continue;

    for(unsigned i_mu=0; i_mu<mu_idx.size(); ++i_mu){
      if (mu_idx[i_mu] == trig_mu_idx[i_trigMatch])
        mu_has_matched_trigger[i_mu] = 1;
    }
  }

  return mu_has_matched_trigger;
}

//this is to get candidate invariant mass with different 
//final state particle mass ipothesis
RVec<Double_t> get_PxPyPzE_newM_HnlInvariantMass(RVec<Double_t> px_1, RVec<Double_t> py_1, RVec<Double_t> pz_1, RVec<Double_t> E_1, Double_t m_1,
		                                 RVec<Double_t> px_2, RVec<Double_t> py_2, RVec<Double_t> pz_2, RVec<Double_t> E_2, Double_t m_2)
{
  size_t n_cand = px_1.size();
  RVec<Double_t> cand_inv_mass(n_cand,0.);

  for(unsigned i=0; i<n_cand; ++i)
  {
    TLorentzVector p4_1;
    TLorentzVector p4_2;
    p4_1.SetPxPyPzE(px_1[i],py_1[i],pz_1[i],E_1[i]);
    p4_2.SetPxPyPzE(px_2[i],py_2[i],pz_2[i],E_2[i]);
    TLorentzVector p4_newM_1;
    TLorentzVector p4_newM_2;
    p4_newM_1.SetPtEtaPhiM(p4_1.Pt(),p4_1.Eta(),p4_1.Phi(),m_1);
    p4_newM_2.SetPtEtaPhiM(p4_2.Pt(),p4_2.Eta(),p4_2.Phi(),m_2);
    cand_inv_mass[i]=(p4_newM_1+p4_newM_2).M();
  }

  return cand_inv_mass;

}

RVec<Double_t> get_PxPyPzE_HnlInvariantMass(RVec<Double_t> px_1, RVec<Double_t> py_1, RVec<Double_t> pz_1, RVec<Double_t> E_1,
		                            RVec<Double_t> px_2, RVec<Double_t> py_2, RVec<Double_t> pz_2, RVec<Double_t> E_2)
{
  size_t n_cand = px_1.size();
  RVec<Double_t> cand_inv_mass(n_cand,0.);

  for(unsigned i=0; i<n_cand; ++i)
  {
    TLorentzVector p4_1;
    TLorentzVector p4_2;
    p4_1.SetPxPyPzE(px_1[i],py_1[i],pz_1[i],E_1[i]);
    p4_2.SetPxPyPzE(px_2[i],py_2[i],pz_2[i],E_2[i]);
    cand_inv_mass[i]=(p4_1+p4_2).M();
  }

  return cand_inv_mass;

}

//for debugging
RVec<short> get_multiple_permutations(RVec<unsigned> mu1_idx, RVec<unsigned> mu2_idx, RVec<unsigned> pi_idx)
{
  RVec<short> cand_has_multiple_permutations(mu1_idx.size(),0);

  for (unsigned i=0; i<mu1_idx.size(); ++i)
  {
    for (unsigned j=+1; j<mu1_idx.size(); ++j)
    {
      if (mu1_idx.at(i) == mu2_idx.at(j) && mu2_idx.at(i) == mu1_idx.at(j) && pi_idx.at(i) == pi_idx.at(j))
      {
        cand_has_multiple_permutations[i] = 1;
        cand_has_multiple_permutations[j] = 1;
      }
    }
  }

  return cand_has_multiple_permutations;

}
