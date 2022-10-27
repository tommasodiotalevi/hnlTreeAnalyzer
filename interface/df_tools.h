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

//select the best candidate based on difference between two variables,
RVec<short> and_3vec(RVec<short> var1, RVec<short> var2,RVec<short> var3)
{

  size_t n_cand = var1.size();
  RVec<Double_t> and_result(n_cand,0.);

  //select best candidate based on cand_var
  for (unsigned i=0; i<n_cand; ++i)
  {
    short res = var1[i]*var2[i]*var3[i];
    and_result[i] = res;
  }

  return and_result;
}

//select the best candidate based on difference between two variables,
size_t get_maxDiff_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2)
{
  float diff_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    float diff = std::abs(cand_var1.at(i)-cand_var2.at(i));
    if(diff>diff_best)
    {
      diff_best = diff;
      idx_best = i;
    }
  }

  return idx_best;
}

size_t get_maxRatio_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2)
{
  float ratio_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    float ratio = cand_var1.at(i)/cand_var2.at(i);
    if(ratio>ratio_best)
    {
      ratio_best = ratio;
      idx_best = i;
    }
  }

  return idx_best;
}

size_t get_maxWeightedRatio_wPosCond_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2, RVec<double> cand_var3,RVec<double> weight)
{
  float ratio_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    if(cand_var3.at(i)<0) continue;
    float ratio = (cand_var1.at(i)/cand_var2.at(i))*weight.at(i);
    if(ratio>ratio_best)
    {
      ratio_best = ratio;
      idx_best = i;
    }
  }

  return idx_best;
}

size_t get_maxRatio_wPosCond_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2, RVec<double> cand_var3)
{
  float ratio_best = -99999.;
  size_t idx_best = 9999;
  //size_t idx_best = 0;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    if(cand_var3.at(i)<0) continue;
    float ratio = cand_var1.at(i)/cand_var2.at(i);
    if(ratio>ratio_best)
    {
      ratio_best = ratio;
      idx_best = i;
    }
  }

  //if no candidates with positive condition has been found,
  //the condition is released
  if (idx_best==9999)
  {
    for (unsigned i=0; i<cand_var1.size(); ++i)
    {
      float ratio = cand_var1.at(i)/cand_var2.at(i);
      if(ratio>ratio_best)
      {
        ratio_best = ratio;
        idx_best = i;
      }
    }
  }

  return idx_best;
}

size_t get_cand_idx_custom(RVec<double> cand_var1, RVec<double> cand_var2,RVec<double> cand_var3,RVec<double> cond,RVec<double> cand_var4,RVec<double> cand_var5)
{
  float sum_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    if (cond.at(i)<0) continue;
    float sum = (cand_var1.at(i)+cand_var2.at(i)+cand_var3.at(i))*(cand_var4.at(i)/cand_var5.at(i));
    if(sum>sum_best)
    {
      sum_best = sum;
      idx_best = i;
    }
  }

  return idx_best;
}

size_t get_maxSum_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2,RVec<double> cand_var3)
{
  float sum_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    float sum = cand_var1.at(i)+cand_var2.at(i)+cand_var3.at(i);
    if(sum>sum_best)
    {
      sum_best = sum;
      idx_best = i;
    }
  }

  return idx_best;
}

size_t get_maxSum_cand_idx(RVec<double> cand_var1, RVec<double> cand_var2)
{
  float sum_best = -99999.;
  size_t idx_best = 9999;

  //select best candidate based on cand_var
  for (unsigned i=0; i<cand_var1.size(); ++i)
  {
    float sum = std::abs(cand_var1.at(i)+cand_var2.at(i));
    if(sum>sum_best)
    {
      sum_best = sum;
      idx_best = i;
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

RVec<short> get_mu_trigger_matching(RVec<short> trigger_match, RVec<float> trigger_dr, RVec<unsigned> trig_mu_idx, RVec<unsigned> mu_idx)
{
  RVec<short> mu_has_matched_trigger(mu_idx.size());

  for(unsigned i_trigMatch=0; i_trigMatch<trigger_match.size(); ++i_trigMatch){

    if (!(trigger_match.at(i_trigMatch)>0)) continue;
    if (trigger_dr.at(i_trigMatch)>0.005) continue;

    for(unsigned i_mu=0; i_mu<mu_idx.size(); ++i_mu){
      if (mu_idx[i_mu] == trig_mu_idx[i_trigMatch])
        mu_has_matched_trigger[i_mu] = 1;
    }
  }

  return mu_has_matched_trigger;
}

RVec<float> get_mu_trigger_dr(RVec<short> trigger_match, RVec<float> trigger_dr, RVec<unsigned> trig_mu_idx, RVec<unsigned> mu_idx)
{
  RVec<float> muTrig_dr(mu_idx.size(),9999.);

  for(unsigned i_trigMatch=0; i_trigMatch<trigger_match.size(); ++i_trigMatch){

    if (!(trigger_match.at(i_trigMatch)>0)) continue;

    for(unsigned i_mu=0; i_mu<mu_idx.size(); ++i_mu){
      if (mu_idx[i_mu] == trig_mu_idx[i_trigMatch])
        muTrig_dr[i_mu] = trigger_dr[i_trigMatch];
    }
  }

  return muTrig_dr;
}

//this is to get candidate invariant mass with different 
//final state particle mass ipothesis
RVec<Double_t> get_PxPyPzE_newM_2PartInvMass(RVec<Double_t> px_1, RVec<Double_t> py_1, RVec<Double_t> pz_1, RVec<Double_t> E_1, Double_t m_1,
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

RVec<Double_t> get_PtEtaPhi_newM_3PartInvMass(RVec<Double_t> pt_1, RVec<Double_t> eta_1, RVec<Double_t> phi_1, Double_t newM_1,
		                              RVec<Double_t> pt_2, RVec<Double_t> eta_2, RVec<Double_t> phi_2, Double_t newM_2,
		                              RVec<Double_t> pt_3, RVec<Double_t> eta_3, RVec<Double_t> phi_3, Double_t newM_3)
{
  size_t n_cand = pt_1.size();
  RVec<Double_t> cand_inv_mass(n_cand,0.);

  for(unsigned i=0; i<n_cand; ++i)
  {
    TLorentzVector p4_1;
    TLorentzVector p4_2;
    TLorentzVector p4_3;
    p4_1.SetPtEtaPhiM(pt_1[i],eta_1[i],phi_1[i],newM_1);
    p4_2.SetPtEtaPhiM(pt_2[i],eta_2[i],phi_2[i],newM_2);
    p4_3.SetPtEtaPhiM(pt_3[i],eta_3[i],phi_3[i],newM_3);
    cand_inv_mass[i]=(p4_1+p4_2+p4_3).M();
  }

  return cand_inv_mass;

}

RVec<Double_t> get_PtEtaPhi_newM_2PartInvMass(RVec<Double_t> pt_1, RVec<Double_t> eta_1, RVec<Double_t> phi_1, Double_t newM_1,
		                              RVec<Double_t> pt_2, RVec<Double_t> eta_2, RVec<Double_t> phi_2, Double_t newM_2)
{
  size_t n_cand = pt_1.size();
  RVec<Double_t> cand_inv_mass(n_cand,0.);

  for(unsigned i=0; i<n_cand; ++i)
  {
    TLorentzVector p4_1;
    TLorentzVector p4_2;
    p4_1.SetPtEtaPhiM(pt_1[i],eta_1[i],phi_1[i],newM_1);
    p4_2.SetPtEtaPhiM(pt_2[i],eta_2[i],phi_2[i],newM_2);
    cand_inv_mass[i]=(p4_1+p4_2).M();
  }

  return cand_inv_mass;

}

RVec<Double_t> get_PxPyPzE_newM_3PartInvMass(RVec<Double_t> px_1, RVec<Double_t> py_1, RVec<Double_t> pz_1, RVec<Double_t> E_1, Double_t m_1,
		                                 RVec<Double_t> px_2, RVec<Double_t> py_2, RVec<Double_t> pz_2, RVec<Double_t> E_2, Double_t m_2,
		                                 RVec<Double_t> px_3, RVec<Double_t> py_3, RVec<Double_t> pz_3, RVec<Double_t> E_3, Double_t m_3)
{
  size_t n_cand = px_1.size();
  RVec<Double_t> cand_inv_mass(n_cand,0.);

  for(unsigned i=0; i<n_cand; ++i)
  {
    TLorentzVector p4_1;
    TLorentzVector p4_2;
    TLorentzVector p4_3;
    p4_1.SetPxPyPzE(px_1[i],py_1[i],pz_1[i],E_1[i]);
    p4_2.SetPxPyPzE(px_2[i],py_2[i],pz_2[i],E_2[i]);
    p4_3.SetPxPyPzE(px_3[i],py_3[i],pz_3[i],E_3[i]);
    TLorentzVector p4_newM_1;
    TLorentzVector p4_newM_2;
    TLorentzVector p4_newM_3;
    p4_newM_1.SetPtEtaPhiM(p4_1.Pt(),p4_1.Eta(),p4_1.Phi(),m_1);
    p4_newM_2.SetPtEtaPhiM(p4_2.Pt(),p4_2.Eta(),p4_2.Phi(),m_2);
    p4_newM_3.SetPtEtaPhiM(p4_3.Pt(),p4_3.Eta(),p4_3.Phi(),m_3);
    cand_inv_mass[i]=(p4_newM_1+p4_newM_2+p4_newM_3).M();
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
