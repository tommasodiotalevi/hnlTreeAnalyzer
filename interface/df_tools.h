#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"

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
