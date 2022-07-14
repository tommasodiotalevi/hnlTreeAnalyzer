#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;

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
