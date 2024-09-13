#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // Inference of the conifer BDT model
  multiclass_bdt_ = new conifer::BDT<bdt_feature_t, bdt_score_t, false>(edm::FileInPath(pset.getParameter<std::string>("model")).fullPath());
  wp_pu = pset.getParameter<double>("wp_pu");
  wp_pion = pset.getParameter<double>("wp_pion");
  wp_eg_loose = pset.getParameter<double>("wp_eg_loose");
  wp_eg_tight = pset.getParameter<double>("wp_eg_tight");
}

void l1tpf::HGC3DClusterID::evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) {
  // Input for the BDT: showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz
  bdt_feature_t showerlength = cl.showerLength();
  bdt_feature_t coreshowerlength = cl.coreShowerLength();
  bdt_feature_t eot = cl.eot();
  bdt_feature_t eta = std::abs(cl.eta()); // take absolute values for eta for BDT input
  bdt_feature_t meanz = std::abs(cl.zBarycenter()) - 320;
  bdt_feature_t seetot = cl.sigmaEtaEtaTot();
  bdt_feature_t spptot = cl.sigmaPhiPhiTot();
  bdt_feature_t szz = cl.sigmaZZ();

  // Run BDT inference
  inputs = {showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz};

  bdt_score = multiclass_bdt_->decision_function(inputs);

  // BDT score
  //float puScore = bdt_score[0];
  //float emScore = bdt_score[2];
  //float piScore = bdt_score[1];
  
  float puRawScore = bdt_score[0];
  float emRawScore = bdt_score[2];
  float piRawScore = bdt_score[1];

  // softmax (for now, let's compute the softmax in this code; this needs to be changed to implement on firmware)
  // Softmax implemented in conifer (standalone) is to be integrated here soon; for now, just do "offline" softmax :(
  float denom = exp(puRawScore) + exp(emRawScore) + exp(piRawScore);
  float puScore = exp(puRawScore) / denom;
  float emScore = exp(emRawScore) / denom;
  float piScore = exp(piRawScore) / denom;

  // max score to ID the cluster -> Deprecated
  //float maxScore = *std::max_element(bdt_score.begin(), bdt_score.end());

  cpf.setPuIDScore(puScore);
  cpf.setEmIDScore(emScore);
  cpf.setPiIDScore(piScore);

  //return maxScore;
}


bool l1tpf::HGC3DClusterID::passPuID(l1t::PFCluster &cpf) {
  return (cpf.puIDScore() > wp_pu);
}

bool l1tpf::HGC3DClusterID::passPFEmID(l1t::PFCluster &cpf) {
  return ((cpf.puIDScore() <= wp_pu) && (cpf.emIDScore() > wp_eg_loose));
}

bool l1tpf::HGC3DClusterID::passEgEmIDLoose(l1t::PFCluster &cpf) {
  return ((cpf.puIDScore() <= wp_pu) && (cpf.emIDScore() > wp_eg_loose));
}

bool l1tpf::HGC3DClusterID::passEgEmIDTight(l1t::PFCluster &cpf) {
  return ((cpf.puIDScore() <= wp_pu) && (cpf.emIDScore() > wp_eg_tight));
}

bool l1tpf::HGC3DClusterID::passPiID(l1t::PFCluster &cpf) {
  return ((cpf.puIDScore() <= wp_pu) && (cpf.piIDScore() > wp_pion));
}
