#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // Inference of the conifer BDT model
  multiclass_bdt_ = new conifer::BDT<bdt_feature_t, bdt_score_t, false>(edm::FileInPath(pset.getParameter<std::string>("model")).fullPath());
}

float l1tpf::HGC3DClusterID::evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) {
  // Input for the BDT: showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz
  bdt_feature_t showerlength = cl.showerLength();
  bdt_feature_t coreshowerlength = cl.coreShowerLength();
  bdt_feature_t eot = cl.eot();
  bdt_feature_t eta = std::abs(cl.eta()); // take absolute values for eta for BDT input
  bdt_feature_t meanz = cl.zBarycenter();
  bdt_feature_t seetot = cl.sigmaEtaEtaTot();
  bdt_feature_t spptot = cl.sigmaPhiPhiTot();
  bdt_feature_t szz = cl.sigmaZZ();

  // Run BDT inference
  inputs = {showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz};
  bdt_score = multiclass_bdt_->decision_function(inputs);

  // BDT score
  float puScore = bdt_score[0];
  float emScore = bdt_score[2];
  float piScore = bdt_score[1];

  // max score to ID the cluster
  float maxScore = *std::max_element(bdt_score.begin(), bdt_score.end());

  cpf.setPuIDScore(puScore);
  cpf.setEmIDScore(emScore);
  cpf.setPiIDScore(piScore);

  return maxScore;
}


bool l1tpf::HGC3DClusterID::passPuID(l1t::PFCluster &cpf, float maxScore) {
  // Using argmax 'WP' + and pass some 'minimal' WP on the max probability
  bool isMax = cpf.puIDScore() == maxScore;
  float puWP = -10; // dummy WP for now
  return isMax & (cpf.puIDScore() > puWP);
}

bool l1tpf::HGC3DClusterID::passPFEmID(l1t::PFCluster &cpf, float maxScore) {
  // Using argmax 'WP' + and pass some 'minimal' WP on the max probability
  bool isMax = cpf.emIDScore() == maxScore;
  float egWP = -10; // dummy one for now
  return isMax & (cpf.emIDScore() > egWP);
}

bool l1tpf::HGC3DClusterID::passEgEmID(l1t::PFCluster &cpf, float maxScore) {
  // Using argmax 'WP' + and pass some 'minimal' WP on the max probability
  bool isMax = cpf.emIDScore() == maxScore;
  float egWP = -10; // dummy one for now
  return isMax & (cpf.emIDScore() > egWP);
}


bool l1tpf::HGC3DClusterID::passPiID(l1t::PFCluster &cpf, float maxScore) {
  // Using argmax 'WP' + and pass some 'minimal' WP on the max probability
  bool isMax = cpf.piIDScore() == maxScore;
  float piWP = -10; // dummy one for now
  return isMax & (cpf.piIDScore() > piWP);
}
