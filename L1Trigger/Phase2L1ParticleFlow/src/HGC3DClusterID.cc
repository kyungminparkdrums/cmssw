#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // first create all the variables
  for (const auto &psvar : pset.getParameter<std::vector<edm::ParameterSet>>("variables")) {
    variables_.emplace_back(psvar.getParameter<std::string>("name"), psvar.getParameter<std::string>("value"));
  }
}

void l1tpf::HGC3DClusterID::evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) const {
  // Inference of the conifer BDT model
  multiclass_bdt_ = new conifer::BDT<bdt_feature_t, bdt_score_t, false>("../data/multiclassID/conifer_bridge_1712763069.so");

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
  std::vector<bdt_feature_t> inputs = {showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz};
  std::vector<bdt_score_t> bdt_score = multiclass_bdt_->decision_function(inputs);

  // BDT score
  float puScore = bdt_score[0];
  float emScore = bdt_score[2];
  float piScore = bdt_score[1];

  cpf.setPuIDScore(puScore);
  cpf.setEmIDScore(emScore);
  cpf.setPiIDScore(piScore);
}


bool l1tpf::HGC3DClusterID::passPuID(l1t::PFCluster &cpf) {
  // here we evaluate the WPs
  return cpf.puIDScore() > -1;
}

bool l1tpf::HGC3DClusterID::passPFEmID(l1t::PFCluster &cpf) {
  // here we evaluate the WPs
  return cpf.emIDScore() > -1;
}

bool l1tpf::HGC3DClusterID::passEgEmID(l1t::PFCluster &cpf) {
  // here we evaluate the WPs
  return cpf.emIDScore() > -1;
}


bool l1tpf::HGC3DClusterID::passPiID(l1t::PFCluster &cpf) {
  // here we evaluate the WPs
  return cpf.piIDScore() > -1;
}
