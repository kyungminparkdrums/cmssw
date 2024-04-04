#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // first create all the variables
  for (const auto &psvar : pset.getParameter<std::vector<edm::ParameterSet>>("variables")) {
    variables_.emplace_back(psvar.getParameter<std::string>("name"), psvar.getParameter<std::string>("value"));
  }
}


void l1tpf::HGC3DClusterID::evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) const {
  //Here we run the model and get the scores
  float puScore = 0; 
  float emScore = 0;
  float piScore = 0;

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
