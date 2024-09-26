#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // Inference of the conifer BDT model
  multiclass_bdt_ = new conifer::BDT<bdt_feature_t, bdt_score_t, false>(edm::FileInPath(pset.getParameter<std::string>("model")).fullPath());

  /*
  wp_PU_lowPt = 0.13562682; // 97% [5,10] GeV
  wp_PU_highPt = 0.36245897; // 90% [10, -] GeV
  wp_Pi_lowPt = 0.07440273; // 95%
  wp_Pi_highPt = 0.15240844; // 90%
  wp_Eg_lowPt = 0.24925998; // 96%
  wp_Eg_highPt = 0.3052651; // 99.3%
  */

  wp_PU = 0.37893048; // 91.3%
  wp_Pi = 0.13260031; // 98%
  wp_eg = 0.18857427; // 99%
}

float l1tpf::HGC3DClusterID::evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) {
  // Input for the BDT: showerlength, coreshowerlength, eot, eta, meanz, seetot, spptot, szz
  bdt_feature_t showerlength = cl.showerLength();
  bdt_feature_t coreshowerlength = cl.coreShowerLength();
  bdt_feature_t eot = cl.eot();
  bdt_feature_t eta = std::abs(cl.eta()); // take absolute values for eta for BDT input
  bdt_feature_t meanz = std::abs(cl.zBarycenter()) - 320;
  bdt_feature_t seetot = cl.sigmaEtaEtaTot();
  bdt_feature_t spptot = cl.sigmaPhiPhiTot();
  bdt_feature_t szz = cl.sigmaZZ();

  //bdt_feature_t srrtot = cl.sigmaRRTot();
  //bdt_feature_t emaxe = cl.eMax() / cl.energy();
  //bdt_feature_t emax3layers = cl.emax3layers(); 

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
  float maxScore = *std::max_element(bdt_score.begin(), bdt_score.end());

  cpf.setPuIDScore(puScore);
  cpf.setEmIDScore(emScore);
  cpf.setPiIDScore(piScore);
  
  cpf.setRho(cl.pt());

  return maxScore;
}


bool l1tpf::HGC3DClusterID::passPuID(l1t::PFCluster &cpf, float maxScore) {

  return (cpf.puIDScore() > wp_PU);
  
  //return (cpf.getRho() < 10 ? (cpf.puIDScore() > wp_PU_lowPt) : (cpf.puIDScore() > wp_PU_highPt));
}

bool l1tpf::HGC3DClusterID::passPFEmID(l1t::PFCluster &cpf, float maxScore) {

  return ((cpf.puIDScore() <= wp_PU) && (cpf.emIDScore() > wp_eg));

  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.emIDScore() > wp_Eg_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.emIDScore() > wp_Eg_highPt)));
}

bool l1tpf::HGC3DClusterID::passEgEmID(l1t::PFCluster &cpf, float maxScore) {
  
  return ((cpf.puIDScore() <= wp_PU) && (cpf.emIDScore() > wp_eg));

  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.emIDScore() > wp_Eg_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.emIDScore() > wp_Eg_highPt)));
}


bool l1tpf::HGC3DClusterID::passPiID(l1t::PFCluster &cpf, float maxScore) {
  
  return ((cpf.puIDScore() <= wp_PU) && (cpf.piIDScore() > wp_Pi));
  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.piIDScore() > wp_Pi_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.piIDScore() > wp_Pi_highPt)));
}
