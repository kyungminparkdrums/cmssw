#include "L1Trigger/Phase2L1ParticleFlow/interface/HGC3DClusterID.h"

l1tpf::HGC3DClusterID::HGC3DClusterID(const edm::ParameterSet &pset) {
  // Inference of the conifer BDT model
  multiclass_bdt_ = new conifer::BDT<bdt_feature_t, bdt_score_t, false>(edm::FileInPath(pset.getParameter<std::string>("model")).fullPath());

  // First pT-depedent WPs
  //wp_PU = {0.31114486, 0.3014869, 0.19625592};
  //wp_Pi = {0.13330694, 0.063722104, 0.04143001};
  //wp_Eg = {0.116029575, 0.061485082, 0.052452337};

  // More tuned version; lower rate at higher pT bin, but also lower eff :( 
  //wp_PU = {0.31114486, 0.3014869, 0.25903508};
  //wp_Pi = {0.13330694, 0.17924163, 0.10514013};
  //wp_Eg = {0.116029575, 0.21594438, 0.5486468};
  
  // Two pT-bin WPs: bin boundary of 20GeV
  wp_PU = {0.33752286, 0.28073967};
  wp_Pi = {0.1621892, 0.10513939};
  wp_Eg = {0.12471038, 0.3620104};
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

  bdt_feature_t emaxe = cl.eMax() / cl.energy();
  bdt_feature_t emax3layers = cl.emax3layers(); 
  bdt_feature_t first5layers = cl.first5layers(); 
  bdt_feature_t ntc67 = cl.triggerCells67percent(); 

  // Run BDT inference
  inputs = {eta, coreshowerlength, first5layers, emax3layers, emaxe, ntc67, seetot, spptot};

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
  if (cpf.getRho() < 20) { return (cpf.puIDScore() > wp_PU[0]); }
  else { return (cpf.puIDScore() > wp_PU[1]); }

  //return (cpf.getRho() < 10 ? (cpf.puIDScore() > wp_PU_lowPt) : (cpf.puIDScore() > wp_PU_highPt));
}

bool l1tpf::HGC3DClusterID::passPFEmID(l1t::PFCluster &cpf, float maxScore) {
  if (cpf.getRho() < 20) { return ((cpf.puIDScore() <= wp_PU[0]) && (cpf.emIDScore() > wp_Eg[0])); }
  else { return ((cpf.puIDScore() <= wp_PU[1]) && (cpf.emIDScore() > wp_Eg[1])); }

  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.emIDScore() > wp_Eg_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.emIDScore() > wp_Eg_highPt)));
}

bool l1tpf::HGC3DClusterID::passEgEmID(l1t::PFCluster &cpf, float maxScore) {
  if (cpf.getRho() < 20) { return ((cpf.puIDScore() <= wp_PU[0]) && (cpf.emIDScore() > wp_Eg[0])); }
  else { return ((cpf.puIDScore() <= wp_PU[1]) && (cpf.emIDScore() > wp_Eg[1])); }

  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.emIDScore() > wp_Eg_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.emIDScore() > wp_Eg_highPt)));
}


bool l1tpf::HGC3DClusterID::passPiID(l1t::PFCluster &cpf, float maxScore) {
  if (cpf.getRho() < 20) { return ((cpf.puIDScore() <= wp_PU[0]) && (cpf.emIDScore() > wp_Pi[0])); }
  else { return ((cpf.puIDScore() <= wp_PU[1]) && (cpf.emIDScore() > wp_Pi[1])); }
 
  //return (cpf.getRho() < 10 ? ((cpf.puIDScore() <= wp_PU_lowPt) && (cpf.piIDScore() > wp_Pi_lowPt)) : ((cpf.puIDScore() <= wp_PU_highPt) && (cpf.piIDScore() > wp_Pi_highPt)));
}
