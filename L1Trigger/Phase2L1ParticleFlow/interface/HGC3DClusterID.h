#ifndef L1Trigger_Phase2L1ParticleFlow_HGC3DClusterID_h
#define L1Trigger_Phase2L1ParticleFlow_HGC3DClusterID_h
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include <vector>
#include <cmath>

namespace l1tpf {
  class HGC3DClusterID {
  public:
    HGC3DClusterID(const edm::ParameterSet &pset);

    void evaluate(const l1t::HGCalMulticluster &cl, l1t::PFCluster &cpf) const;

    bool passPuID(l1t::PFCluster &cpf);
    bool passPFEmID(l1t::PFCluster &cpf);
    bool passEgEmID(l1t::PFCluster &cpf);
    bool passPiID(l1t::PFCluster &cpf);

  private:
    class Var {
    public:
      Var(const std::string &name, const std::string &expr) : name_(name), expr_(expr) {}
      // void declare(TMVA::Reader &r) { r.AddVariable(name_, &val_); }
      void fill(const l1t::HGCalMulticluster &c) { val_ = expr_(c); }

    private:
      std::string name_;
      StringObjectFunction<l1t::HGCalMulticluster> expr_;
      float val_;
    };

    std::vector<Var> variables_;


  };  //class
};    // namespace l1tpf

#endif
