
#include "SimG4Core/CustomPhysics/interface/CMSAntiSIMP.h"
#include "G4PhysicalConstants.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4ParticleTable.hh"

#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

CMSAntiSIMP* CMSAntiSIMP::theInstance = nullptr;

CMSAntiSIMP* CMSAntiSIMP::Definition(double mass) {
  if (theInstance != nullptr)
    return theInstance;
  const G4String name = "chibar";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance == nullptr) {
    // create particle
    //
    //    Arguments for constructor are as follows
    //               name             mass          width         charge
    //             2*spin           parity  C-conjugation
    //          2*Isospin       2*Isospin3       G-parity
    //               type    lepton number  baryon number   PDG encoding
    //             stable         lifetime    decay table
    //             shortlived      subType    anti_encoding

    anInstance = new G4ParticleDefinition(
        name, mass, 0, 0.0, 1, +1, 0, 0, -1, 0, "simp", 0, +1, -9000006, true, -1.0, nullptr, false, "nucleon");
  }
  theInstance = reinterpret_cast<CMSAntiSIMP*>(anInstance);
  return theInstance;
}

CMSAntiSIMP* CMSAntiSIMP::AntiSIMPDefinition(double mass) { return Definition(mass); }

CMSAntiSIMP* CMSAntiSIMP::AntiSIMP() {
  return Definition(1 * CLHEP::GeV);  // will use correct mass if instance exists
}
