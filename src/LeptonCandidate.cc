//______________________________________________________________________________
// Filename: LeptonCandidate.cc
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Definitions of LeptonCandidate class members.
//______________________________________________________________________________

#include "LeptonCandidate.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default Contructor.
LeptonCandidate::LeptonCandidate()
{
  init();
}

// Lepton Candidate Constructor. 
LeptonCandidate::LeptonCandidate(const Particle &lepton,
    const Hep3Vector &cm_boost)
{
  init(lepton, cm_boost);
}

// Copy Constructor.
LeptonCandidate::LeptonCandidate(const LeptonCandidate &that)
{
  init(*that.lepton_, *that.cm_boost_);
}

// Assignment operator.
LeptonCandidate &LeptonCandidate::operator= (const LeptonCandidate &that)
{
  if (this == &that) {
    return *this;
  }
  dispose();
  init(*that.lepton_, *that.cm_boost_);
  return *this;
}

// Lepton Candidate Destructor.
LeptonCandidate::~LeptonCandidate()
{
  dispose();
}

void LeptonCandidate::init()
{
  lepton_ = 0;
  cm_boost_ = 0;
  mdst_charged_ = 0;
  p_ = 0;
  p_cm_ = 0;
}

void LeptonCandidate::init(const Particle &lepton, const Hep3Vector &cm_boost)
{
  lepton_ = new Particle(lepton);
  cm_boost_ = new Hep3Vector(cm_boost);
  mdst_charged_ = new Mdst_charged(lepton.relation().mdstCharged());
  p_ = new HepLorentzVector(lepton.p());
  p_cm_ = new HepLorentzVector(lepton.p());
  p_cm_->boost(cm_boost);
}

void LeptonCandidate::dispose() throw()
{
  if (lepton_) delete lepton_;
  if (cm_boost_) delete cm_boost_;
  if (mdst_charged_) delete mdst_charged_;
  if (p_) delete p_;
  if (p_cm_) delete p_cm_;
  lepton_ = 0;
  cm_boost_ = 0;
  mdst_charged_ = 0;
  p_ = 0;
  p_cm_ = 0;
}

// Accessor for lepton_.
Particle &LeptonCandidate::lepton()
{
  return *lepton_;
}

// Accessor for cm_boost_.
Hep3Vector &LeptonCandidate::cm_boost()
{
  return *cm_boost_;
}

// Accessor for mdst_charged_.
Mdst_charged &LeptonCandidate::mdstCharged()
{
  return *mdst_charged_;
}

// Accessor for p_.
HepLorentzVector &LeptonCandidate::p()
{
  return *p_;
}

// Accessor for p_cm_.
HepLorentzVector &LeptonCandidate::pCm()
{
  return *p_cm_;
}

// Returns the pythia particle ID code of the assignment given to particle
//   lepton_ at its creation.
double
LeptonCandidate::idAssigned()
{
  return lepton().pType().lund();
}

// Returns the mass hypothesis needed for the interaction point parameters.
// mass_hyp = 1 for muons, mass_hyp = 0 for electrons.
int
LeptonCandidate::massHypothesis()
{
  int mass_hyp = 0;
  if (abs(idAssigned()) == 13) {
     mass_hyp = 1;
  }
  return mass_hyp;
}

// Returns the pythia particle ID code of lepton_ as determined by the MC truth
//   table. Returns 0 if truth table is unavailable.
double
LeptonCandidate::idTrue()
{
  double id = 0;
  if (lepton().relation().genHepevt()) {
    if (lepton().relation().genHepevt().idhep()) {
      id = lepton().relation().genHepevt().idhep();
    }
  }
  return id;
}

// Returns the pythia code of the lepton mother as determined from the truth
//   table. Returns 0 if truth table is unavailable.
double
LeptonCandidate::idMom()
{
  double id = 0;
  if (lepton().relation().genHepevt()) {
    if (lepton().relation().genHepevt().mother()) {
      id = lepton().relation().genHepevt().mother().idhep();
    }
  }
  return id;
}

// Returns the eid likelihood of lepton_.
double
LeptonCandidate::electronProbability()
{
  eid leptonEid(mdstCharged());
  return leptonEid.prob(3, -1, 5);
}

// Returns the muid likelihood of lepton_.
double
LeptonCandidate::muonProbability()
{
  Muid_mdst leptonMuid(mdstCharged());
  return leptonMuid.Muon_likelihood();
}

// Returns the Chi^2 of the associated hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
double
LeptonCandidate::klmHitsChi2()
{
  Muid_mdst leptonMuid(mdstCharged());
  return leptonMuid.Chi_2();
}

// Returns the number of the hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
int
LeptonCandidate::klmHits()
{
  Muid_mdst leptonMuid(mdstCharged());
  return leptonMuid.N_layer_hit_brl() + leptonMuid.N_layer_hit_end();
}

// Returns the Chi^2 of the associated hits in the KLM divided by the number of
//   the hits in the KLM assuming the track is a muon. If the track has no
//   associated hits in the KLM layers, the function returns 0 (consistant with 
//   an electron).
double
LeptonCandidate::klmChi2PerHits()
{
  Muid_mdst leptonMuid(mdstCharged());
  if (klmHits() <= 0) {
    return 0;
  } else {
    return klmHitsChi2() / klmHits();
  }
}

// Returns the number of hits in the r-side of the SVD for lepton_.
double
LeptonCandidate::svdRadialHits()
{
  return mdstCharged().trk().mhyp(1).nhits(3);
}

// Returns the number of hits in the z-side of the SVD for lepton_.
double
LeptonCandidate::svdAxialHits()
{
  return mdstCharged().trk().mhyp(1).nhits(4);
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
