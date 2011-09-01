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
    : lepton_(new Particle()),
      cm_boost_(new Hep3Vector()),
      mdst_charged_(new Mdst_charged()),
      p_(new HepLorentzVector()),
      p_cm_(new HepLorentzVector())
{
  // Intentionally blank.
}

// Lepton Candidate Constructor. 
LeptonCandidate::LeptonCandidate(const Particle &lepton,
    const Hep3Vector &cm_boost)
    : lepton_(new Particle(lepton)),
      cm_boost_(new Hep3Vector(cm_boost)),
      mdst_charged_(new Mdst_charged(lepton.mdstCharged())),
      p_(new HepLorentzVector(lepton.p())),
      p_cm_(new HepLorentzVector(lepton.p()))
{
  (*p_cm_).boost(cm_boost);
}

// Copy Constructor.
LeptonCandidate::LeptonCandidate(const LeptonCandidate &that)
    : lepton_(new Particle(*that.lepton_)),
      cm_boost_(new Hep3Vector(*that.cm_boost_)),
      mdst_charged_(new Mdst_charged(*that.mdst_charged_)),
      p_(new HepLorentzVector(*that.p_)),
      p_cm_(new HepLorentzVector(*that.p_cm_))
{
  // Intentionally blank.
}

// Assignment operator.
LeptonCandidate &LeptonCandidate::operator= (const LeptonCandidate &that)
{
  if (this != &that) {
    *lepton_ = *(that.lepton_);
    *cm_boost_ = *(that.cm_boost_);
    *mdst_charged_ = *(that.mdst_charged_);
    *p_ = *(that.p_);
    *p_cm_ = *(that.p_cm_);
  }
  return *this;
}

// Lepton Candidate Destructor.
LeptonCandidate::~LeptonCandidate()
{
  delete p_cm_;
  delete p_;
  delete mdst_charged_;
  delete cm_boost_;
  delete lepton_;
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
  int mass_hyp = -1;
  if (abs(idAssigned()) == 11) {
    mass_hyp = 0;
  } else if (abs(idAssigned()) == 13) {
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
  eid eid_mdst(*mdst_charged_);
  return eid_mdst.prob(3, -1, 5);
}

// Returns the muid likelihood of lepton_.
double
LeptonCandidate::muonProbability()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.Muon_likelihood();
}

// Returns the Chi^2 of the associated hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
double
LeptonCandidate::klmHitsChi2()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.Chi_2();
}

// Returns the number of the hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
int
LeptonCandidate::klmHits()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.N_layer_hit_brl() + muid_mdst.N_layer_hit_end();
}

// Returns the Chi^2 of the associated hits in the KLM divided by the number of
//   the hits in the KLM assuming the track is a muon. If the track has no
//   associated hits in the KLM layers, the function returns 0 (consistant with 
//   an electron).
double
LeptonCandidate::klmChi2PerHits()
{
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
