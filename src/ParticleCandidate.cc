//______________________________________________________________________________
// Filename: ParticleCandidate.cc
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Definitions of ParticleCandidate class members.
//______________________________________________________________________________

#include "ParticleCandidate.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default Contructor.
ParticleCandidate::ParticleCandidate()
    : particle_(new Particle()),
      cm_boost_(new Hep3Vector()),
      mdst_charged_(new Mdst_charged()),
      p_(new HepLorentzVector()),
      p_cm_(new HepLorentzVector()),
      track_(new TrackParameters())
{
  // Intentionally blank.
}

// Useful Constructor.
ParticleCandidate::ParticleCandidate(const Particle &particle,
    const Hep3Vector &cm_boost, const HepPoint3D &interaction_point)
    : particle_(new Particle(particle)),
      cm_boost_(new Hep3Vector(cm_boost)),
      mdst_charged_(new Mdst_charged(particle.mdstCharged())),
      p_(new HepLorentzVector(particle.p())),
      p_cm_(new HepLorentzVector(particle.p())),
      track_(new TrackParameters(particle, interaction_point))
{
  (*p_cm_).boost(cm_boost);
}

// Copy Constructor.
ParticleCandidate::ParticleCandidate(const ParticleCandidate &that)
    : particle_(new Particle(*that.particle_)),
      cm_boost_(new Hep3Vector(*that.cm_boost_)),
      mdst_charged_(new Mdst_charged(*that.mdst_charged_)),
      p_(new HepLorentzVector(*that.p_)),
      p_cm_(new HepLorentzVector(*that.p_cm_)),
      track_(new TrackParameters(*that.track_))
{
  // Intentionally blank.
}

// Assignment operator.
ParticleCandidate &ParticleCandidate::operator= (const ParticleCandidate &that)
{
  if (this != &that) {
    *particle_ = *(that.particle_);
    *cm_boost_ = *(that.cm_boost_);
    *mdst_charged_ = *(that.mdst_charged_);
    *p_ = *(that.p_);
    *p_cm_ = *(that.p_cm_);
    *track_ = *(that.track_);
  }
  return *this;
}

// Lepton Candidate Destructor.
ParticleCandidate::~ParticleCandidate()
{
  delete track_;
  delete p_cm_;
  delete p_;
  delete mdst_charged_;
  delete cm_boost_;
  delete particle_;
}

// Accessor for particle_.
Particle &ParticleCandidate::particle()
{
  return *particle_;
}

// Accessor for cm_boost_.
Hep3Vector &ParticleCandidate::cm_boost()
{
  return *cm_boost_;
}

// Accessor for mdst_charged_.
Mdst_charged &ParticleCandidate::mdstCharged()
{
  return *mdst_charged_;
}

// Accessor for p_.
HepLorentzVector &ParticleCandidate::p()
{
  return *p_;
}

// Accessor for p_cm_.
HepLorentzVector &ParticleCandidate::pCm()
{
  return *p_cm_;
}

// Accessor for track_.
TrackParameters &ParticleCandidate::track()
{
  return *track_;
}

// Returns the pythia particle ID code of the assignment given to particle
//   particle_ at its creation.
double
ParticleCandidate::idAssigned()
{
  return particle().pType().lund();
}

// Returns the mass hypothesis needed for the interaction point parameters.
int
ParticleCandidate::massHypothesis()
{
  int mass_hypothesis = -1;
  if (abs(idAssigned()) == 11) {
    mass_hypothesis = 0;
  } else if (abs(idAssigned()) == 13) {
    mass_hypothesis = 1;
  } else if (abs(idAssigned()) == 211) {
    mass_hypothesis = 2;
  } else if (abs(idAssigned()) == 321) {
    mass_hypothesis = 3;
  } else if (abs(idAssigned()) == 2212) {
    mass_hypothesis = 4;
  }
  return mass_hypothesis;
}

// Returns the pythia particle ID code of particle_ from the MC truth table.
//   Returns 0 if truth table is unavailable.
double
ParticleCandidate::idTrue()
{
  double id = 0;
  if (particle().relation().genHepevt()) {
    if (particle().relation().genHepevt().idhep()) {
      id = particle().relation().genHepevt().idhep();
    }
  }
  return id;
}

// Returns the pythia code of the particle mother as determined from the truth
//   table. Returns 0 if truth table is unavailable.
double
ParticleCandidate::idMom()
{
  double id = 0;
  if (particle().relation().genHepevt()) {
    if (particle().relation().genHepevt().mother()) {
      id = particle().relation().genHepevt().mother().idhep();
    }
  }
  return id;
}

// Returns the eid likelihood of particle_.
double
ParticleCandidate::electronProbability()
{
  eid eid_mdst(*mdst_charged_);
  return eid_mdst.prob(3, -1, 5);
}

// Returns the muid likelihood of particle_.
double
ParticleCandidate::muonProbability()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.Muon_likelihood();
}

// Returns the Chi^2 of the associated hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
double
ParticleCandidate::klmHitsChi2()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.Chi_2();
}

// Returns the number of the hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
int
ParticleCandidate::klmHits()
{
  Muid_mdst muid_mdst(*mdst_charged_);
  return muid_mdst.N_layer_hit_brl() + muid_mdst.N_layer_hit_end();
}

// Returns the Chi^2 of the associated hits in the KLM divided by the number of
//   the hits in the KLM assuming the track is a muon. If the track has no
//   associated hits in the KLM layers, the function returns 0 (consistant with 
//   an electron).
double
ParticleCandidate::klmChi2PerHits()
{
  if (klmHits() <= 0) {
    return 0;
  } else {
    return klmHitsChi2() / klmHits();
  }
}

// Returns the number of hits in the r-phi side of the SVD for particle_.
double
ParticleCandidate::svdRHits()
{
  int mass_hypothesis = massHypothesis();
  if (mass_hypothesis > -1) {
    return mdstCharged().trk().mhyp(mass_hypothesis).nhits(3);
  } else {
    return 0;
  }
}

// Returns the number of hits in the z-side of the SVD for particle_.
double
ParticleCandidate::svdZHits()
{
  int mass_hypothesis = massHypothesis();
  if (mass_hypothesis > -1) {
    return mdstCharged().trk().mhyp(mass_hypothesis).nhits(4);
  } else {
    return 0;
  }
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
