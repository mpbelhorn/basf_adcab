//______________________________________________________________________________
// Filename: ParticleCandidate.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Class for storing and processing charged particle candidates.
//______________________________________________________________________________

#ifndef PARTICLECANDIDATE_H
#define PARTICLECANDIDATE_H

#include "particle/Particle.h"       // The BELLE Particle Class.
#include "eid/eid.h"                 // For electron identification.
#include "mdst/Muid_mdst.h"          // For muon identification.
#include "TrackParameters.h"         // For checking track parameters.
#include HEPEVT_H                    // Panther.
#include MDST_H                      // Panther.

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class 
ParticleCandidate {
 public:
  // Default constructor.
  ParticleCandidate();

  // Copy constructor.
  ParticleCandidate(const ParticleCandidate &that);

  // Assignment operator.
  ParticleCandidate &operator= (const ParticleCandidate &that);

  // Useful constructor.
  ParticleCandidate(const Particle &particle, const Hep3Vector &cm_boost,
      const HepPoint3D &interaction_point);

  // Destructor.
  virtual ~ParticleCandidate();

  // Accessors.
  Particle &particle();
  Hep3Vector &cm_boost();
  Mdst_charged &mdstCharged();
  HepLorentzVector &p();
  HepLorentzVector &pCm();
  TrackParameters &track();

  // Methods.
  double idAssigned();
  double idTrue();
  double idMom();
  int massHypothesis();
  double electronProbability();
  double muonProbability();
  double klmHitsChi2();
  int klmHits();
  double klmChi2PerHits();
  double svdRHits();
  double svdZHits();

 private:
  // Existing objects.
  Particle *particle_;
  Hep3Vector *cm_boost_;
  Mdst_charged *mdst_charged_;
  HepLorentzVector *p_;

  // New objects.
  HepLorentzVector *p_cm_;
  TrackParameters *track_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
