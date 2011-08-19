//______________________________________________________________________________
// Filename: LeptonCandidate.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Class for storing and processing lepton candidates.
//______________________________________________________________________________

#ifndef LEPTONCANDIDATE_H
#define LEPTONCANDIDATE_H

#include "particle/Particle.h"       // The BELLE Particle Class.
#include "eid/eid.h"                 // For electron identification.
#include "mdst/Muid_mdst.h"           // For muon identification. 

#include HEPEVT_H                    // Panther.
#include MDST_H                      // Panther.

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class 
LeptonCandidate {
 protected:
  // Intrinsic constructor and destructor methods.
  LeptonCandidate();
  void init()
  void init(const Particle &lepton, const Hep3Vector &cm_boost);
  void dispose();
  
 private:
  LeptonCandidate(const LeptonCandidate &that);
  LeptonCandidate &operator= (const LeptonCandidate &that);
  
 public:
  LeptonCandidate(const Particle &lepton, const Hep3Vector &cm_boost);
  virtual ~LeptonCandidate();

  // Accessors.
  Particle &lepton();
  Hep3Vector &cm_boost();
  Mdst_charged &mdstCharged();
  HepLorentzVector &p();
  HepLorentzVector &pCm();

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
  double svdRadialHits();
  double svdAxialHits();

 private:
  // Existing objects.
  Particle *lepton_;
  Hep3Vector *cm_boost_;
  Mdst_charged *mdst_charged_;
  HepLorentzVector *p_;

  // New objects.
  HepLorentzVector *p_cm_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
