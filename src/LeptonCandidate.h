//______________________________________________________________________________
// Filename: LeptonCandidate.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Class for storing and processing lepton candidates.
//______________________________________________________________________________

#ifndef LEPTONCANDIDATE_H
#define LEPTONCANDIDATE_H

#include <cmath>                     // Uses cmath functions.

#include "belle.h"                   // BELLE Library.
#include "particle/Particle.h"       // The BELLE Particle Class.
#include "kid/atc_pid.h"             // For particle species separation.
#include "eid/eid.h"                 // For electron identification.
#include "mdst/mdst.h"               // For MDST files.
#include "mdst/Muid_mdst.h"           // For muon identification. 
#include "ip/IpProfile.h"            // Beam Interaction Point (IP) analysis
                                     //   tools. Position unit = cm.

#include <panther/panther.h>         // Panther.
#include BELLETDF_H                  // Panther.
#include HEPEVT_H                    // Panther.
#include MDST_H                      // Panther.

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif

class 
LeptonCandidate {
 public:

  // Constructors and destructor.
  LeptonCandidate();
  LeptonCandidate( Particle lepton, Hep3Vector cm_boost );
  ~LeptonCandidate() {};

  // Mutators.
  void set_lepton( Particle lepton );
  void set_cm_boost_vector( Hep3Vector cm_boost );

  // Accessors.
  Particle lepton();
  Hep3Vector cm_boost();

  // Methods.
  double idAssigned();
  double idTrue();
  double idMom();
  int massHypothesis();
  HepLorentzVector pCm();
  HepLorentzVector p();
  double electronProbability();
  double muonProbability();
  double klmHitsChi2();
  int klmHits();
  double klmChi2PerHits();
  double svdRadialHits();
  double svdAxialHits();

 private:

  // Attributes.
  Particle lepton_;
  Hep3Vector cm_boost_;

};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
