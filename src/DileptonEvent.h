//______________________________________________________________________________
// Filename: DileptonEvent.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing information about dilepton events.
//______________________________________________________________________________

#ifndef DILEPTONEVENT_H
#define DILEPTONEVENT_H

#include <cmath>                     // Uses cmath functions.

#include "particle/Particle.h"       // The BELLE Particle Class.
#include "LeptonCandidate.h"

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif

class
DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent( Particle lepton0, Particle lepton1, Hep3Vector cm_boost );
  ~DileptonEvent() {}

  // Mutators.
  void set_l0( Particle lepton0 );
  void set_l1( Particle lepton1 );
  void set_cm_boost( Hep3Vector cm_boost );
  
  // Accessors.
  LeptonCandidate l0();
  LeptonCandidate l1();
  Hep3Vector cm_boost();
  
  // Methods
  double eventType();
  double cosThetaLL();
  double pSum();
  double pDifference();

 private:
  // Attributes.
  LeptonCandidate l0_;
  LeptonCandidate l1_;
  Hep3Vector cm_boost_;

};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
