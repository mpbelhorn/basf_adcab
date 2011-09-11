//______________________________________________________________________________
// Filename: DileptonEvent.h
// Version: 2011.08.11.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing information about dilepton events.
//______________________________________________________________________________

#ifndef DILEPTONEVENT_H
#define DILEPTONEVENT_H

#include <cmath>                     // Uses cmath functions.

#include "ParticleCandidate.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class
DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent(ParticleCandidate &lepton0, ParticleCandidate &lepton1);
  ~DileptonEvent() {}
  
  // Accessors.
  ParticleCandidate &l0() {return *l0_;}
  ParticleCandidate &l1() {return *l1_;}
  
  // Methods
  double eventTypeAssigned();
  double eventTypeTrue();
  double eventSign();
  double cosThetaLL();
  double invariantMass();
  
 private:
  // Attributes.
  ParticleCandidate *l0_;
  ParticleCandidate *l1_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
