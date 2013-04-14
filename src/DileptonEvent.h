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

#include "particle/Particle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class
DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent(Particle &lepton0, Particle &lepton1);
  ~DileptonEvent() {}

  // Accessors.
  Particle &l0() {return *l0_;}
  Particle &l1() {return *l1_;}

  // Methods
  double eventTypeAssigned();
  double eventTypeTrue();
  double eventSign();
  double cosThetaLL();
  double invariantMass();

 private:
  // Attributes.
  Particle *l0_;
  Particle *l1_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
