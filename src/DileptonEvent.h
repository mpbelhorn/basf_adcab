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

#include "LeptonCandidate.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class
DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent(LeptonCandidate &lepton0, LeptonCandidate &lepton1);
  ~DileptonEvent() {}
  
  // Accessors.
  LeptonCandidate &l0() {return *l0_;}
  LeptonCandidate &l1() {return *l1_;}
  
  // Methods
  double eventType();
  double eventSign();
  double cosThetaLL();
  
 private:
  // Attributes.
  LeptonCandidate *l0_;
  LeptonCandidate *l1_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
