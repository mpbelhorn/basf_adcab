//______________________________________________________________________________
// Filename: DileptonEvent.cc
// Version: 2011.08.11.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing information about dilepton events.
//______________________________________________________________________________

#include "DileptonEvent.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default constructor.
DileptonEvent::DileptonEvent()
{
  // Function is blank.
}

// Dilepton constructor.
DileptonEvent::DileptonEvent(ParticleCandidate &lepton0,
                             ParticleCandidate &lepton1)
{
  l0_ = &lepton0;
  l1_ = &lepton1;
}

// Determines the type of dilepton event and returns an integer for the
//   cases of 0: Unidentified/Error - includes false lepton events,
//            1: Dielectron,
//            2: Dimuon,
//            3: Electron-muon.
// The Lund values as of 2010/08/12 are e(+/-) = (+/-)11, mu(+/-) = (+/-)13.
double
DileptonEvent::eventTypeAssigned()
{
  double lund_sum = abs(l0().idAssigned()) + abs(l1().idAssigned());
  double type = 0;
  if (lund_sum == 22) {
    type = 1;
  } else if (lund_sum == 24) {
    type = 2;
  } else if (lund_sum == 26) {
    type = 3;
  }
  return type;
}

// Same as eventTypeAssigned except uses truth-table 
double
DileptonEvent::eventTypeTrue()
{
  double lund_sum = abs(l0().idTrue()) + abs(l1().idTrue());
  double type = 0;
  if (lund_sum == 22) {
    type = 1;
  } else if (lund_sum == 24) {
    type = 2;
  } else if (lund_sum == 26) {
    type = 3;
  }
  return type;
}

// Returns the sign of the dilepton pair. -1 == --, 0 == +-/-+, and 1 == ++.
double
DileptonEvent::eventSign()
{
  return (l0().particle().charge() + l1().particle().charge()) / 2;
}

// Calculates the cosine of the opening angle theta_ll between the two leptons
//   in the CM frame.
double
DileptonEvent::cosThetaLL()
{
  Hep3Vector lepton0_p3_cm = l0().pCm().vect();
  Hep3Vector lepton1_p3_cm = l1().pCm().vect();
  
  double cosine_theta_ll;
  double total_p_squared = lepton0_p3_cm.mag2() * lepton1_p3_cm.mag2();
  if (total_p_squared <= 0) {
    // This case can only happen if a lepton 3-momenta is not real or the 0
    // vector. Thus, I set the cosine to a value that will be cut.
    cosine_theta_ll = 1.0;
  } else {
    cosine_theta_ll = lepton0_p3_cm.dot(lepton1_p3_cm) / sqrt(total_p_squared);
    if (cosine_theta_ll >  1.0) cosine_theta_ll =  1.0;
    if (cosine_theta_ll < -1.0) cosine_theta_ll = -1.0;
  }
  return cosine_theta_ll;
}

// Calculates the invariant mass of the two leptons.
double
DileptonEvent::invariantMass()
{
  HepLorentzVector l0_p_cm = l0().p();
  HepLorentzVector l1_p_cm = l1().p();

  return (l0_p_cm + l1_p_cm).m();
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
