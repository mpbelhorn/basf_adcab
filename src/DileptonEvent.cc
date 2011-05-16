//______________________________________________________________________________
// Filename: DileptonEvent.cc
// Version: 2011.05.15.A
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
DileptonEvent::DileptonEvent( Particle lepton0, Particle lepton1,
    Hep3Vector cm_boost )
{
  l0_ = LeptonCandidate( lepton0, cm_boost );
  l1_ = LeptonCandidate( lepton1, cm_boost );
  cm_boost_ = cm_boost;
}

// Mutator for l0_.
void
DileptonEvent::set_l0( Particle lepton0 )
{
  l0_ = LeptonCandidate( lepton0, cm_boost_ );
}

// Mutator for l1_.
void
DileptonEvent::set_l1( Particle lepton1 )
{
  l1_ = LeptonCandidate( lepton1, cm_boost_ );
}

// Mutator for cm_boost_.
void
DileptonEvent::set_cm_boost( Hep3Vector cm_boost )
{
  cm_boost_ = cm_boost;
}

// Accessor for l0_.
LeptonCandidate
DileptonEvent::l0()
{
  return l0_;
}

// Accessor for l1_.
LeptonCandidate
DileptonEvent::l1()
{
  return l1_;
}

// Determines the type of dilepton event and returns an integer for the
//   cases of 22:(e+e+)/(e-e-);     -22:e+e-;
//            26:(mu+e+)/(mu-e-);   -26:mu+mu-
//            24:(mu+mu+)/(mu-mu-); -24:(mu+e-)/(mu-e+).
// Negative codes imply opposite-sign whereas positive codes imply same-sign.
// The magnitude of the code is the sum of the magnitudes of the Lund ID values.
// The Lund values as of 2010/08/12 are e(+/-) = (+/-)11, mu(+/-) = (+/-)13.
double
DileptonEvent::eventType()
{
  double sign = l0_.lepton().charge() * l1_.lepton().charge();
  double lundSum = abs( l0_.idAssigned() ) + abs( l1_.idAssigned() );
  return sign * lundSum;
}

// Calculates the cosine of the opening angle theta_ll between the two leptons
//   in the CM frame.
double
DileptonEvent::cosThetaLL()
{
  Hep3Vector lepton0P3Cm = l0_.pCm().vect();
  Hep3Vector lepton1P3Cm = l1_.pCm().vect();
  
  double cosThetaLL;
  double pTotal2 = lepton0P3Cm.mag2() * lepton1P3Cm.mag2();
  if ( pTotal2 <= 0 ) {
    // This case can only happen if a lepton 3-momenta is not real or the 0
    // vector. Thus, I set the cosine to a value that will be cut.
    cosThetaLL = 1.0;
  } else {
    cosThetaLL = lepton0P3Cm.dot( lepton1P3Cm ) / sqrt( pTotal2 );
    if ( cosThetaLL >  1.0 ) cosThetaLL =  1.0;
    if ( cosThetaLL < -1.0 ) cosThetaLL = -1.0;
  }
  return cosThetaLL;
}

// Returns the scalar sum of the lepton CM 3-momenta.
double
DileptonEvent::pSum()
{
  return l0_.pCm().vect().mag() + l1_.pCm().vect().mag();
}

// Returns the positive scalar difference of the lepton CM 3-momenta.
double
DileptonEvent::pDifference()
{
  return abs( l0_.pCm().vect().mag() - l1_.pCm().vect().mag() );
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
