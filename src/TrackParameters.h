//______________________________________________________________________________
// Filename: TrackParameters.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing the track parameters.
//______________________________________________________________________________

#ifndef IPPARAMETERS_H
#define IPPARAMETERS_H

#include "particle/Particle.h"       // The BELLE Particle Class.
#include MDST_H                      // Panther.

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif

// Class for altering track parameterization.
class
TrackParameters {
 public:

  // Constructors and destructor.
  TrackParameters();
  TrackParameters(const Particle& particle, HepPoint3D new_pivot);
  TrackParameters(const TrackParameters &that);
  TrackParameters &operator= (const TrackParameters &that);
  virtual ~TrackParameters();

  // Accessors
  double &dr();
  double &dz();
  int &svdHitsR();
  int &svdHitsZ();
  
 private:
  // Attributes.
  double *dr_;
  double *dz_;
  int *svd_hits_r_;
  int *svd_hits_z_;
  int *mass_hypothesis_;
  Helix *helix_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
