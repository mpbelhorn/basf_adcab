//______________________________________________________________________________
// Filename: IpParameters.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing the interaction point parameters dr and dz.
//______________________________________________________________________________

#ifndef IPPARAMETERS_H
#define IPPARAMETERS_H

#include "particle/Particle.h"       // The BELLE Particle Class.
#include MDST_H                      // Panther.

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif


//______________________________________________________________________________
// IPdrdz class definition and prototypes.

// Class for Impact Parameters "dr" and "dz"
class
IpParameters {
 public:

  // Constructors and destructor.
  IpParameters();
  IpParameters(const Mdst_charged& chg, HepPoint3D ip, int massHyp);
  ~IpParameters() {}

  // Mutators.
  void init(const Mdst_charged& chg, HepPoint3D ip, int massHyp);

  // Accessors
  double dr();
  double dz();
  int svdHitsR();
  int svdHitsZ();
  
 private:
  // Attributes.
  double dr_;
  double dz_;
  int svd_hits_r_;
  int svd_hits_z_;
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
