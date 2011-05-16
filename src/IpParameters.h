//______________________________________________________________________________
// Filename: IpParameters.h
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: A class for managing the interaction point parameters dr and dz.
//______________________________________________________________________________

#ifndef IPPARAMETERS_H
#define IPPARAMETERS_H

#include <cmath>                     // Uses cmath functions.

#include "belle.h"                   // BELLE Library.
#include "particle/Particle.h"       // The BELLE Particle Class.
#include "mdst/mdst.h"               // For MDST files.
#include "ip/IpProfile.h"            // Beam Interaction Point (IP) analysis
                                     //   tools. Position unit = cm.

#include <panther/panther.h>         // Panther.
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
  IpParameters( const Mdst_charged&, HepPoint3D, int );
  ~IpParameters() {}

  // Mutators.
  void init( const Mdst_charged&, HepPoint3D, int );

  // Accessors
  double dr();
  double dz();
  
 private:
  // Attributes.
  double dr_;
  double dz_;
  
};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
