//
//******************************************************************************
// Filename: AnalysisTools.h
// Version: 2010.09.20.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Declaration of custom analysis classes and functions. 
//******************************************************************************

#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

#include <cmath>                     // Uses cmath functions.

#include "belle.h"                   // BELLE Library.
#include "event/BelleEvent.h"        // For managing BELLE events.
#include "tuple/BelleTupleManager.h" // For managing BELLE nTuples.
#include "benergy/BeamEnergy.h"      // For determining run beam energy.
#include "basf/module.h"             // For interfacing with BASF.
#include "basf/module_descr.h"       // ???
#include "particle/Particle.h"       // The BELLE Particle Class.
#include "kid/atc_pid.h"             // For particle species separation.
#include "eid/eid.h"                 // For electron identification.
#include "mdst/mdst.h"               // For MDST files.
#include "ip/IpProfile.h"            // Beam Interaction Point (IP) analysis
                                     //   tools. Position unit = cm.

#include <panther/panther.h>         // Panther.
#include BELLETDF_H                  // Panther.
#include HEPEVT_H                    // Panther.
#include MDST_H                      // Panther.

#include "HEPconstants.h"      // PDG masses and constants.
#include "BsdileptonCuts.h"    // Analysis specfic selection cut constants.

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif


//******************************************************************************
// IPdrdz class definition and prototypes.
//******************************************************************************

// Class for Impact Parameters "dr" and "dz"
class IpDrDz {
 public:

  // Constructors and destructor.
  IpDrDz();
  IpDrDz( const Mdst_charged&, HepPoint3D, int );
  ~IpDrDz() {}

  // Mutators.
  void setDrDz( const Mdst_charged&, HepPoint3D, int );

  // Accessors
  double dr();
  double dz();
  
 private:
  // Attributes.
  double dr_;
  double dz_;
  
};


//******************************************************************************
// DileptonEvent class definition and prototypes.
//******************************************************************************

class DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent( Particle lepton0, Particle lepton1 );
  ~DileptonEvent() {}

  // Mutators.
  void setLepton0( Particle lepton0 );
  void setLepton1( Particle lepton1 );
  
  // Accessors.
  Particle lepton0();
  Particle lepton1();
  
  // Methods
  //
  // Calculates the cosine of the opening angle theta_ll between the two leptons
  //   in the CM frame.
  double cosThetaLLCm( Hep3Vector cmBoostVector );
  double eventType();
  double l0Id();
  double l1Id();
  double l0MotherId();
  double l1MotherId();
  //
  // Checks if the mothers of both leptons are MC signal Bs's (Pythia ID 531).
  // Returns 1 if true, 0 if not true, and -1 if not both MC info available. 
  int mcBsParents();

private:
  // Attributes.
  Particle lepton0_;
  Particle lepton1_;

};


//******************************************************************************
// General function prototypes.
//******************************************************************************


#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
