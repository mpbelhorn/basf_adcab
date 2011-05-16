//______________________________________________________________________________
// Filename: AnalysisTools.h
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Declaration of custom analysis classes and functions. 
//______________________________________________________________________________

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
#include "mdst/Muid_mdst.h"           // For muon identification. 
#include "ip/IpProfile.h"            // Beam Interaction Point (IP) analysis
                                     //   tools. Position unit = cm.

#include <panther/panther.h>         // Panther.
#include BELLETDF_H                  // Panther.
#include HEPEVT_H                    // Panther.
#include MDST_H                      // Panther.

#include "HEPconstants.h"      // PDG masses and constants.
#include "AdcabCuts.h"    // Analysis specfic selection cut constants.

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

//______________________________________________________________________________
// LeptonCandidate class definition and prototypes.

class 
LeptonCandidate {
 public:

  // Constructors and destructor.
  LeptonCandidate();
  LeptonCandidate( Particle lepton, Hep3Vector cm_boost );
  ~LeptonCandidate() {};

  // Mutators.
  void set_lepton( Particle lepton );
  void set_cm_boost_vector( Hep3Vector cm_boost );

  // Accessors.
  Particle lepton();
  Hep3Vector cm_boost();

  // Methods.
  double idAssigned();
  double idTrue();
  double idMom();
  int massHypothesis();
  HepLorentzVector pCm();
  HepLorentzVector p();
  double electronProbability();
  double muonProbability();
  double klmHitsChi2();
  int klmHits();
  double klmChi2PerHits();
  double svdRadialHits();
  double svdAxialHits();

 private:

  // Attributes.
  Particle lepton_;
  Hep3Vector cm_boost_;

};

//______________________________________________________________________________
// DileptonEvent class definition and prototypes.

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
