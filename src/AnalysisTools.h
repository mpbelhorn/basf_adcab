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
IpDrDz {
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

//______________________________________________________________________________
// LeptonCandidate class definition and prototypes.

class 
LeptonCandidate {
 public:

  // Constructors and destructor.
  LeptonCandidate();
  LeptonCandidate( Particle lepton, Hep3Vector cmBoostVector );
  ~LeptonCandidate() {};

  // Mutators.
  void setLepton( Particle lepton );
  void setCmBoostVector( Hep3Vector cmBoostVector );

  // Accessors.
  Particle lepton();
  Hep3Vector cmBoostVector();
  Particle particle();

  // Methods.
  double idAssigned();
  int massHypothesis();
  double idTruth();
  double idMother();
  HepLorentzVector pCm();
  HepLorentzVector p();
  double likelihoodE();
  double likelihoodMu();
  double klmHitsChi2();
  int    klmHitsN();
  double klmHitsChi2PerN();
  double svdHitsR();
  double svdHitsZ();

 private:

  // Attributes.
  Particle lepton_;
  Hep3Vector cmBoostVector_;

};

//______________________________________________________________________________
// DileptonEvent class definition and prototypes.

class
DileptonEvent {
 public:

  // Constructors and destructor.
  DileptonEvent();
  DileptonEvent( Particle lepton0, Particle lepton1, Hep3Vector cmBoostVector );
  ~DileptonEvent() {}

  // Mutators.
  void setL0( Particle lepton0 );
  void setL1( Particle lepton1 );
  void setCmBoostVector( Hep3Vector cmBoostVector );
  
  // Accessors.
  LeptonCandidate l0();
  LeptonCandidate l1();
  Hep3Vector cmBoostVector();
  
  // Methods
  double eventType();
  double cosThetaLL();
  double pSum();
  double pDifference();
  // int goodEvent();

 private:
  // Attributes.
  LeptonCandidate l0_;
  LeptonCandidate l1_;
  Hep3Vector cmBoostVector_;

};

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
