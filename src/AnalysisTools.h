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
  double id_assigned();
  int mass_hypothesis();
  double id_true();
  double id_mother();
  HepLorentzVector p_cm();
  HepLorentzVector p();
  double electron_probability();
  double muon_probability();
  double klm_hits_chi2();
  int number_of_klm_hits();
  double klm_chi2_per_hits();
  double svd_radial_hits();
  double svd_axial_hits();

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
