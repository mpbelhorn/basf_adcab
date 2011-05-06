//______________________________________________________________________________
// Filename: Adcab.h
// Version: 2011.02.09.A
// Author: M.P. Belhorn
// Original Date: 2011.02.09
// Description: Analysis of the (A)nomalous (D)ilepton (C)harge (A)symmetry in
//   (B)s0 decays.
//______________________________________________________________________________

#include <cmath>                      // Uses cmath functions.
  
#include "belle.h"                    // BELLE Library.
#include <utility>                    // For std::pair.
#include "event/BelleEvent.h"         // For managing BELLE events.
#include "tuple/BelleTupleManager.h"  // For managing BELLE nTuples.
#include "benergy/BeamEnergy.h"       // For determining run beam energy.
#include "basf/module.h"              // For interfacing with BASF.
#include "basf/module_descr.h"        // ???
#include "particle/Particle.h"        // The BELLE Particle Class.
#include "particle/utility.h"         // Particle class add-on fucntions.
#include "kid/atc_pid.h"              // For particle species separation.
#include "mdst/Muid_mdst.h"           // For muon identification.
#include "eid/eid.h"                  // For electron identification.
#include "mdst/mdst.h"                // For MDST files.
#include "ip/IpProfile.h"             // Beam Interaction Point (IP) analysis
                                      //    tools. Position unit = cm.
  
#include <panther/panther.h>    // Panther - general header?
#include BELLETDF_H             // Panther - runhead/belle_event tables.
#include HEPEVT_H               // Panther - not sure?
#include MDST_H                 // Panther - mdst tables.
#include EVTCLS_H               // Panther - event classification tables.
  
#include "HEPconstants.h"       // PDG masses and constants.
#include "AnalysisTools.h"      // General analysis functions and utilities.
#include "AdcabCuts.h"          // Analysis specfic selection cut constants.
#include "geninfo.h"            // Custom analysis functions to use Zupanc's MC.
#include "userinfo.h"           // Custom analysis functions to use Zupanc's MC.
  
#if defined(BELLE_NAMESPACE)    // Namespace container for backwards
namespace Belle {               //  compatibility with older versions of
#endif                          //  BELLE Library (used for b200611xx onward).
                                //  Must be in all files.

//______________________________________________________________________________
// Adcab Class Definition

// Declare analysis class, inheriting Module from BASF.
class
Adcab : public Module
{
 public:
  // Class functions.
  Adcab();                             // Constructor.
  ~Adcab() {}                          // Destructor.

  // BASF Module functions.
  // These functions are called by BASF via the BASF interface script.
  //                                        Each function is run...
  void init( int * );                       // Once by BASF "initialize".
  void disp_stat( const char* ) {}          // Not used.
  void hist_def();                          // Once by BASF "histogram define".
  void begin_run( BelleEvent*, int* );      // At beginning of each run.
  void event( BelleEvent*, int* );          // For each event in run.
  void end_run( BelleEvent*, int* );        // At end of each run.
  void other( int*, BelleEvent*, int* ) {}  // Not used.
  void term();                              // Once by BASF "terminate".
  
  // BASF passable parameters.
  int basf_parameter_allow_charge_bias;
  int basf_parameter_verbose_log;

  // Particle type (Ptype) constants.
  Ptype ptypeElectron;
  Ptype ptypeEMinus;
  Ptype ptypeEPlus;
  Ptype ptypeNuE;
  Ptype ptypeNuEBar;
  Ptype ptypeMuMinus;
  Ptype ptypeMuPlus;
  Ptype ptypeNuMu;
  Ptype ptypeNuMuBar;
  Ptype ptypeKMinus;
  Ptype ptypeKPlus;
  Ptype ptypeDsMinus;
  Ptype ptypeDsPlus;
  Ptype ptypeDsStarMinus;
  Ptype ptypeDsStarPlus;

  // Runhead / run analysis information.
  int experimentNumber;
  int runNumber;
  int eventNumber;
  int numberOfEvents;
  int numDileptonEvents;

  // Interaction point information.
  HepPoint3D ip;
  HepSymMatrix ipErr;
  int ipUsable;

  // Beam Information.
  HepLorentzVector cmFourMomentum;
  Hep3Vector cmBoostVector;
  double beamEnergyCMFrame;
  double beamEnergyError;
  double lerBeamEnergy;
  double herBeamEnergy;
  double kekbBeamEnergy;     // Uncalibrated energy reported by KEKB.
  double kekbLerBeamEnergy;  // Uncalibrated energy reported by KEKB.
  double kekbHerBeamEnergy;  // Uncalibrated energy reported by KEKB.
  double beamCrossingAngle;

  // Bs - Bsbar Diagnostic Tests.
  std::pair<int,int> num_bs_after_lepton_level;
  std::pair<int,int> num_bs_after_pair_removal;
  std::pair<int,int> num_bs_at_after_event_selection;

  // Flags           
  bool flag_mc;  // Data type flag. 'true' = MC, 'false' = Real Data.

 private:
  BelleTuple *nTuple_;  // Pointer for writing to the n-tuple. 
};

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
