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
#include "IpParameters.h"       // Custom class managing the IP.
#include "LeptonCandidate.h"    // Class for managing lepton candidate info.
#include "DileptonEvent.h"      // Class for managing dilepton event info.
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
  //                                      Each function is run...
  void init(int *);                       // Once by BASF "initialize".
  void disp_stat(const char*) {}          // Not used.
  void hist_def();                        // Once by BASF "histogram define".
  void begin_run(BelleEvent*, int*);      // At beginning of each run.
  void event(BelleEvent*, int*);          // For each event in run.
  void end_run(BelleEvent*, int*);        // At end of each run.
  void other(int*, BelleEvent*, int*) {}  // Not used.
  void term();                            // Once by BASF "terminate".
  
  // BASF passable parameters.
  int basf_parameter_allow_charge_bias_;
  int basf_parameter_verbose_log_;

  // Particle type (Ptype) constants.
  Ptype particle_e_minus_;
  Ptype particle_e_plus_;
  Ptype particle_mu_minus_;
  Ptype particle_mu_plus_;

  // Runhead / run analysis information.
  int experiment_number_;
  int run_number_;
  int event_number_;

  // Interaction point information.
  HepPoint3D interaction_point_;
  HepSymMatrix interaction_point_error_;
  int flag_good_interaction_point_;
  
  typedef std::vector<Mdst_charged>::const_iterator MdstChargedConstIterator;
  typedef std::vector<Particle>::iterator ParticleIterator;
  typedef std::vector<DileptonEvent>::iterator DileptonEventIterator;

  // Beam Information.
  Hep3Vector cm_boost_;
  double beam_energy_cm_frame_;
  double beam_energy_error_;
  double ler_beam_energy_;
  double her_beam_energy_;
  double kekb_beam_energy_;      // Uncalibrated energy reported by KEKB.
  double kekb_ler_beam_energy_;  // Uncalibrated energy reported by KEKB.
  double kekb_her_beam_energy_;  // Uncalibrated energy reported by KEKB.
  double beam_crossing_angle_;

  // Flags           
  bool flag_mc_;  // Data type flag. 'true' = MC, 'false' = Real Data.

 private:
  BelleTuple *nTuple_;  // Pointer for writing to the n-tuple. 
};

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
