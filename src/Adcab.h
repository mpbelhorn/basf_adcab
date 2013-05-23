#ifndef ADCAB_H
#define ADCAB_H
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
#include "TrackParameters.h"    // Custom class managing the track parameters.
#include "DileptonEvent.h"      // Class for managing dilepton event info.
#include "EntryTypes.h"         // Structure for managing ntuple line data.
#include "AdcabCuts.h"          // Analysis specfic selection cut constants.
#include "geninfo.h"            // Custom analysis functions to use Zupanc's MC.
#include "userinfo.h"           // For saving special info to Particle instances.
#include "VertexFit.h"          // Custom analysis functions to for vertexing.

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
  int basf_parameter_is_continuum_; // 1 = continuum data, 0 = on resonance.
  int basf_parameter_mc_stream_number_;
  int basf_parameter_scale_momentum_;


  // Particle type (Ptype) constants.
  Ptype particle_e_minus_;
  Ptype particle_e_plus_;
  Ptype particle_mu_minus_;
  Ptype particle_mu_plus_;
  Ptype particle_k_plus_;
  Ptype particle_k_minus_;

  // Runhead / run analysis information.
  int experiment_number_;
  int run_number_;
  int event_number_;
  int stream_number_;

  // Interaction point information.
  HepPoint3D interaction_point_;
  HepSymMatrix interaction_point_error_;
  int flag_good_interaction_point_;

  typedef std::vector<Mdst_charged>::const_iterator MdstChargedIterator;
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
  // Pointers for writing to the n-tuples.
  BelleTuple* nTuple_charged_;
  BelleTuple* nTuple_phi_;
  BelleTuple* nTuple_dileptons_;

  // Pointers to the BASF level histograms.
  // Muon ID cuts.
  BelleHistogram* generated_muon_multiplicity_;
  BelleHistogram* muon_muid_histogram_;
  BelleHistogram* muon_muid_multiplicity_;
  BelleHistogram* muon_svd_r_histogram_;
  BelleHistogram* muon_svd_r_multiplicity_;
  BelleHistogram* muon_svd_z_histogram_;
  BelleHistogram* muon_svd_z_multiplicity_;
  BelleHistogram* muon_klm_signature_histogram_;
  BelleHistogram* muon_klm_signature_multiplicity_;

  // Electron ID cuts.
  BelleHistogram* generated_electron_multiplicity_;
  BelleHistogram* electron_eid_histogram_;
  BelleHistogram* electron_eid_multiplicity_;
  BelleHistogram* electron_svd_r_histogram_;
  BelleHistogram* electron_svd_r_multiplicity_;
  BelleHistogram* electron_svd_z_histogram_;
  BelleHistogram* electron_svd_z_multiplicity_;

  // Kaon ID cuts.
  BelleHistogram* generated_kaon_multiplicity_;
  BelleHistogram* kaon_kid_histogram_;
  BelleHistogram* kaon_kid_multiplicity_;

  // Prompt Muon cuts.
  std::vector<BelleHistogram*> muon_ip_dr_histograms_; // 0: BG, 1: Bs, 2: Bd, 3: Secondary
  BelleHistogram* muon_ip_dr_multiplicities_;
  std::vector<BelleHistogram*> muon_ip_dz_histograms_;
  BelleHistogram* muon_ip_dz_multiplicities_;
  std::vector<BelleHistogram*> muon_polar_angle_cosine_histograms_;
  BelleHistogram* muon_polar_angle_cosine_multiplicities_;
  std::vector<BelleHistogram*> muon_cm_momentum_histograms_;
  BelleHistogram* muon_cm_momentum_multiplicities_;
  BelleHistogram* muon_jpsi_multiplicities_;
  BelleHistogram* muon_accepted_multiplicities_;

  // Prompt Electron cuts.
  std::vector<BelleHistogram*> electron_ip_dr_histograms_; // 0: BG, 1: Bs, 2: Bd, 3: Secondary
  BelleHistogram* electron_ip_dr_multiplicities_;
  std::vector<BelleHistogram*> electron_ip_dz_histograms_;
  BelleHistogram* electron_ip_dz_multiplicities_;
  std::vector<BelleHistogram*> electron_polar_angle_cosine_histograms_;
  BelleHistogram* electron_polar_angle_cosine_multiplicities_;
  std::vector<BelleHistogram*> electron_cm_momentum_histograms_;
  BelleHistogram* electron_cm_momentum_multiplicities_;
  BelleHistogram* electron_pair_multiplicities_;
  BelleHistogram* electron_jpsi_multiplicities_;
  BelleHistogram* electron_accepted_multiplicities_;

  // Phi cuts.
  BelleHistogram* generated_phi_multiplicity_;
  BelleHistogram* generated_phi_to_dikaon_multiplicity_;

};

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
#endif
