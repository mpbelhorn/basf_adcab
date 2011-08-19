//______________________________________________________________________________
// Filename: Adcab.cc
// Version: 2011.02.09.A
// Author: M.P. Belhorn
// Original Date: 2010-07-19
// Description: Analysis of the (A)nomalous (D)ilepton (C)harge (A)symmetry in
//   (B)s0 decays.
//______________________________________________________________________________

#include "Adcab.h"
#include <iostream>

using std::cout;
using std::endl;

#if defined(BELLE_NAMESPACE)    // Namespace container for backwards
namespace Belle {               //  compatibility with older versions of
#endif                          //  BELLE Library (used for b200611xx onward).
                                //  Must be in all files.

//______________________________________________________________________________
// BASF Module interface.

// Registers analysis module in BASF.
extern "C" Module_descr
*mdcl_Adcab() 
{ 
  // Creates pointer to allocated Adcab class object.
  Adcab *module = new Adcab;
  
  // Creates pointer "dscr" to description of Adcab module.
  Module_descr *dscr = new Module_descr("Adcab", module);
  
  // Set up module parameters.
  dscr->define_param ("JPsi_Veto_OS_Only",
      "Apply J/Psi veto to opposite-sign pairs only",
      &module->basf_parameter_allow_charge_bias_);
  dscr->define_param ("Verbose_Log",
      "Writes diagnotic information to the log",
      &module->basf_parameter_verbose_log_);

  // Provide path to pass paramaters to BeamEnergy class.
  BeamEnergy::define_global(dscr);
  return dscr;
}

//______________________________________________________________________________
// Module class member definitions.

// Adcab constructor definition. Loaded at the line
//   "path add_module analysis Adcab" in a BASF script. Place BASF passable
//   parameter initializations here, as they will apply for the entire 
//   analysis run.
Adcab::Adcab()
{ 
  // Define particle types (Ptype) constants for particle class objects.
  // Valid names are those in the qq98 decay file located at
  // $BELLE_TOP_DIR/share/data-files/qq98/decay.dec
  particle_e_minus_ = (Ptype("E-"));
  particle_e_plus_ = (Ptype("E+"));
  particle_mu_minus_ = (Ptype("MU-"));
  particle_mu_plus_ = (Ptype("MU+"));

  // Initialize BASF parameters.
  basf_parameter_allow_charge_bias_ = 0;
  basf_parameter_verbose_log_ = 0;

  cout << "\n\n"
       << "____________________________________________________________\n"
       << " Adcab Analysis Module loaded successfully \n\n"
       << " BASF Parameter Flags Settings:\n"
       << "   basf_parameter_allow_charge_bias_ = " 
       << basf_parameter_allow_charge_bias_ << "\n"
       << "   basf_parameter_verbose_log_ = " 
       << basf_parameter_verbose_log_ << "\n"
       << endl;
  
  return;
}

// Function run at loading of BASF module.
void
Adcab::init(int *)
{
  // Body is intentionally blank.
}

// Function run at termination of BASF module.
void
Adcab::term()
{
  cout << "\n\n"
       << "____________________________________________________________\n"
       << " Adcab Analysis Module terminated successfully \n\n"
       << endl;

  return;
}


//______________________________________________________________________________
// Event Analysis Functions

// This function is exectuted once per data run. Thus it resets the class flags
//   and event run counters, it initializes the runhead information (data type,
//   experiment number, run number, mdst access), it retrieves the run-dependent
//   beam and interaction point information. A basic run summery is dumped to
//   the log for the record.
void
Adcab::begin_run(BelleEvent* evptr, int *status) 
{
  (void)evptr;
  (void)status;
  
  // Set run information to default values.
  experiment_number_ = 0;
  run_number_ = 0;

  // Set default flags.
  flag_mc_ = false;
  
  // Set interaction point and error to default values.
  interaction_point_ = HepPoint3D(0, 0, 0);
  interaction_point_error_ = HepSymMatrix(3, 0);
  flag_good_interaction_point_ = 0;

  // Get BELLE_EVENT runhead manager. This is the data stored in the 
  //   belletdf.tdf panther table.
  Belle_runhead_Manager &runhead_manager = Belle_runhead_Manager::get_manager();
  Belle_runhead &runhead = runhead_manager((Panther_ID) 1);
  experiment_number_ = runhead.ExpNo();
  run_number_ = runhead.RunNo();
  
  // Initialise BeamEnergy class.
  BeamEnergy::begin_run();
  
  // Set beam energy related class variables.
  beam_energy_cm_frame_ = BeamEnergy::E_beam_corr();  
  beam_energy_error_ = BeamEnergy::E_beam_err();
  ler_beam_energy_  = BeamEnergy::E_LER();
  her_beam_energy_ = BeamEnergy::E_HER();
  kekb_beam_energy_ = BeamEnergy::E_beam_orig();
  kekb_ler_beam_energy_ = BeamEnergy::E_LER_orig();
  kekb_her_beam_energy_ = BeamEnergy::E_HER_orig();
  beam_crossing_angle_ = BeamEnergy::Cross_angle();
  cm_boost_     = -BeamEnergy::CMBoost();

  // Initialize PID functions.
  eid::init_data();
  
  // Get interaction point profile data from $BELLE_POSTGRES_SERVER. 
  IpProfile::begin_run();
  
  // Set interaction point and error to run values.
  if (IpProfile::usable()) {
    interaction_point_ = IpProfile::e_position();
    interaction_point_error_ = IpProfile::e_position_err_b_life_smeared();
    flag_good_interaction_point_ = 1;
  } else {
    interaction_point_ = HepPoint3D(0, 0, 0);
    interaction_point_error_ = HepSymMatrix(3, 0);
  }
  
  // Print run information to the log.
  cout << "\n\n"
       << "____________________________________________________________\n"
       << " New Run: "
       << "Experiment " << experiment_number_ 
       << ", Run " << run_number_ 
       << "\n" 
       << endl;

  if (runhead.ExpMC() == 1) {
    flag_mc_ = false; // Set Data type flag to Real Data.
    cout << " Data is Real." << endl;
  } else {
    flag_mc_ = true;  // Set Data type flag to Monte Carlo.
    cout << " Data is Monte Carlo." << endl;
  }
  
  cout << " Actual Beam Energy: " << beam_energy_cm_frame_ 
       << " +/- " << beam_energy_error_ << " GeV\n"
       << " Reported Beam Energy: " << kekb_beam_energy_ << " GeV\n"
       << " BE Class cm_boost_: " << cm_boost_ << "\n"
       << "____________________________________________________________\n"
       << endl;

  return;
}

// This function is run once at the end of a data run. It does nothing but
//   write a message to somewhere. This message does not appear in log
//   although MPB believes it should.
void
Adcab::end_run(BelleEvent* evptr, int *status)
{ 
  (void) evptr;
  (void) status;
  
  cout << "\n\n"
       << "***************************************************\n"
       << "* WHERE IS THIS FUNCTION end_run() EXECUTED?!?!?! *\n"
       << "***************************************************\n"
       << endl;

  return;
}

// This is function contains the core of the analysis module and is run
//   for each event in a data run. Each event is analyized to find 
//   signal dilepton event candidates and records information
//   about those events to the n-tuple.
void
Adcab::event(BelleEvent* evptr, int* status)
{
  (void)evptr;
  (void)status;

  // Get the event number. 
  // The bitwise operation "& ~(0 << 28)" forces the event number to count from 
  //   zero again after EvtNo() reaches 2^28. Not sure why this is necessary?
  Belle_event_Manager& EvMgr = Belle_event_Manager::get_manager();
  event_number_ = EvMgr[0].EvtNo() & ~(~0 << 28);

  if (basf_parameter_verbose_log_) {
    cout << "____________________________________________________________\n"
         << "New Event #" << event_number_ << "(MC: " << flag_mc_ << ")"
         << endl;
  }

  // Check the event classification information for HadronB criteria.
  Evtcls_hadronic_flag_Manager &hadronFlagManager
      = Evtcls_hadronic_flag_Manager::get_manager();
  Evtcls_hadronic_flag &hadronFlags(hadronFlagManager(Panther_ID(1)));
  float hadronb_code = hadronFlags.hadronic_flag(2);

  // Print hadron flag to the log.
  if (basf_parameter_verbose_log_ && hadronb_code < 10) {
    cout << "  BAD EVENT :: Fails hadronB criteria." << endl;
  }

  // Get the Fox-wolfram R2 value for the event. Spherical events accepted (Low
  //   R2 values). If R2 is not calculated in the hadron info table, event will
  //   be accepted by default with an R2 of -1 to point out in the ntuple that
  //   the cut is "effectively off".
  double fox_wolfram_r2 = -1;
  Evtcls_hadron_info_Manager& 
      hadron_info_manager = Evtcls_hadron_info_Manager::get_manager();
  if (hadron_info_manager.count()) {
    fox_wolfram_r2 = hadron_info_manager[0].R2();
  }
  
  // Instantiate constant data structures stored in external header files.
  // Might as well make them static since
  //     they will be needed again and they never change.
  static PDGmasses masses;
  static AdcabCuts cuts;
    
  // Define lists (vector template) to store event particles.
  // Need a list for all mother and daughter particle species.
  // Note that these vectors are static and will persist until BASF is closed.
  //     They therefore must(!) be cleared for each call of Adcab::event().
  static std::vector<LeptonCandidate> lepton_candidates(5);
  static std::vector<DileptonEvent> dilepton_event_candidates(10);
  lepton_candidates.clear();
  dilepton_event_candidates.clear();
  
  // Alias the MDST charged manager, which contains the measured charged tracks
  //     for each event and make pointers to it's boundaries.
  Mdst_charged_Manager
      &mdst_charged_particles = Mdst_charged_Manager::get_manager();
  MdstChargedIterator first_mdst_charged = mdst_charged_particles.begin();
  MdstChargedIterator last_mdst_charged = mdst_charged_particles.end();
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "  Passed event initialization." << endl;
  }

  // Write unique event-level data to the n-tuple. This information is the
  //     first of many rows associated with a physical event.
  // Column names can be no greater than eight (8) characters long.
  nTuple_events_->column("exp_no"  , experiment_number_);
  nTuple_events_->column("run_no"  , run_number_);
  nTuple_events_->column("evt_no"  , event_number_);
  nTuple_events_->column("fw_r2"   , fox_wolfram_r2);
  nTuple_events_->column("hadronb" , hadronb_code);
  nTuple_events_->dumpData();
  
  // Populate the lepton candidate lists.
  for (MdstChargedIterator lepton_scan_iterator = first_mdst_charged;
      lepton_scan_iterator != last_mdst_charged; ++lepton_scan_iterator) {

    if (basf_parameter_verbose_log_) {
      cout << "  >>> New Lepton Candidate <<<" << endl;
    }
    // Alias the current particle as "charged_particle".
    const Mdst_charged &charged_particle = *lepton_scan_iterator;

    // Get electron and muon likelihoods.
    eid charged_particle_eid(charged_particle);
    double eid_probability = charged_particle_eid.prob(3, -1, 5);
    Muid_mdst charged_particle_muid(charged_particle);
    double muid_probability = charged_particle_muid.Muon_likelihood();

    // Reject particle if below both electron and muon likelihood cuts.
    bool good_eid = eid_probability >= cuts.minEidProb;
    bool good_muid = ((muid_probability >= cuts.minMuidProb) &&
        (charged_particle_muid.Chi_2() != 0));
    if (!(good_eid || good_muid)) continue;
    if (basf_parameter_verbose_log_) {
      cout << "    Passed PID check." << endl;
    }

    // Cut on IP dr and dz and SVD hits.
    // This is to make sure that particles were created near the IP.
    // TODO - 2010.08.11 - Is mass hypothesis = 3 (kaon) appropriate? Check with
    //                       authorities!
    IpParameters ip_parameters(charged_particle, interaction_point_, 3);
    if (abs(ip_parameters.dr()) > cuts.maxIpDr) continue;
    if (abs(ip_parameters.dz()) > cuts.maxIpDz) continue;
    if (basf_parameter_verbose_log_) {
      cout << "    Passed dr/dz check." << endl;
    }
    if (ip_parameters.svdHitsR() < cuts.minSvdRHits) continue;
    if (ip_parameters.svdHitsZ() < cuts.minSvdZHits) continue;
    if (basf_parameter_verbose_log_) {
      cout << "    Passed SVD hits check." << endl;
    }

    // Reject if the particle track has polar angle pointing outside the
    //   barrel (p is given closest to coord. origin - see mdst table).
    Hep3Vector momentum3_lab(
        charged_particle.px(), charged_particle.py(), charged_particle.pz());
    double track_cosine_polar_angle = momentum3_lab.cosTheta();
    if (track_cosine_polar_angle < cuts.minLeptonCosTheta) continue;
    if (track_cosine_polar_angle > cuts.maxLeptonCosTheta) continue;
    if (basf_parameter_verbose_log_) {
      cout << "    Passed Barrel PID check." << endl;
    }

    // Create an electron candidate...
    double electric_charge = charged_particle.charge();
    Particle electron_candidate(charged_particle,
        electric_charge > 0 ? particle_e_plus_ : particle_e_minus_);
    // ... and a muon candidate.
    Particle muon_candidate(charged_particle,
        electric_charge > 0 ? particle_mu_plus_ : particle_mu_minus_);

    // Check if CM momentum is good assuming muon and electron masses.
    // First, get the lab-frame 4-momentum assuming each mass
    //     and boost it to CM frame.
    HepLorentzVector electron_candidate_p4_cm(electron_candidate.p());
    HepLorentzVector muon_candidate_p4_cm(muon_candidate.p());
    electron_candidate_p4_cm.boost(cm_boost_);
    muon_candidate_p4_cm.boost(cm_boost_);
    
    double electron_candidate_p4_cm_mag = electron_candidate_p4_cm.vect().mag();
    double muon_candidate_p4_cm_mag = muon_candidate_p4_cm.vect().mag();
    bool good_electron_momentum = !(
        electron_candidate_p4_cm_mag < cuts.minLeptonMomentumCm ||
        electron_candidate_p4_cm_mag > cuts.maxLeptonMomentumCm);
    bool good_muon_momentum = !(
        muon_candidate_p4_cm_mag < cuts.minLeptonMomentumCm ||
        electron_candidate_p4_cm_mag > cuts.maxLeptonMomentumCm);
    bool good_muon(good_muid && good_muon_momentum ? true : false);
    bool good_electron(good_eid && good_electron_momentum ? true : false);

    if (!(good_muon || good_electron)) continue;
    if (basf_parameter_verbose_log_) {
      cout << "    Passed CM-frame momentum check." << endl;
    }
    
    // Remove possible pair production and J/Psi daughters.
    for (MdstChargedIterator jpsi_pair_iterator = first_mdst_charged;
        jpsi_pair_iterator != last_mdst_charged; ++jpsi_pair_iterator) {
      const Mdst_charged &jpsi_pair_particle = *jpsi_pair_iterator;

      // Reject case where pointers point to same object.
      if (lepton_scan_iterator == jpsi_pair_iterator) continue;

      // If allowing a charge bias and pair is SS, skip to the next pair.
      if (basf_parameter_allow_charge_bias_ &&
          (electric_charge == jpsi_pair_particle.charge())) {
        continue;
      }

      if (good_muon) {
        // Check to see if the candidate makes a bad muon. First, assume
        //     the pair particle is a muon.
        Particle other_muon(jpsi_pair_particle,
            jpsi_pair_particle.charge() > 0 ?
            particle_mu_plus_ : particle_mu_minus_);

        // Calculate the invariant mass of the two particles.
        double pair_invariant_mass = abs(
            (muon_candidate.p() + other_muon.p()).m());

        // If at any time candidate proves to be likely from a J/Psi, the flag
        //   is switched.
        double delta_mass = pair_invariant_mass - cuts.massJPsi;
        if (cuts.minMuMuJPsiCandidate < delta_mass &&
            delta_mass < cuts.maxMuMuJPsiCandidate) {
          if (basf_parameter_verbose_log_) {
            cout << "    BAD MU CANDIDATE :: Likely J/Psi daughter." << endl;
          }
          good_muon = false;
        }
      }

      if (good_electron) {
        // Check to see if the candidate makes a bad electron. First, assume
        //     the pair particle is an electron.
        Particle other_electron(jpsi_pair_particle,
            jpsi_pair_particle.charge() > 0 ?
            particle_e_plus_ : particle_e_minus_);

        // Calculate the invariant mass of the electron
        //   candidate and the other charged track.
        double pair_invariant_mass = abs(
            (electron_candidate.p() + other_electron.p()).m());
        double delta_mass = pair_invariant_mass - cuts.massJPsi;
  
        // Cut possible pair production electrons or J/Psi daughters.
        if (pair_invariant_mass < cuts.minEPlusEMinusMass) {
          if (basf_parameter_verbose_log_) {
            cout << "    BAD E CANDIDATE :: Likely pair production." << endl;
          }
          good_electron = false;
        } else if (cuts.minElElJPsiCandidate < delta_mass &&
            delta_mass < cuts.maxElElJPsiCandidate) {
          if (basf_parameter_verbose_log_) {
            cout << "    BAD E CANDIDATE :: Likely J/Psi daughter." << endl;
          }
          good_electron = false;
        }
      }
    }

    // While the candidate is still good, add it to the lepton list. Only add
    //     the candidate to the lepton list once! If it is a good muon
    //     candidate, then consider it a muon, otherwise it is deemed an
    //     electron, but flag the verbose log if the particle passes as both.
    Particle *good_candidate = NULL;
    if (good_muon) {
      setMCtruth(muon_candidate);
      good_candidate = &muon_candidate;
    } else if (good_electron) {
      setMCtruth(electron_candidate);
      good_candidate = &electron_candidate;
    }

    if (good_electron && good_muon && basf_parameter_verbose_log_) {
      cout << "    WARNING :: Good muon and good electron candidate." << endl;
    }

    // TODO - Dump lepton candidate information to the ntuple.
    if (good_candidate) {
      LeptonCandidate lepton((*good_candidate), cm_boost_);
      lepton_candidates.push_back(lepton);
      
      nTuple_leptons_->column("charge"   , electric_charge);
      nTuple_leptons_->column("mass"     , lepton.p().mag());
      nTuple_leptons_->column("good_mu"  , good_muon);
      nTuple_leptons_->column("good_el"  , good_electron);
      nTuple_leptons_->column("id_asn"   , lepton.idAssigned());
      nTuple_leptons_->column("id_tru"   , lepton.idTrue());
      nTuple_leptons_->column("id_mom"   , lepton.idMom());
      nTuple_leptons_->column("eid_prob" , eid_probability);
      nTuple_leptons_->column("muid_prb" , muid_probability);
      nTuple_leptons_->column("muid_rto" , lepton.klmChi2PerHits());
      nTuple_leptons_->column("p_lb_mag" , lepton.p().vect().mag());
      nTuple_leptons_->column("p_cm_mag" , lepton.pCm().vect().mag());
      nTuple_leptons_->column("e_cm"     , lepton.pCm().e());
      nTuple_leptons_->column("p_cm_x"   , lepton.pCm().px());
      nTuple_leptons_->column("p_cm_y"   , lepton.pCm().py());
      nTuple_leptons_->column("p_cm_z"   , lepton.pCm().pz());
      nTuple_leptons_->column("cos_pol"  , lepton.p().cosTheta());
      nTuple_leptons_->column("ip_dr"    , ip_parameters.dr());
      nTuple_leptons_->column("ip_dz"    , ip_parameters.dz());
      nTuple_leptons_->column("svd_hitr" , ip_parameters.svdHitsR());
      nTuple_leptons_->column("svd_hitz" , ip_parameters.svdHitsZ());
      nTuple_leptons_->dumpData();
    }
  }

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "  Passed lepton candidate selection." << endl;
    cout << "    Total lepton candidates: " << lepton_candidates.size() << endl;
    cout << "  Searching for dilepton candidates." << endl;
  }

  // TODO - Find the number of kaons in the event.
  // TODO - Find the number of phi mesons in the event.
 
  // Find dilepton event candidates.
  // Loop over the lepton list.
  for (LeptonCandidateIterator j = lepton_candidates.begin();
      j != lepton_candidates.end(); ++j) {
    // Loop over remaining leptons in the list to check all lepton pairs.
    for (LeptonCandidateIterator i = j; i != lepton_candidates.end(); ++i) {
      // Exclude the case where both iterators point to the same particle.
      if (i == j) continue;
      if (basf_parameter_verbose_log_) {
        cout << "  >>> New Event Candidate <<<" << endl;
      }

      // Determine higher momentum lepton and add it to an event candidate.
      // Reference the higher momentum lepton as "lepton0"
      //  and the lower momentum lepton as "lepton1".
      LeptonCandidate *lepton0 = &(*i);
      LeptonCandidate *lepton1 = &(*j);
      if ((*j).pCm().mag() > (*i).pCm().mag()) {
        lepton0 = &(*j);
        lepton1 = &(*i);
      }
      DileptonEvent event_candidate(*lepton0, *lepton1);
      
      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosine_cm_opening_angle = event_candidate.cosThetaLL();
      if (cosine_cm_opening_angle < cuts.minCosThetaLLCm ||
          cosine_cm_opening_angle > cuts.maxCosThetaLLCm) {
        continue;
      }
      
      if (basf_parameter_verbose_log_) {
        cout << "    Passed lepton opening angle check." << endl;
        cout << "    Commiting event candidate to n-tuple." << endl;
      }

      // Add event candidate information to the n-tuple.
      nTuple_dileptons_->column("evt_type", event_candidate.eventType());
      nTuple_dileptons_->column("evt_sign", event_candidate.eventSign());
      nTuple_dileptons_->column("llcostha", event_candidate.cosThetaLL());
      nTuple_dileptons_->dumpData();
    }
  }
  
  return;
}

// Specifies n-tuples to write and names for root histograms.
void
Adcab::hist_def()
{
  extern BelleTupleManager *BASF_Histogram;   // Define a BASF Histogram
  
  BelleTupleManager *tm = BASF_Histogram;
  const char *event_variables = "exp_no "
                                "run_no "
                                "evt_no "
                                "fw_r2 "
                                "hadronb";

  const char *lepton_variables = "charge "
                                 "mass "
                                 "id_asn "
                                 "id_tru "
                                 "id_mom "
                                 "eid_p "
                                 "muid_p "
                                 "muid_r "
                                 "p_cm_mag "
                                 "p_lb_mag "
                                 "e_cm "
                                 "p_cm_x "
                                 "p_cm_y "
                                 "p_cm_z "
                                 "costheta "
                                 "ip_dr "
                                 "ip_dz "
                                 "svdr_hit "
                                 "svdz_hit";
                                 
                                
  const char *kaon_variables = "charge "
                               "mass "
                               "id_asn "
                               "id_tru "
                               "id_mom "
                               "eid_p "
                               "muid_p "
                               "muid_r "
                               "p_cm_mag "
                               "p_lb_mag "
                               "e_cm "
                               "p_cm_x "
                               "p_cm_y "
                               "p_cm_z "
                               "costheta "
                               "ip_dr "
                               "ip_dz "
                               "svdr_hit "
                               "svdz_hit";
                                 
  const char *dilepton_variables = "evt_type "
                                   "evt_sign "
                                   "llcostha";

  nTuple_events_ = tm->ntuple("Events", event_variables, 1);
  nTuple_leptons_ = tm->ntuple("Leptons", lepton_variables, 2);
  nTuple_kaons_ = tm->ntuple("Kaons", kaon_variables, 3);
  nTuple_dileptons_ = tm->ntuple("Dilepton", dilepton_variables, 4);
  
  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
