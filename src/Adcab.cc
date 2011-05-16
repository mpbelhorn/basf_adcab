//______________________________________________________________________________
// Filename: Adcab.cc
// Version: 2011.02.09.A
// Author: M.P. Belhorn
// Original Date: 2010-07-19
// Description: Analysis of the (A)nomalous (D)ilepton (C)harge (A)symmetry in
//   (B)s0 decays.
//______________________________________________________________________________

#include "Adcab.h"              // Adcab Analysis header.
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
  (void)evptr;
  (void)status;
  
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
    cout << " Bad event: Fails hadronB criteria." << endl;
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
  PDGmasses masses;
  AdcabCuts cuts;
    
  // Define lists (vector template) to store event particles.
  // Need a list for all mother and daughter particle species.
  // Reserve a large amount of memory (50 elements) to avoid 
  //   relocating list to a new block of memory in case of large
  //   events.
  static std::vector<Particle> initial_electron_candidates(50);
  static std::vector<Particle> electron_candidates(50);
  static std::vector<Particle> initial_muon_candidates(50);
  static std::vector<Particle> muon_candidates(50);
  static std::vector<Particle> lepton_candidates(100);

  // Define list to store candidate dilepton events.
  static std::vector<DileptonEvent> dilepton_event_candidates(50);
  
  // Ensure lists are empty so as not to consume too much memory.
  initial_electron_candidates.clear();
  electron_candidates.clear();
  initial_muon_candidates.clear();
  muon_candidates.clear();
  lepton_candidates.clear();
  dilepton_event_candidates.clear();
  
  // Alias the MDST charged manager, which contains
  // the measured charged tracks for each event.
  Mdst_charged_Manager
      &mdst_charged_particles = Mdst_charged_Manager::get_manager();
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed initialization." << endl;
  }
  
  // Populate the lepton candidate lists.
  for (MdstChargedConstIterator i = mdst_charged_particles.begin();
      i != mdst_charged_particles.end(); ++i) {

    // Alias the current particle as "charged_particle".
    const Mdst_charged &charged_particle = *i;
    
    // Get electron and muon liklihoods.
    eid charged_particle_eid(charged_particle);
    double eid_probability = charged_particle_eid.prob(3, -1, 5);
    Muid_mdst charged_particle_muid(charged_particle);
    double muid_probability = charged_particle_muid.Muon_likelihood();

    // Reject particle if below both electron and muon liklihood cuts.
    if (eid_probability < cuts.minEidProb &&
        muid_probability < cuts.minMuidProb) {
      continue;
    }
    
    // Cut on IP dr and dz.
    // This is to make sure that particles were created near the IP.
    // TODO - 2010.08.11 - Is mass hypothesis = 3 (kaon) appropriate? Check with
    //                       authorities!
    IpParameters ip_parameters(charged_particle, interaction_point_, 3);
    if (abs(ip_parameters.dr()) > cuts.maxIpDr ||
        abs(ip_parameters.dz()) > cuts.maxIpDz) {
      continue;
    }
      
    // Mdst_trk member function mhyp(int hypID) returns the fitted
    //   track parameters assuming certain particle mass hypotheses set
    //   by hypID. Possible mass hypotheses are:
    //   hypID = 0:e; 1:mu; 2:pi; 3:K; 4:p.

    // Mdst_trk_fit member function nhits(int detID) returns the
    //   number of associated hits in CDC or SVD detector elements
    //   as per the value of detID. Possible detID values are:
    //   detID = 0:axial-wire; 1:stereo-wire; 2:cathode; 
    //   3:SVD-rphi; 4:SVD-z.
    // Note that nhits(SVD)=0 indicates that the track fit is performed using
    //   only CDC info.
    // See mdst.tdf for information about what information is contained in the
    //   mdst files.
    
    // Reject particles with too few hits in the SVD.
    // For the purpose of checking the number of hits in the SVD,
    //   it is not necessary to use a specific mass hypothesis.
    Mdst_trk_fit &charged_particle_track = charged_particle.trk().mhyp(1);
    if (charged_particle_track.nhits(3) < cuts.minSvdRHits) continue;
    if (charged_particle_track.nhits(4) < cuts.minSvdZHits) continue;

    // Reject particle if likelihoods are the same from eid and muid.
    if (eid_probability == muid_probability) continue;

    // Assume the particle species is that of the highest PID liklihood.
    if (eid_probability > muid_probability) {
      // Cut on minimum EID.
      if (eid_probability < cuts.minEidProb) continue;
      
      // Assuming particle is an electron.
      Particle electron_candidate(charged_particle,
          charged_particle.charge() > 0 ? particle_e_plus_ : particle_e_minus_);

      // Reject if the particle track has polar angle pointing outside the
      //   barrel (p is given closest to coord. origin - see mdst table).
      if (electron_candidate.p().cosTheta() < cuts.minLeptonCosTheta ||
          electron_candidate.p().cosTheta() > cuts.maxLeptonCosTheta) {
        continue;
      }

      // Cut on lepton CM momentum.
      // Get the lab-frame momentum and boost it to CM frame.
      HepLorentzVector candidate_p4_cm(electron_candidate.p());
      candidate_p4_cm.boost(cm_boost_);
      double candidate_p4_cm_mag = candidate_p4_cm.vect().mag();
      if (candidate_p4_cm_mag < cuts.minLeptonMomentumCm ||
          candidate_p4_cm_mag > cuts.maxLeptonMomentumCm) {
        continue;
      }

      // Store the MC truth to the lepton candidate.
      setMCtruth(electron_candidate);

      // Add electron_candidate to list of e+/- candidates.
      initial_electron_candidates.push_back(electron_candidate);
    } else {
      // Cut on minimum MUID.
      if (muid_probability < cuts.minMuidProb) continue;
      
      // Assuming particle is a muon.
      Particle muon_candidate(charged_particle,
          charged_particle.charge() > 0 ?
          particle_mu_plus_ : particle_mu_minus_);

      // Reject if the particle track has polar angle pointing outside the
      //   barrel (p is given closest to coord. origin - see mdst table).
      if (muon_candidate.p().cosTheta() < cuts.minLeptonCosTheta ||
          muon_candidate.p().cosTheta() > cuts.maxLeptonCosTheta) {
        continue;
      }

      // Cut on lepton CM momentum.
      // Get the lab-frame momentum and boost it to CM frame.
      HepLorentzVector candidate_p4_cm(muon_candidate.p());
      candidate_p4_cm.boost(cm_boost_);
      double candidate_p4_cm_mag = candidate_p4_cm.vect().mag();
      if (candidate_p4_cm_mag < cuts.minLeptonMomentumCm ||
          candidate_p4_cm_mag > cuts.maxLeptonMomentumCm) {
        continue;
      }

      // Store the MC truth to the lepton candidate.
      setMCtruth(muon_candidate);
      
      // Add muon_candidate to list of mu+/- candidates.
      muon_candidates.push_back(muon_candidate);
    }
  } // End for() loop populating lepton lists.

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed lepton candidate selection..." << endl;
    cout << "  electron candidates: " << initial_electron_candidates.size()
        << "." << endl;
    cout << "      muon candidates: " << initial_muon_candidates.size()
        << "." << endl;
  }

  // Remove possible pair production electrons and J/Psi candidates.
  for (ParticleIterator j = initial_electron_candidates.begin();
      j != initial_electron_candidates.end(); ++j) {
    const Particle &electron_candidate = *j;
    
    // Flag for electrons that are not candidates for pair production daughters.
    // By default, assume all electrons are good.
    bool flag_good_electron = true;

    // If the invariant mass of an electron candidate any other charged
    //   track is smaller than the cut value (nominally 100 MeV), the
    //   electron candidate is rejected.
    
    for (MdstChargedConstIterator i = mdst_charged_particles.begin();
        i != mdst_charged_particles.end(); ++i) {
      const Mdst_charged &mdst_particle = *i;

      // Reject case where pointers point to same object.
      if (electron_candidate.relation().mdstCharged() == mdst_particle) {
        continue;
      }

      Particle other_particle(mdst_particle,
          mdst_particle.charge() > 0 ? particle_e_plus_ : particle_e_minus_);
      
      // We need to know if the pair is same sign or opposite sign.
      bool ss_pair = (electron_candidate.charge() == other_particle.charge());
      bool allow_charge_bias = basf_parameter_allow_charge_bias_ && ss_pair;
      
      // If allowing a charge bias and pair is SS, skip to the next pair.
      if (allow_charge_bias) {
        continue;
      }

      // Calculate the invariant mass of the electron candidate and the other
      //   charged track.
      HepLorentzVector candidate_p4 = electron_candidate.p();
      HepLorentzVector other_particle_p4 = other_particle.p();
      double pair_invariant_mass = abs((candidate_p4 + other_particle_p4).m());
      double delta_mass = pair_invariant_mass - cuts.massJPsi;

      // Print diagnostic information to the log.
      if (basf_parameter_verbose_log_ == 2) {
        cout << " e-chg mass: " << pair_invariant_mass << endl;
        cout << " delta mass: " << delta_mass << muon_candidates.size() << endl;
      }

      // Cut possible pair production electrons or J/Psi daughters.
      if (pair_invariant_mass < cuts.minEPlusEMinusMass) {
        flag_good_electron = false;
      } else if (cuts.minElElJPsiCandidate < delta_mass &&
          delta_mass < cuts.maxElElJPsiCandidate) {
        flag_good_electron = false;
      }
    }
    
    // While electron_candidate is still a good candidate, add it to the list.
    if (flag_good_electron) {
      electron_candidates.push_back(electron_candidate);
    }
  }
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed electron pair production and J/Psi vetoes." << endl;
    cout << "   Remaining electrons: " << electron_candidates.size() << "."
         << endl;
  }
  
  // Remove possible J/Psi daughter muons.
  for (ParticleIterator j = initial_muon_candidates.begin();
      j != initial_muon_candidates.end(); ++j) {
    const Particle &muon_candidate = *j;
    
    // Flag for muons that are not candidates J/psi daughters.
    // By default, assume all muons are good.
    bool flagGoodMuon = true;

    // If the difference of the invariant mass of a muon candidate and every
    //   other opposite charged tracks is within the cut range, the candidate
    //   is rejected.
    for (MdstChargedConstIterator i = mdst_charged_particles.begin();
        i != mdst_charged_particles.end(); ++i) {
      const Mdst_charged &mdst_particle = *i;

      // Reject case where pointers point to same object.
      if (muon_candidate.relation().mdstCharged() == mdst_particle) continue;

      Particle other_particle(mdst_particle,
          mdst_particle.charge() > 0 ? particle_mu_plus_ : particle_mu_minus_);
      
      // We need to know if the pair is same sign or opposite sign.
      bool ss_pair = (muon_candidate.charge() == other_particle.charge());
      bool allow_charge_bias = basf_parameter_allow_charge_bias_ && ss_pair;
      
      // If allowing a charge bias and pair is SS, skip to the next pair.
      if (allow_charge_bias) {
        continue;
      }
      
      // Calculate the invariant mass of the two particles.
      HepLorentzVector candidate_p4 = muon_candidate.p();
      HepLorentzVector other_particle_p4 = other_particle.p();
      double pair_invariant_mass = abs((candidate_p4 + other_particle_p4).m());

      // If at any time candidate proves to be likely from a J/Psi, the flag
      //   is switched.
      double delta_mass = pair_invariant_mass - cuts.massJPsi;
      if (cuts.minMuMuJPsiCandidate < delta_mass &&
          delta_mass < cuts.maxMuMuJPsiCandidate) {
        flagGoodMuon = false;
      }
    }

    // While candidate is still good, add it to the muon list.
    if (flagGoodMuon) {
      muon_candidates.push_back(muon_candidate);
    }
  }
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed muon J/Psi veto." << endl;
    cout << "   Remaining muons: " << muon_candidates.size() << "." << endl;
  }
  
  // Combine all leptons into a single list.
  lepton_candidates.reserve(
      electron_candidates.size() + muon_candidates.size());
  lepton_candidates.insert(lepton_candidates.end(), electron_candidates.begin(),
      electron_candidates.end());
  lepton_candidates.insert(lepton_candidates.end(), muon_candidates.begin(),
      muon_candidates.end());
 
  // Find dilepton event candidates.
  // Loop over the lepton list.
  for (ParticleIterator j = lepton_candidates.begin();
      j != lepton_candidates.end(); ++j) {
    
    // Loop over remaining leptons in the list to check all lepton pairs.
    for (ParticleIterator i = j;
        i != lepton_candidates.end(); ++i) {
      
      // Exclude the case where both iterators point to the same particle.
      if (i == j) continue;
      
      // Alias the first lepton as "lepton0"
      //  and the second electron as "lepton1".
      Particle &lepton0 = *j;
      Particle &lepton1 = *i;
      
      // Cut on each lepton's CM momentum.
      // TODO - 2010.08.10 - Perhaps this should be done when populating the
      //                     lepton lists.
      // Get the lab-frame momentum for each lepton and boost it to CM frame.
      HepLorentzVector lepton0_p4_cm(lepton0.p());
      HepLorentzVector lepton1_p4_cm(lepton1.p());
      lepton0_p4_cm.boost(cm_boost_);
      lepton1_p4_cm.boost(cm_boost_);

      // Make the CM momentum cut.
      double lepton0_p3_cm_mag = lepton0_p4_cm.vect().mag();
      double lepton1_p3_cm_mag = lepton1_p4_cm.vect().mag();
      if (lepton0_p3_cm_mag < cuts.minLeptonMomentumCm ||
          lepton0_p3_cm_mag > cuts.maxLeptonMomentumCm) {
        continue;
      }

      if (lepton1_p3_cm_mag < cuts.minLeptonMomentumCm ||
          lepton1_p3_cm_mag > cuts.maxLeptonMomentumCm) {
        continue;
      }
      
      // Determine higher momentum lepton and add it to an event candidate.
      //   lepton0 is henceforth considered the higher momentum lepton.
      Particle &greater_p_lepton = (lepton0_p3_cm_mag > lepton1_p3_cm_mag ?
          lepton0 : lepton1);
      Particle &lower_p_lepton = (lepton0_p3_cm_mag > lepton1_p3_cm_mag ?
          lepton1 : lepton0);
      DileptonEvent event_candidate(greater_p_lepton, lower_p_lepton,
          cm_boost_);
      
      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosine_cm_opening_angle = event_candidate.cosThetaLL();
      if (cosine_cm_opening_angle < cuts.minCosThetaLLCm ||
          cosine_cm_opening_angle > cuts.maxCosThetaLLCm) {
        continue;
      }

      // Add event candidate to the list of dilepton events.
      dilepton_event_candidates.push_back(event_candidate);
    }
  }
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed event candidate selection." << endl;
    cout << "   " << dilepton_event_candidates.size() << " candidates found."
         << endl;
  }

  // Write the event candidate information to n-tuple.
  for (DileptonEventIterator i = dilepton_event_candidates.begin();
      i != dilepton_event_candidates.end(); ++i) {
    DileptonEvent &event_candidate = *i;

    // Rebuild the dr and dz track parameters for each lepton.
    IpParameters l0_ip_parameters(
        event_candidate.l0().lepton().relation().mdstCharged(),
        interaction_point_, event_candidate.l0().massHypothesis());
    IpParameters l1_ip_parameters(
        event_candidate.l1().lepton().relation().mdstCharged(),
        interaction_point_, event_candidate.l1().massHypothesis());

    // Write data to the n-tuple. As of 2010.08.25, each row in the n-tuple will
    //   consist of a single dilepton event candidate.
    // Column names can be no greater than eight (8) characters long.
    // Format:      "name    ", value)
    // Event singular information
    nTuple_->column("exp_no"  , experiment_number_);
    nTuple_->column("run_no"  , run_number_);
    nTuple_->column("evt_no"  , event_number_);
    nTuple_->column("fw_r2"   , fox_wolfram_r2);
    nTuple_->column("hadronb" , hadronb_code);
    
    // Event candidate information.
    nTuple_->column("evnttype", event_candidate.eventType());
    nTuple_->column("llcostha", event_candidate.cosThetaLL());
    nTuple_->column("pcm_sum" , event_candidate.pSum());
    nTuple_->column("pcm_dif" , event_candidate.pDifference());

    // Lepton information.
    nTuple_->column("l0_chg"  , event_candidate.l0().lepton().charge());
    nTuple_->column("l1_chg"  , event_candidate.l1().lepton().charge());
    nTuple_->column("l0_m"    , event_candidate.l0().p().mag());
    nTuple_->column("l1_m"    , event_candidate.l1().p().mag());
    nTuple_->column("l0_idasn", event_candidate.l0().idAssigned());
    nTuple_->column("l1_idasn", event_candidate.l1().idAssigned());
    nTuple_->column("l0_idtru", event_candidate.l0().idTrue());
    nTuple_->column("l1_idtru", event_candidate.l1().idTrue());
    nTuple_->column("l0_idmom", event_candidate.l0().idMom());
    nTuple_->column("l1_idmom", event_candidate.l1().idMom());
    nTuple_->column("l0_eidp" , event_candidate.l0().electronProbability());
    nTuple_->column("l1_eidp" , event_candidate.l1().electronProbability());
    nTuple_->column("l0_muidp", event_candidate.l0().muonProbability());
    nTuple_->column("l1_muidp", event_candidate.l1().muonProbability());
    nTuple_->column("l0_muidr", event_candidate.l0().klmChi2PerHits());
    nTuple_->column("l1_muidr", event_candidate.l1().klmChi2PerHits());
    nTuple_->column("l0_pcm"  , event_candidate.l0().pCm().vect().mag());
    nTuple_->column("l1_pcm"  , event_candidate.l1().pCm().vect().mag());
    nTuple_->column("l0_plab" , event_candidate.l0().p().vect().mag());
    nTuple_->column("l1_plab" , event_candidate.l1().p().vect().mag());
    nTuple_->column("l0_cme"  , event_candidate.l0().pCm().e());
    nTuple_->column("l1_cme"  , event_candidate.l1().pCm().e());
    nTuple_->column("l0_cmpx" , event_candidate.l0().pCm().px());
    nTuple_->column("l1_cmpx" , event_candidate.l1().pCm().px());
    nTuple_->column("l0_cmpy" , event_candidate.l0().pCm().py());
    nTuple_->column("l1_cmpy" , event_candidate.l1().pCm().py());
    nTuple_->column("l0_cmpz" , event_candidate.l0().pCm().pz());
    nTuple_->column("l1_cmpz" , event_candidate.l1().pCm().pz());
    nTuple_->column("l0_costh", event_candidate.l0().p().cosTheta());
    nTuple_->column("l1_costh", event_candidate.l1().p().cosTheta());
    nTuple_->column("l0_dr"   , l0_ip_parameters.dr());
    nTuple_->column("l1_dr"   , l1_ip_parameters.dr());
    nTuple_->column("l0_dz"   , l0_ip_parameters.dz());
    nTuple_->column("l1_dz"   , l1_ip_parameters.dz());
    nTuple_->column("l0_svdr" , event_candidate.l0().svdRadialHits());
    nTuple_->column("l1_svdr" , event_candidate.l1().svdRadialHits());
    nTuple_->column("l0_svdz" , event_candidate.l0().svdAxialHits());
    nTuple_->column("l1_svdz" , event_candidate.l1().svdAxialHits());
    nTuple_->dumpData();
  }
  
  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << " Passed commit to ntuple." << endl;
  }
  
  return;
}

// Specifies n-tuples to write and names for root histograms.
void
Adcab::hist_def()
{
  extern BelleTupleManager *BASF_Histogram;   // Define a BASF Histogram
  BelleTupleManager *tm = BASF_Histogram;   
  const char *ntList = "exp_no "
                       "run_no "
                       "evt_no "
                       "fw_r2 "
                       "hadronb "
                       "evnttype "
                       "llcostha "
                       "pcm_sum "
                       "pcm_dif "
                       "l0_chg "
                       "l1_chg "
                       "l0_m "
                       "l1_m "
                       "l0_idasn "
                       "l1_idasn "
                       "l0_idtru "
                       "l1_idtru "
                       "l0_idmom "
                       "l1_idmom "
                       "l0_eidp "
                       "l1_eidp "
                       "l0_muidp "
                       "l1_muidp "
                       "l0_muidr "
                       "l1_muidr "
                       "l0_pcm "
                       "l1_pcm "
                       "l0_plab "
                       "l1_plab "
                       "l0_cme "
                       "l1_cme "
                       "l0_cmpx "
                       "l1_cmpx "
                       "l0_cmpy "
                       "l1_cmpy "
                       "l0_cmpz "
                       "l1_cmpz "
                       "l0_costh "
                       "l1_costh "
                       "l0_dr "
                       "l1_dr "
                       "l0_dz "
                       "l1_dz "
                       "l0_svdr "
                       "l1_svdr "
                       "l0_svdz "
                       "l1_svdz";

  nTuple_ = tm->ntuple("Dilepton", ntList, 1);
  
  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
