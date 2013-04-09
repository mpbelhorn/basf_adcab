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
  dscr->define_param ("MC_Stream_Number",
      "Allows the MC stream number to be used to avoid name collisions",
      &module->basf_parameter_mc_stream_number_);
  dscr->define_param ("Is_Continuum",
      "Flags data as being continuum",
      &module->basf_parameter_is_continuum_);
  dscr->define_param ("Scale_Momentum",
      "Flag the analysis to scale momentum to the Y(5S) energy",
      &module->basf_parameter_scale_momentum_);

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
  particle_k_minus_ = (Ptype("K-"));
  particle_k_plus_ = (Ptype("K+"));

  // Initialize BASF parameters.
  basf_parameter_allow_charge_bias_ = 0;
  basf_parameter_verbose_log_ = 0;
  basf_parameter_mc_stream_number_ = 0;
  basf_parameter_is_continuum_ = 0;
  basf_parameter_scale_momentum_ = 0;

  return;
}

// Function run at loading of BASF module.
void
Adcab::init(int *)
{
  cout << "\n\n"
       << "____________________________________________________________\n"
       << " Adcab Analysis Module loaded successfully \n\n"
       << " BASF Parameter Flags Settings:\n"
       << "   basf_parameter_allow_charge_bias_ = "
       << basf_parameter_allow_charge_bias_ << "\n"
       << "   basf_parameter_verbose_log_ = "
       << basf_parameter_verbose_log_ << "\n"
       << "   basf_parameter_mc_stream_number_ = "
       << basf_parameter_mc_stream_number_ << "\n"
       << "   basf_parameter_is_continuum_ = "
       << basf_parameter_is_continuum_ << "\n"
       << "   basf_parameter_scale_momentum_ = "
       << basf_parameter_scale_momentum_ << "\n"
       << endl;
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
    cout << " Data is Monte Carlo - Stream " << basf_parameter_mc_stream_number_
         << endl;
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
    cout << "    BAD EVENT :: Fails hadronB criteria." << endl;
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

  // Instantiate cuts as stored in external header file.
  // Might as well make them static since
  //     they will be needed again and they never change.
  static AdcabCuts cuts;

  // Instantiate PID classes. These are constant and
  //     should only be initialized once.
  static atc_pid pid_kaon_to_pi(3, 1, 5, 3, 2);
  static atc_pid pid_kaon_to_pr(3, 1, 5, 3, 4);

  // Define lists (vector template) to store event particles.
  // Need a list for all mother and daughter particle species.
  // Note that these vectors are static and will persist until BASF is closed.
  //     They therefore must(!) be cleared for each call of Adcab::event().
  if (basf_parameter_verbose_log_ > 1) {
    cout << "    Initializing particle containers." << endl;
  }
  static std::vector<ParticleCandidate> lepton_candidates(5);
  static std::vector<ParticleCandidate> kaon_candidates(10);
  static std::vector<ParticleCandidate> true_kaons(10);
  static std::vector<ParticleCandidate> pi_candidates(10);
  static std::vector<Particle> phi_candidates(5);
  static std::vector<Particle> d_candidates(5);
  static std::vector<DileptonEvent> dilepton_event_candidates(10);
  lepton_candidates.clear();
  kaon_candidates.clear();
  true_kaons.clear();
  pi_candidates.clear();
  phi_candidates.clear();
  d_candidates.clear();
  dilepton_event_candidates.clear();

  // Alias the MDST charged manager, which contains the measured charged tracks
  //     for each event and make pointers to it's boundaries.
  Mdst_charged_Manager
      &mdst_charged_particles = Mdst_charged_Manager::get_manager();
  MdstChargedIterator first_mdst_charged = mdst_charged_particles.begin();
  MdstChargedIterator last_mdst_charged = mdst_charged_particles.end();

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "    Passed event initialization." << endl;
  }

  // Populate the lepton candidate lists.
  for (MdstChargedIterator lepton_scan_iterator = first_mdst_charged;
      lepton_scan_iterator != last_mdst_charged; ++lepton_scan_iterator) {

    if (basf_parameter_verbose_log_ > 1) {
      cout << "    New charged Track" << endl;
    }
    // Alias the current particle as "charged_particle".
    const Mdst_charged &charged_particle = *lepton_scan_iterator;

    // Get electron and muon likelihoods.
    eid charged_particle_eid(charged_particle);
    Muid_mdst charged_particle_muid(charged_particle);

    double muid_probability = charged_particle_muid.Muon_likelihood();
    double eid_probability = charged_particle_eid.prob(3, -1, 5);
    double kaon_to_pion_likelihood = pid_kaon_to_pi.prob(charged_particle);
    double kaon_to_proton_likelihood = pid_kaon_to_pr.prob(charged_particle);

    // Reject particle if below both electron and muon likelihood cuts.
    bool good_muon = ((muid_probability >= cuts.minMuidProb) &&
        (charged_particle_muid.Chi_2() != 0));
    bool good_electron = eid_probability >= cuts.minEidProb;
    bool good_kaon = (
        (kaon_to_pion_likelihood > cuts.minKaonToPionLikelihood) &&
        (kaon_to_proton_likelihood > cuts.minKaonToProtonLikelihood));

    if (basf_parameter_verbose_log_ > 1) {
      cout << "        ["<< good_muon << good_electron << good_kaon << "] "
           << "PID check." << endl;
    }

    // Create particle class objects of each species.
    double electric_charge = charged_particle.charge();
    Particle muon_particle(charged_particle,
        electric_charge > 0 ? particle_mu_plus_ : particle_mu_minus_);
    Particle electron_particle(charged_particle,
        electric_charge > 0 ? particle_e_plus_ : particle_e_minus_);
    Particle kaon_particle(charged_particle,
        electric_charge > 0 ? particle_k_plus_ : particle_k_minus_);

    ParticleCandidate muon_candidate(muon_particle, cm_boost_,
        interaction_point_, basf_parameter_scale_momentum_);
    ParticleCandidate electron_candidate(electron_particle, cm_boost_,
        interaction_point_, basf_parameter_scale_momentum_);

    // Cut on IP dr and dz and SVD hits.
    // This is to make sure that particles were created near the IP.
    if ((abs(muon_candidate.track().dr()) > cuts.maxIpDr) ||
        (abs(muon_candidate.track().dz()) > cuts.maxIpDz) ||
        (muon_candidate.svdRHits() < cuts.minSvdRHits) ||
        (muon_candidate.svdZHits() < cuts.minSvdZHits)) {
      good_muon = false;
    }
    if ((abs(electron_candidate.track().dr()) > cuts.maxIpDr) ||
        (abs(electron_candidate.track().dz()) > cuts.maxIpDz) ||
        (electron_candidate.svdRHits() < cuts.minSvdRHits) ||
        (electron_candidate.svdZHits() < cuts.minSvdZHits)) {
      good_electron = false;
    }
    if (basf_parameter_verbose_log_ > 1) {
      cout << "        ["<< good_muon << good_electron << good_kaon << "] "
           << "dr/dz and SVD check." << endl;
    }
    // Reject if the particle track has polar angle pointing outside the
    //   barrel (p is given closest to coord. origin - see mdst table).
    //   The lab momentum is the same for each particle candidate, so we use the
    //   muon's particle instance for the job.
    double track_cosine_polar_angle = muon_particle.p().cosTheta();
    if ((track_cosine_polar_angle < cuts.minLeptonCosTheta) ||
        (track_cosine_polar_angle > cuts.maxLeptonCosTheta)) {
      good_muon = false;
      good_electron = false;
    }
    if (basf_parameter_verbose_log_ > 1) {
      cout << "        ["<< good_muon << good_electron << good_kaon << "] "
           << "Barrel intersection check." << endl;
    }

    // Check if CM momentum is good assuming muon and electron masses.
    if ((muon_candidate.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (muon_candidate.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      good_muon = false;
    }
    if ((electron_candidate.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (electron_candidate.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      good_electron = false;
    }
    if (basf_parameter_verbose_log_ > 1) {
      cout << "        ["<< good_muon << good_electron << good_kaon << "] "
           << "CM-frame momentum check." << endl;
    }

    // Tag possible pair production and J/Psi daughters.
    bool vetoed_muon = false;
    bool vetoed_electron = false;
    for (MdstChargedIterator jpsi_pair_iterator = first_mdst_charged;
        (good_muon || good_electron) && !(vetoed_muon && vetoed_electron) &&
        (jpsi_pair_iterator != last_mdst_charged); ++jpsi_pair_iterator) {
      const Mdst_charged &jpsi_pair_particle = *jpsi_pair_iterator;

      // Reject case where pointers point to same object.
      if (lepton_scan_iterator == jpsi_pair_iterator) continue;

      // If allowing a charge bias and pair is SS, skip to the next pair.
      if (basf_parameter_allow_charge_bias_ &&
          (electric_charge == jpsi_pair_particle.charge())) {
        continue;
      }

      if (good_muon && !vetoed_muon) {
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
          if (basf_parameter_verbose_log_ > 1) {
            cout << "            Bad muon: Likely J/Psi daughter." << endl;
          }
          vetoed_muon = true;
        }
      }

      if (good_electron && !vetoed_electron) {
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
          if (basf_parameter_verbose_log_ > 1) {
            cout << "            Bad electron: Likely pair production." << endl;
          }
          vetoed_electron = true;
        } else if (cuts.minElElJPsiCandidate < delta_mass &&
            delta_mass < cuts.maxElElJPsiCandidate) {
          if (basf_parameter_verbose_log_ > 1) {
            cout << "            Bad electron: Likely J/Psi daughter." << endl;
          }
          vetoed_electron = true;
        }
      }
    }
    if (basf_parameter_verbose_log_ > 1) {
      cout << "        ["<< good_muon << good_electron << good_kaon << "] "
           << "J/Psi and pair production veto." << endl;
    }

    // While the candidate is still good, add it to the lepton list. Only add
    //     the candidate to the lepton list once! If it is a good muon
    //     candidate, then consider it a muon, otherwise it is deemed an
    //     electron, but flag the verbose log if the particle passes as both.
    ParticleCandidate *good_lepton = NULL;
    bool vetoed_lepton = false;
    if (good_muon) {
      if (vetoed_muon) vetoed_lepton = true;
      setMCtruth(muon_candidate.particle());
      good_lepton = &muon_candidate;
    } else if (good_electron) {
      if (vetoed_electron) vetoed_lepton = true;
      setMCtruth(electron_candidate.particle());
      good_lepton = &electron_candidate;
    }

    // Dump lepton candidate information to the ntuple for all lepton
    //   candidates including those that fail the J/psi or pair production veto,
    //   but only pass non-vetoed leptons on to the good lepton list.
    if (good_lepton) {
      if (basf_parameter_verbose_log_ > 1) {
        cout << "        Committing lepton to ntuple." << endl;
      }
      ParticleCandidate lepton((*good_lepton));
      if (!vetoed_lepton) {
        lepton_candidates.push_back(lepton);
      }

      nTuple_leptons_->column("stm_no"  , basf_parameter_mc_stream_number_);
      nTuple_leptons_->column("exp_no"  , experiment_number_);
      nTuple_leptons_->column("run_no"  , run_number_);
      nTuple_leptons_->column("evt_no"  , event_number_);
      nTuple_leptons_->column("is_mc"   , flag_mc_);
      nTuple_leptons_->column("is_cntnm", basf_parameter_is_continuum_);
      nTuple_leptons_->column("cm_enrgy", beam_energy_cm_frame_);
      nTuple_leptons_->column("charge"  , electric_charge);
      nTuple_leptons_->column("mass"    , lepton.p().mag());
      nTuple_leptons_->column("good_mu" , good_muon);
      nTuple_leptons_->column("veto_mu" , vetoed_muon);
      nTuple_leptons_->column("good_el" , good_electron);
      nTuple_leptons_->column("veto_el" , vetoed_electron);
      nTuple_leptons_->column("good_k"  , good_kaon);
      nTuple_leptons_->column("pid_k_pi", kaon_to_pion_likelihood);
      nTuple_leptons_->column("pid_k_pr", kaon_to_proton_likelihood);
      nTuple_leptons_->column("id_asn"  , lepton.idAssigned());
      nTuple_leptons_->column("id_tru"  , lepton.idTrue());
      nTuple_leptons_->column("id_mom"  , lepton.idMom());
      nTuple_leptons_->column("eid_prob", eid_probability);
      nTuple_leptons_->column("muid_prb", muid_probability);
      nTuple_leptons_->column("muid_rto", lepton.klmChi2PerHits());
      nTuple_leptons_->column("p_lb_mag", lepton.p().rho());
      nTuple_leptons_->column("p_cm_mag", lepton.pCm().rho());
      nTuple_leptons_->column("e_cm"    , lepton.pCm().e());
      nTuple_leptons_->column("p_cm_x"  , lepton.pCm().px());
      nTuple_leptons_->column("p_cm_y"  , lepton.pCm().py());
      nTuple_leptons_->column("p_cm_z"  , lepton.pCm().pz());
      nTuple_leptons_->column("cos_pol" , lepton.p().cosTheta());
      nTuple_leptons_->column("ip_dr"   , lepton.track().dr());
      nTuple_leptons_->column("ip_dz"   , lepton.track().dz());
      nTuple_leptons_->column("svd_hitr", lepton.svdRHits());
      nTuple_leptons_->column("svd_hitz", lepton.svdZHits());
      nTuple_leptons_->dumpData();
    }

    // Treat the particle as a charged kaon.
    setMCtruth(kaon_particle);
    ParticleCandidate kaon(kaon_particle, cm_boost_, interaction_point_);
    if (good_kaon || abs(kaon.idTrue()) == 321) {
      if (good_kaon && !good_lepton) {
        if (basf_parameter_verbose_log_ > 1) {
          cout << "        Committing kaon to list." << endl;
        }
        kaon_candidates.push_back(kaon);
      }
      if (abs(kaon.idTrue()) == 321) {
        true_kaons.push_back(kaon);
      }

      nTuple_kaons_->column("stm_no"  , basf_parameter_mc_stream_number_);
      nTuple_kaons_->column("exp_no"  , experiment_number_);
      nTuple_kaons_->column("run_no"  , run_number_);
      nTuple_kaons_->column("evt_no"  , event_number_);
      nTuple_kaons_->column("is_mc"   , flag_mc_);
      nTuple_kaons_->column("is_cntnm", basf_parameter_is_continuum_);
      nTuple_kaons_->column("cm_enrgy", beam_energy_cm_frame_);
      nTuple_kaons_->column("charge"  , electric_charge);
      nTuple_kaons_->column("mass"    , kaon.p().mag());
      nTuple_kaons_->column("good_mu" , good_muon);
      nTuple_kaons_->column("veto_mu" , vetoed_muon);
      nTuple_kaons_->column("good_el" , good_electron);
      nTuple_kaons_->column("veto_el" , vetoed_electron);
      nTuple_kaons_->column("good_k"  , good_kaon);
      nTuple_kaons_->column("pid_k_pi", kaon_to_pion_likelihood);
      nTuple_kaons_->column("pid_k_pr", kaon_to_proton_likelihood);
      nTuple_kaons_->column("id_asn"  , kaon.idAssigned());
      nTuple_kaons_->column("id_tru"  , kaon.idTrue());
      nTuple_kaons_->column("id_mom"  , kaon.idMom());
      nTuple_kaons_->column("eid_prob", eid_probability);
      nTuple_kaons_->column("muid_prb", muid_probability);
      nTuple_kaons_->column("muid_rto", kaon.klmChi2PerHits());
      nTuple_kaons_->column("p_lb_mag", kaon.p().rho());
      nTuple_kaons_->column("p_cm_mag", kaon.pCm().rho());
      nTuple_kaons_->column("e_cm"    , kaon.pCm().e());
      nTuple_kaons_->column("p_cm_x"  , kaon.pCm().px());
      nTuple_kaons_->column("p_cm_y"  , kaon.pCm().py());
      nTuple_kaons_->column("p_cm_z"  , kaon.pCm().pz());
      nTuple_kaons_->column("cos_pol" , kaon.p().cosTheta());
      nTuple_kaons_->column("ip_dr"   , kaon.track().dr());
      nTuple_kaons_->column("ip_dz"   , kaon.track().dz());
      nTuple_kaons_->column("svd_hitr", kaon.svdRHits());
      nTuple_kaons_->column("svd_hitz", kaon.svdZHits());
      nTuple_kaons_->dumpData();
    }
  }

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "    Passed track selection." << endl;
    cout << "      Lepton candidates: " << lepton_candidates.size() << endl;
    cout << "      Kaon candidates: " << kaon_candidates.size() << endl;
  }

  // Find Phi canididates.
  for (ParticleCandidateIterator kaon1 = kaon_candidates.begin();
      kaon1 != kaon_candidates.end(); ++kaon1) {
    for (ParticleCandidateIterator kaon2 = kaon1;
        kaon2 != kaon_candidates.end(); ++kaon2) {
      if (kaon1 == kaon2) continue;
      Particle &k1 = kaon1->particle();
      Particle &k2 = kaon2->particle();
      if (k1.charge() == k2.charge()) continue;
      HepLorentzVector phi_momentum = k1.p() + k2.p();
      Particle phi_candidate(phi_momentum, Ptype("PHI"));
      int error = fitPhiVertex(phi_candidate);
      if (!error) {
        phi_candidates.push_back(phi_candidate);
      }
    }
  }
  for (ParticleIterator phi_iter = phi_candidates.begin();
      phi_iter != phi_candidates.end(); ++ phi_iter) {
    setMCtruth(*phi_iter);
    cout << "Phi vertex chisq: "
         << dynamic_cast<UserInfo&>(phi_iter->userInfo()).chisq()
         << endl;
    cout << "Phi truth: "
         << IDhep(*phi_iter)
         << endl;
  }

  // Find good dilepton event candidates.
  // Check all lepton pairs.
  for (ParticleCandidateIterator outer_lepton = lepton_candidates.begin();
      outer_lepton != lepton_candidates.end(); ++outer_lepton) {
    for (ParticleCandidateIterator inner_lepton = outer_lepton;
        inner_lepton != lepton_candidates.end(); ++inner_lepton) {
      // Exclude the case where both iterators point to the same particle.
      if (outer_lepton == inner_lepton) continue;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "    New dilepton event candidate" << endl;
      }

      // Determine higher momentum lepton and add it to an event candidate.
      // Reference the higher momentum lepton as "lepton0"
      //  and the lower momentum lepton as "lepton1". First, assume the outer
      //  loop lepton has the greater momentum, then check that.
      ParticleCandidate *greater_p_lepton;
      ParticleCandidate *lesser_p_lepton;
      if ((*inner_lepton).pCm().mag() > (*outer_lepton).pCm().mag()) {
        greater_p_lepton = &(*inner_lepton);
        lesser_p_lepton = &(*outer_lepton);
      } else {
        greater_p_lepton = &(*outer_lepton);
        lesser_p_lepton = &(*inner_lepton);
      }
      ParticleCandidate &l0 = *greater_p_lepton;
      ParticleCandidate &l1 = *lesser_p_lepton;

      DileptonEvent event_candidate(l0, l1);

      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosine_cm_opening_angle = event_candidate.cosThetaLL();
      if (cosine_cm_opening_angle < cuts.minCosThetaLLCm ||
          cosine_cm_opening_angle > cuts.maxCosThetaLLCm) {
        continue;
      }

      if (basf_parameter_verbose_log_ > 1) {
        cout << "        Found good dilepton event candidate." << endl;
      }
      // Save the good dilepton event candidates the list.
      dilepton_event_candidates.push_back(event_candidate);
    }
  }

  // TODO - Scan through the dilepton event list and see if two kaon candidates
  //   "vertex" well with the lepton candidates. Different kaons with smallest
  //   vertex error for each lepton should be added to the dilepton event and
  //   dumped out with the dilepton event below.

  // Write all dilepton event candidates to the ntuple.
  for (DileptonEventIterator i = dilepton_event_candidates.begin();
      i != dilepton_event_candidates.end(); i++) {
    if (basf_parameter_verbose_log_) {
      cout << "    Commiting Event Candidate to the ntuple" << endl;
    }
    DileptonEvent &event_candidate = *i;
    ParticleCandidate &l0 = event_candidate.l0();
    ParticleCandidate &l1 = event_candidate.l1();

    int number_of_candidate_k_plus = 0;
    int number_of_candidate_k_minus = 0;
    for (ParticleCandidateIterator kaon = kaon_candidates.begin();
        kaon < kaon_candidates.end(); kaon++) {
      if (kaon->particle().charge() > 0) {
        ++number_of_candidate_k_plus;
      } else {
        ++number_of_candidate_k_minus;
      }
    }
    int number_of_true_k_plus = 0;
    int number_of_true_k_minus = 0;
    for (ParticleCandidateIterator kaon = true_kaons.begin();
        kaon < true_kaons.end(); kaon++) {
      if (kaon->particle().charge() > 0) {
        ++number_of_true_k_plus;
      } else {
        ++number_of_true_k_minus;
      }
    }

    // Column names can be no greater than eight (8) characters long.
    // Write event-level data to the n-tuple.
    nTuple_dileptons_->column("stm_no"  , basf_parameter_mc_stream_number_);
    nTuple_dileptons_->column("exp_no"  , experiment_number_);
    nTuple_dileptons_->column("run_no"  , run_number_);
    nTuple_dileptons_->column("evt_no"  , event_number_);
    nTuple_dileptons_->column("is_mc"   , flag_mc_);
    nTuple_dileptons_->column("is_cntnm", basf_parameter_is_continuum_);
    nTuple_dileptons_->column("fw_r2"   , fox_wolfram_r2);
    nTuple_dileptons_->column("hadronb" , hadronb_code);
    nTuple_dileptons_->column("cm_enrgy", beam_energy_cm_frame_);
    nTuple_dileptons_->column("n_can_kp", number_of_candidate_k_plus);
    nTuple_dileptons_->column("n_can_km", number_of_candidate_k_minus);
    nTuple_dileptons_->column("n_tru_kp", number_of_true_k_plus);
    nTuple_dileptons_->column("n_tru_km", number_of_true_k_minus);

    // Write dilepton-level data to the n-tuple.
    nTuple_dileptons_->column("typ_asn" , event_candidate.eventTypeAssigned());
    nTuple_dileptons_->column("typ_tru" , event_candidate.eventTypeTrue());
    nTuple_dileptons_->column("evt_sign", event_candidate.eventSign());
    nTuple_dileptons_->column("cos_thta", event_candidate.cosThetaLL());
    nTuple_dileptons_->column("inv_mass", event_candidate.invariantMass());

    // Write lepton-level data to the n-tuple.
    nTuple_dileptons_->column("l0_chrge", l0.particle().charge());
    nTuple_dileptons_->column("l0_idasn", l0.idAssigned());
    nTuple_dileptons_->column("l0_idtru", l0.idTrue());
    nTuple_dileptons_->column("l0_idmom", l0.idMom());
    nTuple_dileptons_->column("l0_plab" , l0.p().rho());
    nTuple_dileptons_->column("l0_pcm"  , l0.pCm().rho());
    nTuple_dileptons_->column("l0_e_cm" , l0.pCm().e());
    nTuple_dileptons_->column("l0_pcm_x", l0.pCm().px());
    nTuple_dileptons_->column("l0_pcm_y", l0.pCm().py());
    nTuple_dileptons_->column("l0_pcm_z", l0.pCm().pz());
    nTuple_dileptons_->column("l0_ip_dr", l0.track().dr());
    nTuple_dileptons_->column("l0_ip_dz", l0.track().dz());
    nTuple_dileptons_->column("l0_svdr" , l0.svdRHits());
    nTuple_dileptons_->column("l0_svdz" , l0.svdZHits());

    nTuple_dileptons_->column("l1_chrge", l1.particle().charge());
    nTuple_dileptons_->column("l1_idasn", l1.idAssigned());
    nTuple_dileptons_->column("l1_idtru", l1.idTrue());
    nTuple_dileptons_->column("l1_idmom", l1.idMom());
    nTuple_dileptons_->column("l1_plab" , l1.p().rho());
    nTuple_dileptons_->column("l1_pcm"  , l1.pCm().rho());
    nTuple_dileptons_->column("l1_e_cm" , l1.pCm().e());
    nTuple_dileptons_->column("l1_pcm_x", l1.pCm().px());
    nTuple_dileptons_->column("l1_pcm_y", l1.pCm().py());
    nTuple_dileptons_->column("l1_pcm_z", l1.pCm().pz());
    nTuple_dileptons_->column("l1_ip_dr", l1.track().dr());
    nTuple_dileptons_->column("l1_ip_dz", l1.track().dz());
    nTuple_dileptons_->column("l1_svdr" , l1.svdRHits());
    nTuple_dileptons_->column("l1_svdz" , l1.svdZHits());

    nTuple_dileptons_->dumpData();
  }

  return;
}

// Specifies n-tuples to write and names for root histograms.
void
Adcab::hist_def()
{
  extern BelleTupleManager *BASF_Histogram;   // Define a BASF Histogram

  BelleTupleManager *tm = BASF_Histogram;
  const char *charged_particle_variables = "stm_no "
                                           "exp_no "
                                           "run_no "
                                           "evt_no "
                                           "is_mc "
                                           "is_cntnm "
                                           "cm_enrgy "
                                           "charge "
                                           "mass "
                                           "good_mu "
                                           "veto_mu "
                                           "good_el "
                                           "veto_el "
                                           "good_k "
                                           "pid_k_pi "
                                           "pid_k_pr "
                                           "id_asn "
                                           "id_tru "
                                           "id_mom "
                                           "eid_prob "
                                           "muid_prb "
                                           "muid_rto "
                                           "p_lb_mag "
                                           "p_cm_mag "
                                           "e_cm "
                                           "p_cm_x "
                                           "p_cm_y "
                                           "p_cm_z "
                                           "cos_pol "
                                           "ip_dr "
                                           "ip_dz "
                                           "svd_hitr "
                                           "svd_hitz";

  const char *dilepton_variables = "stm_no "
                                   "exp_no "
                                   "run_no "
                                   "evt_no "
                                   "is_mc "
                                   "is_cntnm "
                                   "fw_r2 "
                                   "hadronb "
                                   "cm_enrgy "
                                   "n_can_kp "
                                   "n_can_km "
                                   "n_tru_kp "
                                   "n_tru_km "
                                   "typ_asn "
                                   "typ_tru "
                                   "evt_sign "
                                   "cos_thta "
                                   "inv_mass "
                                   "l0_chrge " "l1_chrge "
                                   "l0_idasn " "l1_idasn "
                                   "l0_idtru " "l1_idtru "
                                   "l0_idmom " "l1_idmom "
                                   "l0_plab "  "l1_plab "
                                   "l0_pcm "   "l1_pcm "
                                   "l0_e_cm "  "l1_e_cm "
                                   "l0_pcm_x " "l1_pcm_x "
                                   "l0_pcm_y " "l1_pcm_y "
                                   "l0_pcm_z " "l1_pcm_z "
                                   "l0_ip_dr " "l1_ip_dr "
                                   "l0_ip_dz " "l1_ip_dz "
                                   "l0_svdr "  "l1_svdr "
                                   "l0_svdz "  "l1_svdz";

  nTuple_leptons_ = tm->ntuple("Leptons", charged_particle_variables, 1);
  nTuple_kaons_ = tm->ntuple("Kaons", charged_particle_variables, 2);
  nTuple_dileptons_ = tm->ntuple("Dilepton", dilepton_variables, 3);

  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
