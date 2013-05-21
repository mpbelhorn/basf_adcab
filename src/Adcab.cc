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
  static std::vector<Particle> lepton_candidates(5);
  static std::vector<Particle> kaon_candidates(10);
  static std::vector<Particle> phi_candidates(5);
  static std::vector<DileptonEvent> dilepton_event_candidates(10);
  lepton_candidates.clear();
  kaon_candidates.clear();
  phi_candidates.clear();
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

  // Populate the final state particle (FSP) candidate lists.
  for (MdstChargedIterator lepton_scan_iterator = first_mdst_charged;
      lepton_scan_iterator != last_mdst_charged; ++lepton_scan_iterator) {

    if (basf_parameter_verbose_log_ > 1) {
      cout << "    New charged Track" << endl;
    }
    // Alias the current particle as "charged_particle".
    const Mdst_charged &charged_particle = *lepton_scan_iterator;

    // Set up an instantance of UserInfo to store PID data.
    UserInfo pid_info;

    // Get final state particle (FSP) ID likelihoods and track quality.
    eid charged_particle_eid(charged_particle);
    Muid_mdst charged_particle_muid(charged_particle);
    pid_info.muonId(charged_particle_muid);
    pid_info.electronLikelihood(charged_particle_eid.prob(3, -1, 5));
    pid_info.kaonToPionLikelihood(pid_kaon_to_pi.prob(charged_particle));
    pid_info.kaonToProtonLikelihood(pid_kaon_to_pr.prob(charged_particle));
    pid_info.svdHits(charged_particle);

    if (basf_parameter_verbose_log_ > 2) {
      const Gen_hepevt & hep(get_hepevt(charged_particle));
      if (hep) {
        cout << "      PID=" << hep.idhep();
        if (hep.mother()) {
          cout << ", mom=" << hep.mother().idhep();
        }
        cout << endl;
      }
    }

    // Check that track has good likelihoods to be interesting FSP.
    bool good_muon_likelihood = pid_info.muonLikelihood() >= cuts.minMuidProb;
    bool good_electron_likelihood = pid_info.electronLikelihood() >= cuts.minEidProb;
    bool good_kaon_likelihood = (
        (pid_info.kaonToPionLikelihood() > cuts.minKaonToPionLikelihood) &&
        (pid_info.kaonToProtonLikelihood() > cuts.minKaonToProtonLikelihood));

    if (basf_parameter_verbose_log_ > 2) {
      cout << "      PID mu:" << pid_info.muonLikelihood() << "|"
                              << pid_info.klmSignature()
           << " el:" << pid_info.electronLikelihood()
           << " k-pi:" << pid_info.kaonToPionLikelihood()
           << " k-pr:" << pid_info.kaonToProtonLikelihood()
           << endl;
    }

    // Check that track has enough hits in the SVD for reliable vertexing.
    bool good_svd_electron = (
        (pid_info.svdRHits(0) >= cuts.minSvdRHits) &&
        (pid_info.svdZHits(0) >= cuts.minSvdZHits));
    bool good_svd_muon = (
        (pid_info.svdRHits(1) >= cuts.minSvdRHits) &&
        (pid_info.svdZHits(1) >= cuts.minSvdZHits));
    // bool good_svd_kaon = (
    //    (pid_info.svdRHits(3) >= cuts.minSvdRHits) &&
    //    (pid_info.svdZHits(3) >= cuts.minSvdZHits));

    if (basf_parameter_verbose_log_ > 2) {
      cout << "      SVD mu:" << pid_info.svdRHits(1) << "|"
                              << pid_info.svdZHits(1)
           << " el:" << pid_info.svdRHits(0) << "|" << pid_info.svdZHits(0)
           << " k:" << pid_info.svdRHits(3) << "|" << pid_info.svdZHits(3)
           << endl;
    }

    // Check that the KLM signature is good for muon candidates.
    bool good_klm_signature = (
        (pid_info.klmSignature() > 0) &&
        (pid_info.klmSignature() < cuts.maxKlmChi2PerHits));

    // Add up the cuts to determine species candidacy.
    bool good_muon = (good_muon_likelihood && good_svd_muon && good_klm_signature);
    bool good_electron = (good_electron_likelihood && good_svd_electron);
    bool good_kaon = (good_kaon_likelihood /*&& good_svd_kaon*/);

    // Create a particle instance for the candidate final state particle.
    double particle_charge(charged_particle.charge());
    Particle muon_candidate(charged_particle,
        particle_charge > 0 ? Ptype("MU+") : Ptype("MU-"));
    muon_candidate.userInfo(pid_info);
    setMCtruth(muon_candidate);
    UserInfo &muon_info = dynamic_cast<UserInfo&>(muon_candidate.userInfo());
    muon_info.svdHits(pid_info.svdRHits(1), pid_info.svdZHits(1));

    Particle electron_candidate(charged_particle,
        particle_charge > 0 ? Ptype("E+") : Ptype("E-"));
    electron_candidate.userInfo(pid_info);
    setMCtruth(electron_candidate);
    UserInfo &electron_info = dynamic_cast<UserInfo&>(electron_candidate.userInfo());
    electron_info.svdHits(pid_info.svdRHits(0), pid_info.svdZHits(0));

    Particle kaon_candidate(charged_particle,
        particle_charge > 0 ? Ptype("K+") : Ptype("K-"));
    kaon_candidate.userInfo(pid_info);
    setMCtruth(kaon_candidate);
    UserInfo &kaon_info = dynamic_cast<UserInfo&>(kaon_candidate.userInfo());
    kaon_info.svdHits(pid_info.svdRHits(3), pid_info.svdZHits(3));

    if (! (good_muon || good_electron || good_kaon)) {
      // Print diagnostic information to the log.
      if (basf_parameter_verbose_log_ > 1) {
        cout << "        Rejected" << endl;
      }
      continue;
    }

    TrackParameters muon_ip_profile(muon_candidate, interaction_point_);
    muon_info.ipDeltaR(muon_ip_profile.dr());
    muon_info.ipDeltaZ(muon_ip_profile.dz());

    TrackParameters electron_ip_profile(electron_candidate, interaction_point_);
    electron_info.ipDeltaR(electron_ip_profile.dr());
    electron_info.ipDeltaZ(electron_ip_profile.dz());

    TrackParameters kaon_ip_profile(kaon_candidate, interaction_point_);
    kaon_info.ipDeltaR(kaon_ip_profile.dr());
    kaon_info.ipDeltaZ(kaon_ip_profile.dz());

    // All leptons are considered prompt signal candidates unless shown to
    // be otherwise.
    bool prompt_muon = true;
    bool prompt_electron = true;

    // Cut on IP proximity.
    if (abs(muon_info.ipDeltaR()) > cuts.maxIpDr ||
        abs(muon_info.ipDeltaZ()) > cuts.maxIpDz) {
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut muon on IP dz or dr." << endl;
      }
      prompt_muon = false;
    }
    if (abs(electron_info.ipDeltaR()) > cuts.maxIpDr ||
        abs(electron_info.ipDeltaZ()) > cuts.maxIpDz) {
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut electron on IP dz or dr." << endl;
      }
      prompt_electron = false;
    }

    // Cut on barrel intersection.
    double muon_polar_angle_cosine = muon_candidate.p().cosTheta();
    if ((muon_polar_angle_cosine < cuts.minLeptonCosTheta) ||
        (muon_polar_angle_cosine > cuts.maxLeptonCosTheta)) {
      prompt_muon = false;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut muon on barrel intersection." << endl;
      }
    }
    double electron_polar_angle_cosine = muon_candidate.p().cosTheta();
    if ((electron_polar_angle_cosine < cuts.minLeptonCosTheta) ||
        (electron_polar_angle_cosine > cuts.maxLeptonCosTheta)) {
      prompt_electron = false;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut electron on barrel intersection." << endl;
      }
    }

    // Calculate the CM momentum, scale momentum to Y(5S) CM if necessary.
    if (basf_parameter_scale_momentum_) {
      // Set scaled P(CM) in particle info.
      // TODO - Possibly change lab momentum in particle instance.
      double muon_scale_factor(1.03291);
      HepLorentzVector scaled_muon_p(
          muon_candidate.p().vect() * muon_scale_factor,
          muon_candidate.p().t() * muon_scale_factor);
      muon_info.pCm(scaled_muon_p, cm_boost_);

      double electron_scale_factor(1.03290);
      HepLorentzVector scaled_electron_p(
          electron_candidate.p().vect() * electron_scale_factor,
          electron_candidate.p().t() * electron_scale_factor);
      electron_info.pCm(scaled_electron_p, cm_boost_);
    } else {
      muon_info.pCm(muon_candidate.p(), cm_boost_);
      electron_info.pCm(electron_candidate.p(), cm_boost_);
    }

    // CM momentum cut.
    if ((muon_info.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (muon_info.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      prompt_muon = false;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut muon on CM momentum." << endl;
      }
    }
    if ((electron_info.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (electron_info.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      prompt_electron = false;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "          Cut electron on CM momentum." << endl;
      }
    }

    // J/Psi veto.
    for (MdstChargedIterator mdst_charged_iterator = first_mdst_charged;
        mdst_charged_iterator != last_mdst_charged; ++mdst_charged_iterator) {
      const Mdst_charged &sister_mdst = *mdst_charged_iterator;

      // Reject case where pointers point to same object.
      if (charged_particle == sister_mdst) continue;

      // If allowing a charge bias and pair is SS, skip to the next pair.
      if (basf_parameter_allow_charge_bias_ &&
          (particle_charge == sister_mdst.charge())) {
        continue;
      }

      // Check to see if the candidate lepton and the other charged track
      //   make a good J/Psi. If candidate and sister are likely from a J/Psi
      //   (or pair production for the electron case), reject the candidate.
      Particle sister_electron(sister_mdst,
          sister_mdst.charge() > 0 ? Ptype("E+") : Ptype("E-"));
      double dielectron_invariant_mass = abs(
          (electron_candidate.p() + sister_electron.p()).m());
      if (dielectron_invariant_mass < cuts.minEPlusEMinusMass) {
        prompt_electron = false;
        if (basf_parameter_verbose_log_ > 1) {
          cout << "          Cut electron on pair-production veto." << endl;
        }
      }
      double dielectron_jpsi_mass_difference = dielectron_invariant_mass - cuts.massJPsi;
      if (cuts.minElElJPsiCandidate < dielectron_jpsi_mass_difference &&
          dielectron_jpsi_mass_difference < cuts.maxElElJPsiCandidate) {
        prompt_electron = false;
        if (basf_parameter_verbose_log_ > 1) {
          cout << "          Cut electron on J/Psi veto." << endl;
        }
      }

      Particle sister_muon(sister_mdst,
          sister_mdst.charge() > 0 ? Ptype("MU+") : Ptype("MU-"));
      double dimuon_invariant_mass = abs(
          (muon_candidate.p() + sister_muon.p()).m());
      double dimuon_jpsi_mass_difference = dimuon_invariant_mass - cuts.massJPsi;
      if (cuts.minMuMuJPsiCandidate < dimuon_jpsi_mass_difference &&
          dimuon_jpsi_mass_difference < cuts.maxMuMuJPsiCandidate) {
        prompt_muon = false;
        if (basf_parameter_verbose_log_ > 1) {
          cout << "          Cut muon on J/Psi veto." << endl;
        }
      }
    }
    // If still a good prompt candidate after above cuts, add to lepton_candidates.
    if (prompt_muon) {
      lepton_candidates.push_back(muon_candidate);
    } else if (prompt_electron) {
      lepton_candidates.push_back(electron_candidate);
    }

    // Regardless if the particle makes a good lepton candidate, if it is a
    // good K candidate, add it to the list of kaons.
    if (good_kaon) {
      kaon_candidates.push_back(kaon_candidate);
    }
    if (good_muon || good_electron || good_kaon) {
      nTuple_charged_->column("stm_no"  , basf_parameter_mc_stream_number_);
      nTuple_charged_->column("exp_no"  , experiment_number_);
      nTuple_charged_->column("run_no"  , run_number_);
      nTuple_charged_->column("evt_no"  , event_number_);
      nTuple_charged_->column("is_mc"   , flag_mc_);
      nTuple_charged_->column("is_cntnm", basf_parameter_is_continuum_);
      nTuple_charged_->column("cm_enrgy", beam_energy_cm_frame_);
      nTuple_charged_->column("cm_bst_x", cm_boost_.x());
      nTuple_charged_->column("cm_bst_y", cm_boost_.y());
      nTuple_charged_->column("cm_bst_z", cm_boost_.z());
      nTuple_charged_->column("charge"  , particle_charge);

      nTuple_charged_->column("mu_enrgy", muon_info.pCm().e());  // Scaled CM frame
      nTuple_charged_->column("mu_p_x"  , muon_info.pCm().px());
      nTuple_charged_->column("mu_p_y"  , muon_info.pCm().py());
      nTuple_charged_->column("mu_p_z"  , muon_info.pCm().pz());
      nTuple_charged_->column("el_enrgy", electron_info.pCm().e());  // Scaled CM frame
      nTuple_charged_->column("el_p_x"  , electron_info.pCm().px());
      nTuple_charged_->column("el_p_y"  , electron_info.pCm().py());
      nTuple_charged_->column("el_p_z"  , electron_info.pCm().pz());
      nTuple_charged_->column("k_enrgy" , kaon_candidate.p().e());  // Lab frame
      nTuple_charged_->column("k_p_x"   , kaon_candidate.p().px());
      nTuple_charged_->column("k_p_y"   , kaon_candidate.p().py());
      nTuple_charged_->column("k_p_z"   , kaon_candidate.p().pz());

      nTuple_charged_->column("id_tru"  , IDhep(muon_candidate)); // Same for all candidates.
      nTuple_charged_->column("id_mom"  , IDmom(muon_candidate)); // Same for all candidates.
      nTuple_charged_->column("good_mu" , good_muon);
      nTuple_charged_->column("good_el" , good_electron);
      nTuple_charged_->column("good_k"  , good_kaon);
      nTuple_charged_->column("prmpt_mu", prompt_muon);
      nTuple_charged_->column("prmpt_el", prompt_electron);
      nTuple_charged_->column("muid_prb", pid_info.muonLikelihood());
      nTuple_charged_->column("eid_prob", pid_info.electronLikelihood());
      nTuple_charged_->column("pid_k_pi", pid_info.kaonToPionLikelihood());
      nTuple_charged_->column("pid_k_pr", pid_info.kaonToProtonLikelihood());
      nTuple_charged_->column("muid_rto", pid_info.klmSignature());
      nTuple_charged_->column("mu_svd_r", muon_info.svdRHits());
      nTuple_charged_->column("mu_svd_z", muon_info.svdZHits());
      nTuple_charged_->column("el_svd_r", electron_info.svdRHits());
      nTuple_charged_->column("el_svd_z", electron_info.svdZHits());
      nTuple_charged_->column("k_svd_r" , kaon_info.svdRHits());
      nTuple_charged_->column("k_svd_z" , kaon_info.svdZHits());
      nTuple_charged_->column("mu_ip_dr", muon_info.ipDeltaR());
      nTuple_charged_->column("mu_ip_dz", muon_info.ipDeltaZ());
      nTuple_charged_->column("el_ip_dr", electron_info.ipDeltaR());
      nTuple_charged_->column("el_ip_dz", electron_info.ipDeltaZ());
      nTuple_charged_->column("k_ip_dr" , kaon_info.ipDeltaR());
      nTuple_charged_->column("k_ip_dz" , kaon_info.ipDeltaZ());
      nTuple_charged_->dumpData();
    }
  }

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "    Passed track selection." << endl;
    cout << "      Lepton candidates: " << lepton_candidates.size() << endl;
    cout << "      Kaon candidates: " << kaon_candidates.size() << endl;
  }

  // Find Phi canididates.
  for (ParticleIterator kaon1 = kaon_candidates.begin();
      kaon1 != kaon_candidates.end(); ++kaon1) {
    for (ParticleIterator kaon2 = kaon1;
        kaon2 != kaon_candidates.end(); ++kaon2) {
      if (kaon1 == kaon2) continue;
      double kaon1_charge = kaon1->charge();
      double kaon2_charge = kaon2->charge();
      if (kaon1_charge == kaon2_charge) continue;
      Particle &kaon_minus = kaon1_charge > 0 ? *kaon2 : *kaon1;
      Particle &kaon_plus = kaon1_charge > 0 ? *kaon1 : *kaon2;
      HepLorentzVector phi_momentum = kaon_minus.p() + kaon_plus.p();
      Particle phi_candidate(phi_momentum, Ptype("PHI"));
      phi_candidate.relation().append(kaon_minus);
      phi_candidate.relation().append(kaon_plus);
      int error = fitPhiVertex(phi_candidate);
      setMCtruth(phi_candidate);
      UserInfo &phi_info = dynamic_cast<UserInfo&>(phi_candidate.userInfo());
      if (!error) {
        if (basf_parameter_verbose_log_ > 2) {
          cout << "Phi vertex chisq/dof: " << phi_info.chisq()
               << "/" << phi_info.ndf()
               << " | truth: " << IDhep(phi_candidate)
               << "[" << phi_info.genHepevtLink()
               << "] (" << IDhep(phi_candidate.child(0))
               << " : " << IDhep(phi_candidate.child(1)) << ")"
               << endl;
        }

        phi_candidates.push_back(phi_candidate);
        nTuple_phi_->column("cm_enrgy", beam_energy_cm_frame_);
        nTuple_phi_->column("cm_bst_x", cm_boost_.x());
        nTuple_phi_->column("cm_bst_y", cm_boost_.y());
        nTuple_phi_->column("cm_bst_z", cm_boost_.z());
        nTuple_phi_->column("km_e"    , kaon_minus.p().e());
        nTuple_phi_->column("km_px"   , kaon_minus.p().x());
        nTuple_phi_->column("km_py"   , kaon_minus.p().y());
        nTuple_phi_->column("km_pz"   , kaon_minus.p().z());
        nTuple_phi_->column("km_pid"  , IDhep(kaon_minus));
        nTuple_phi_->column("kp_e"    , kaon_plus.p().e());
        nTuple_phi_->column("kp_px"   , kaon_plus.p().x());
        nTuple_phi_->column("kp_py"   , kaon_plus.p().y());
        nTuple_phi_->column("kp_pz"   , kaon_plus.p().z());
        nTuple_phi_->column("kp_pid"  , IDhep(kaon_plus));
        nTuple_phi_->column("ph_pid"  , IDhep(phi_candidate));
        nTuple_phi_->column("ph_recon", phi_info.genHepevtLink());
        nTuple_phi_->column("decay_x" , phi_candidate.momentum().decayVertex().x());
        nTuple_phi_->column("decay_y" , phi_candidate.momentum().decayVertex().y());
        nTuple_phi_->column("decay_z" , phi_candidate.momentum().decayVertex().z());
        nTuple_phi_->column("vx_chi2" , phi_info.chisq());
        nTuple_phi_->column("vx_dof"  , phi_info.ndf());
        nTuple_phi_->column("trk_prx" , 0); // TODO - Get closest approach of K tracks.
        nTuple_phi_->dumpData();
      }
    }
  }

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "    Phi selection." << endl;
    cout << "      phi candidates: " << phi_candidates.size() << endl;
  }

  // Find good dilepton event candidates.
  // Check all lepton pairs.
  for (ParticleIterator outer_lepton = lepton_candidates.begin();
      outer_lepton != lepton_candidates.end(); ++outer_lepton) {
    UserInfo &outer_lepton_info =
        dynamic_cast<UserInfo&>(outer_lepton->userInfo());
    for (ParticleIterator inner_lepton = outer_lepton;
        inner_lepton != lepton_candidates.end(); ++inner_lepton) {
      // Exclude the case where both iterators point to the same particle.
      if (outer_lepton == inner_lepton) continue;
      if (basf_parameter_verbose_log_ > 1) {
        cout << "    New dilepton event candidate" << endl;
      }

      UserInfo &inner_lepton_info =
          dynamic_cast<UserInfo&>(inner_lepton->userInfo());

      // Determine higher momentum lepton and add it to an event candidate.
      // Reference the higher momentum lepton as "lepton0"
      //  and the lower momentum lepton as "lepton1". First, assume the outer
      //  loop lepton has the greater momentum, then check that.
      Particle *greater_p_lepton;
      Particle *lesser_p_lepton;
      if (inner_lepton_info.pCm().mag() > outer_lepton_info.pCm().mag()) {
        greater_p_lepton = &(*inner_lepton);
        lesser_p_lepton = &(*outer_lepton);
      } else {
        greater_p_lepton = &(*outer_lepton);
        lesser_p_lepton = &(*inner_lepton);
      }
      Particle &l0 = *greater_p_lepton;
      Particle &l1 = *lesser_p_lepton;
      UserInfo &l0_info = dynamic_cast<UserInfo&>(l0.userInfo());
      UserInfo &l1_info = dynamic_cast<UserInfo&>(l1.userInfo());

      DileptonEvent event_candidate(l0, l1);

      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosine_cm_opening_angle = event_candidate.cosThetaLL();
      if (cosine_cm_opening_angle < cuts.minCosThetaLLCm ||
          cosine_cm_opening_angle > cuts.maxCosThetaLLCm) {
        continue;
      }

      int l0_id(l0.mdstCharged().get_ID());
      int l1_id(l1.mdstCharged().get_ID());
      int phi_candidates_in_mass_range = 0;
      for (ParticleIterator phi_candidate = phi_candidates.begin();
          phi_candidate != phi_candidates.end(); ++phi_candidate) {
        int kaon_minus_id(phi_candidate->relation().child(0).mdstCharged().get_ID());
        int kaon_plus_id(phi_candidate->relation().child(1).mdstCharged().get_ID());
        if (kaon_minus_id == l0_id ||
            kaon_minus_id == l1_id ||
            kaon_plus_id == l0_id ||
            kaon_plus_id == l1_id) {
          // Bad phi+dilepton: double counted tracks.
          continue;
        }
        if (phi_candidate->p().m() < 1.045) {
          phi_candidates_in_mass_range++;
        }
      }

      if (basf_parameter_verbose_log_ > 1) {
        cout << "        Found good dilepton event candidate." << endl;
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
      nTuple_dileptons_->column("cm_bst_x", cm_boost_.x());
      nTuple_dileptons_->column("cm_bst_y", cm_boost_.y());
      nTuple_dileptons_->column("cm_bst_z", cm_boost_.z());

      // Write dilepton-level data to the n-tuple.
      nTuple_dileptons_->column("cos_thta", event_candidate.cosThetaLL());
      nTuple_dileptons_->column("phi_cnds", phi_candidates.size());
      nTuple_dileptons_->column("good_phi", phi_candidates_in_mass_range);

      // Write lepton-level data to the n-tuple.
      nTuple_dileptons_->column("l0_chrge", l0.charge());
      nTuple_dileptons_->column("l0_idasn", l0.pType().lund());
      nTuple_dileptons_->column("l0_idtru", IDhep(l0));
      nTuple_dileptons_->column("l0_idmom", IDmom(l0));
      nTuple_dileptons_->column("l0_e_cm" , l0_info.pCm().e());
      nTuple_dileptons_->column("l0_pcm_x", l0_info.pCm().px());
      nTuple_dileptons_->column("l0_pcm_y", l0_info.pCm().py());
      nTuple_dileptons_->column("l0_pcm_z", l0_info.pCm().pz());
      nTuple_dileptons_->column("l0_ip_dr", l0_info.ipDeltaR());
      nTuple_dileptons_->column("l0_ip_dz", l0_info.ipDeltaZ());

      nTuple_dileptons_->column("l1_chrge", l1.charge());
      nTuple_dileptons_->column("l1_idasn", l1.pType().lund());
      nTuple_dileptons_->column("l1_idtru", IDhep(l1));
      nTuple_dileptons_->column("l1_idmom", IDmom(l1));
      nTuple_dileptons_->column("l1_e_cm" , l1_info.pCm().e());
      nTuple_dileptons_->column("l1_pcm_x", l1_info.pCm().px());
      nTuple_dileptons_->column("l1_pcm_y", l1_info.pCm().py());
      nTuple_dileptons_->column("l1_pcm_z", l1_info.pCm().pz());
      nTuple_dileptons_->column("l1_ip_dr", l1_info.ipDeltaR());
      nTuple_dileptons_->column("l1_ip_dz", l1_info.ipDeltaZ());

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
  const char *charged_particle_variables =
      "stm_no "
      "exp_no "
      "run_no "
      "evt_no "
      "is_mc "
      "is_cntnm "
      "cm_enrgy "
      "cm_bst_x "
      "cm_bst_y "
      "cm_bst_z "
      "charge "
      "mu_enrgy "
      "mu_p_x "
      "mu_p_y "
      "mu_p_z "
      "el_enrgy "
      "el_p_x "
      "el_p_y "
      "el_p_z "
      "k_enrgy "
      "k_p_x "
      "k_p_y "
      "k_p_z "
      "id_tru "
      "id_mom "
      "good_mu "
      "good_el "
      "good_k "
      "prmpt_mu "
      "prmpt_el "
      "muid_prb "
      "eid_prob "
      "pid_k_pi "
      "pid_k_pr "
      "muid_rto "
      "mu_svd_r "
      "mu_svd_z "
      "el_svd_r "
      "el_svd_z "
      "k_svd_r "
      "k_svd_z "
      "mu_ip_dr "
      "mu_ip_dz "
      "el_ip_dr "
      "el_ip_dz "
      "k_ip_dr "
      "k_ip_dz";

  const char *phi_variables =
      "cm_enrgy "
      "cm_bst_x "
      "cm_bst_y "
      "cm_bst_z "
      "km_e "     "kp_e "
      "km_px "    "kp_px "
      "km_py "    "kp_py "
      "km_pz "    "kp_pz "
      "km_pid "   "kp_pid "
      "ph_pid "
      "ph_recon "
      "decay_x "
      "decay_y "
      "decay_z "
      "vx_chi2 "
      // "vx_dof "
      "trk_prx ";

  const char *dilepton_variables =
      "stm_no "
      "exp_no "
      "run_no "
      "evt_no "
      "is_mc "
      "is_cntnm "
      "fw_r2 "
      "hadronb "
      "cm_enrgy "
      "cm_bst_x "
      "cm_bst_y "
      "cm_bst_z "
      "cos_thta "
      "phi_cnds "
      "good_phi "
      "l0_chrge " "l1_chrge "
      "l0_idasn " "l1_idasn "
      "l0_idtru " "l1_idtru "
      "l0_idmom " "l1_idmom "
      "l0_e_cm "  "l1_e_cm "
      "l0_pcm_x " "l1_pcm_x "
      "l0_pcm_y " "l1_pcm_y "
      "l0_pcm_z " "l1_pcm_z "
      "l0_ip_dr " "l1_ip_dr "
      "l0_ip_dz " "l1_ip_dz";

  nTuple_charged_ = tm->ntuple("Charged", charged_particle_variables, 1);
  nTuple_phi_ = tm->ntuple("Phi", phi_variables, 2);
  nTuple_dileptons_ = tm->ntuple("Dilepton", dilepton_variables, 3);

  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
