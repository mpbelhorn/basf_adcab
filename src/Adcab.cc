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
  particle_k_minus_ = (Ptype("K-"));
  particle_k_plus_ = (Ptype("K+"));

  // Initialize BASF parameters.
  basf_parameter_allow_charge_bias_ = 0;
  basf_parameter_verbose_log_ = 0;
  basf_parameter_is_continuum_ = 0;
  basf_parameter_mc_stream_number_ = 0;

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

  // Instantiate PID classes. These are constant and
  //     should only be initialized once.
  static atc_pid pid_kaon_to_pi(3, 1, 5, 3, 2);
  static atc_pid pid_kaon_to_pr(3, 1, 5, 3, 4);

  // Define lists (vector template) to store event particles.
  // Need a list for all mother and daughter particle species.
  // Note that these vectors are static and will persist until BASF is closed.
  //     They therefore must(!) be cleared for each call of Adcab::event().
  if (basf_parameter_verbose_log_) {
    cout << "  Initializing particle containers." << endl;
  }
  static std::vector<LeptonCandidate> lepton_candidates(5);
  static std::vector<Particle> kaon_candidates(10);
  static std::vector<DileptonEvent> dilepton_event_candidates(10);
  lepton_candidates.clear();
  kaon_candidates.clear();
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

  // Populate the lepton candidate lists.
  for (MdstChargedIterator lepton_scan_iterator = first_mdst_charged;
      lepton_scan_iterator != last_mdst_charged; ++lepton_scan_iterator) {

    if (basf_parameter_verbose_log_) {
      cout << "  >>> NEW charged Track <<<" << endl;
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

    if (basf_parameter_verbose_log_) {
      cout << "    COMPLETED["<< good_muon << good_electron << good_kaon << "] "
           << "PID check." << endl;
    }
    if (!(good_muon || good_electron || good_kaon)) continue;

    // Create particle class objects of each species.
    double electric_charge = charged_particle.charge();
    Particle muon_particle(charged_particle,
        electric_charge > 0 ? particle_mu_plus_ : particle_mu_minus_);
    Particle electron_particle(charged_particle,
        electric_charge > 0 ? particle_e_plus_ : particle_e_minus_);
    LeptonCandidate muon_candidate(muon_particle, cm_boost_);
    LeptonCandidate electron_candidate(electron_particle, cm_boost_);

    // Cut on IP dr and dz and SVD hits.
    // This is to make sure that particles were created near the IP.
    // Seems there is no difference in dr or dz using either mass hypothesis = 1
    //     (muon) vs mass hypothesis = 0 (electron). Mass hypothesis = n (!=0,1)
    //     does produce different dr and dz values.
    IpParameters ip_parameters(charged_particle, interaction_point_, 1);
    if ((abs(ip_parameters.dr()) > cuts.maxIpDr) ||
        (abs(ip_parameters.dz()) > cuts.maxIpDz) ||
        (ip_parameters.svdHitsR() < cuts.minSvdRHits) ||
        (ip_parameters.svdHitsZ() < cuts.minSvdZHits)) {
      good_muon = false;
      good_electron = false;
    }
    if (basf_parameter_verbose_log_) {
      cout << "    COMPLETED["<< good_muon << good_electron << good_kaon << "] "
           << "dr/dz and SVD check." << endl;
    }
    if (!(good_muon || good_electron || good_kaon)) continue;

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
    if (basf_parameter_verbose_log_) {
      cout << "    COMPLETED["<< good_muon << good_electron << good_kaon << "] "
           << "Barrel intersection check." << endl;
    }
    if (!(good_muon || good_electron || good_kaon)) continue;

    // Check if CM momentum is good assuming muon and electron masses.
    if ((muon_candidate.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (muon_candidate.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      good_muon = false;
    }
    if ((electron_candidate.pCm().rho() < cuts.minLeptonMomentumCm) ||
        (electron_candidate.pCm().rho() > cuts.maxLeptonMomentumCm)) {
      good_electron = false;
    }
    if (basf_parameter_verbose_log_) {
      cout << "    COMPLETED["<< good_muon << good_electron << good_kaon << "] "
           << "Passed CM-frame momentum check." << endl;
    }
    if (!(good_muon || good_electron || good_kaon)) continue;

    // Remove possible pair production and J/Psi daughters.
    for (MdstChargedIterator jpsi_pair_iterator = first_mdst_charged;
        (good_muon || good_electron) &&
        (jpsi_pair_iterator != last_mdst_charged); ++jpsi_pair_iterator) {
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
    if (basf_parameter_verbose_log_) {
      cout << "    COMPLETED["<< good_muon << good_electron << good_kaon << "] "
           << "J/Psi and pair production veto." << endl;
    }
    if (!(good_muon || good_electron || good_kaon)) continue;

    // While the candidate is still good, add it to the lepton list. Only add
    //     the candidate to the lepton list once! If it is a good muon
    //     candidate, then consider it a muon, otherwise it is deemed an
    //     electron, but flag the verbose log if the particle passes as both.
    LeptonCandidate *good_lepton = NULL;
    if (good_muon) {
      setMCtruth(muon_candidate.lepton());
      good_lepton = &muon_candidate;
    } else if (good_electron) {
      setMCtruth(electron_candidate.lepton());
      good_lepton = &electron_candidate;
    }

    // Dump lepton candidate information to the ntuple.
    if (good_lepton) {
      if (basf_parameter_verbose_log_) {
        cout << "    Committing lepton to ntuple." << endl;
      }
      LeptonCandidate lepton((*good_lepton));
      lepton_candidates.push_back(lepton);

      nTuple_leptons_->column("charge"   , electric_charge);
      nTuple_leptons_->column("mass"     , lepton.p().mag());
      nTuple_leptons_->column("good_mu"  , good_muon);
      nTuple_leptons_->column("good_el"  , good_electron);
      nTuple_leptons_->column("good_k"   , good_kaon);
      nTuple_leptons_->column("id_asn"   , lepton.idAssigned());
      nTuple_leptons_->column("id_tru"   , lepton.idTrue());
      nTuple_leptons_->column("id_mom"   , lepton.idMom());
      nTuple_leptons_->column("eid_prob" , eid_probability);
      nTuple_leptons_->column("muid_prb" , muid_probability);
      nTuple_leptons_->column("muid_rto" , lepton.klmChi2PerHits());
      nTuple_leptons_->column("p_lb_mag" , lepton.p().rho());
      nTuple_leptons_->column("p_cm_mag" , lepton.pCm().rho());
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

    if (good_kaon && !good_lepton) {
      if (basf_parameter_verbose_log_) {
        cout << "    Committing kaon to ntuple." << endl;
      }
      // Treat the particle as a charged kaon.
      Particle kaon_particle(charged_particle,
          electric_charge > 0 ? particle_k_plus_ : particle_k_minus_);
      setMCtruth(kaon_particle);
      // TODO - Need to fix this by either implementing a KaonCandidate
      //   class or generalizing the LeptonCandidate class.
      LeptonCandidate kaon(kaon_particle, cm_boost_);
      kaon_candidates.push_back(kaon);

      nTuple_kaons_->column("charge"   , electric_charge);
      nTuple_kaons_->column("mass"     , kaon.p().mag());
      nTuple_kaons_->column("good_mu"  , good_muon);
      nTuple_kaons_->column("good_el"  , good_electron);
      nTuple_kaons_->column("good_k"   , good_kaon);
      nTuple_kaons_->column("id_asn"   , kaon.idAssigned());
      nTuple_kaons_->column("id_tru"   , kaon.idTrue());
      nTuple_kaons_->column("id_mom"   , kaon.idMom());
      nTuple_kaons_->column("pid_k_pi" , kaon_to_pion_likelihood);
      nTuple_kaons_->column("pid_k_pr" , kaon_to_proton_likelihood);
      nTuple_kaons_->column("p_lb_mag" , kaon.p().rho());
      nTuple_kaons_->column("p_cm_mag" , kaon.pCm().rho());
      nTuple_kaons_->column("e_cm"     , kaon.pCm().e());
      nTuple_kaons_->column("p_cm_x"   , kaon.pCm().px());
      nTuple_kaons_->column("p_cm_y"   , kaon.pCm().py());
      nTuple_kaons_->column("p_cm_z"   , kaon.pCm().pz());
      nTuple_kaons_->column("cos_pol"  , kaon.p().cosTheta());
      nTuple_kaons_->column("ip_dr"    , ip_parameters.dr());
      nTuple_kaons_->column("ip_dz"    , ip_parameters.dz());
      nTuple_kaons_->column("svd_hitr" , ip_parameters.svdHitsR());
      nTuple_kaons_->column("svd_hitz" , ip_parameters.svdHitsZ());
      nTuple_kaons_->dumpData();
    }
  }

  // Print diagnostic information to the log.
  if (basf_parameter_verbose_log_) {
    cout << "  Passed lepton candidate selection." << endl;
    cout << "    Total lepton candidates: " << lepton_candidates.size() << endl;
    cout << "  Searching for dilepton candidates." << endl;
  }
 
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
      LeptonCandidate &l0 = *lepton0;
      LeptonCandidate &l1 = *lepton1;
      
      DileptonEvent event_candidate(l0, l1);

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
      
      IpParameters l0_ip_parameters(l0.lepton().relation().mdstCharged(),
          interaction_point_, l0.massHypothesis());
      IpParameters l1_ip_parameters(l1.lepton().relation().mdstCharged(),
          interaction_point_, l1.massHypothesis());

      // Write run-level data to the n-tuple.
      // Column names can be no greater than eight (8) characters long.
      nTuple_dileptons_->column("exp_no"  , experiment_number_);
      nTuple_dileptons_->column("run_no"  , run_number_);
      nTuple_dileptons_->column("stream"  , basf_parameter_mc_stream_number_);
      nTuple_dileptons_->column("is_mc"   , flag_mc_);
      nTuple_dileptons_->column("is_cntnm", basf_parameter_is_continuum_);

      // Write event-level data to the n-tuple.
      nTuple_dileptons_->column("evt_no"  , event_number_);
      nTuple_dileptons_->column("fw_r2"   , fox_wolfram_r2);
      nTuple_dileptons_->column("hadronb" , hadronb_code);

      // Write event-candidate data to the n-tuple.
      nTuple_dileptons_->column("evt_type", event_candidate.eventType());
      nTuple_dileptons_->column("evt_sign", event_candidate.eventSign());
      nTuple_dileptons_->column("llcostha", event_candidate.cosThetaLL());
      nTuple_dileptons_->column("n_kaons" , kaon_candidates.size());

      nTuple_dileptons_->column("l0_chrge", l0.lepton().charge());
      nTuple_dileptons_->column("l0_mass" , l0.p().mag());
      nTuple_dileptons_->column("l0_idasn", l0.idAssigned());
      nTuple_dileptons_->column("l0_idtru", l0.idTrue());
      nTuple_dileptons_->column("l0_idmom", l0.idMom());
      nTuple_dileptons_->column("l0_eidp" , l0.electronProbability());
      nTuple_dileptons_->column("l0_muidp", l0.muonProbability());
      nTuple_dileptons_->column("l0_muidr", l0.klmChi2PerHits());
      nTuple_dileptons_->column("l0_plab" , l0.p().rho());
      nTuple_dileptons_->column("l0_pcm"  , l0.pCm().rho());
      nTuple_dileptons_->column("l0_e_cm" , l0.pCm().e());
      nTuple_dileptons_->column("l0_px_cm", l0.pCm().px());
      nTuple_dileptons_->column("l0_py_cm", l0.pCm().py());
      nTuple_dileptons_->column("l0_pz_cm", l0.pCm().pz());
      nTuple_dileptons_->column("l0_cospl", l0.p().cosTheta());
      nTuple_dileptons_->column("l0_ip_dr", l0_ip_parameters.dr());
      nTuple_dileptons_->column("l0_ip_dz", l0_ip_parameters.dz());
      nTuple_dileptons_->column("l0_svd_r", l0_ip_parameters.svdHitsR());
      nTuple_dileptons_->column("l0_svd_z", l0_ip_parameters.svdHitsZ());

      nTuple_dileptons_->column("l1_chrge", l1.lepton().charge());
      nTuple_dileptons_->column("l1_mass" , l1.p().mag());
      nTuple_dileptons_->column("l1_idasn", l1.idAssigned());
      nTuple_dileptons_->column("l1_idtru", l1.idTrue());
      nTuple_dileptons_->column("l1_idmom", l1.idMom());
      nTuple_dileptons_->column("l1_eidp" , l1.electronProbability());
      nTuple_dileptons_->column("l1_muidp", l1.muonProbability());
      nTuple_dileptons_->column("l1_muidr", l1.klmChi2PerHits());
      nTuple_dileptons_->column("l1_plab" , l1.p().rho());
      nTuple_dileptons_->column("l1_pcm"  , l1.pCm().rho());
      nTuple_dileptons_->column("l1_e_cm" , l1.pCm().e());
      nTuple_dileptons_->column("l1_px_cm", l1.pCm().px());
      nTuple_dileptons_->column("l1_py_cm", l1.pCm().py());
      nTuple_dileptons_->column("l1_pz_cm", l1.pCm().pz());
      nTuple_dileptons_->column("l1_cospl", l1.p().cosTheta());
      nTuple_dileptons_->column("l1_ip_dr", l1_ip_parameters.dr());
      nTuple_dileptons_->column("l1_ip_dz", l1_ip_parameters.dz());
      nTuple_dileptons_->column("l1_svd_r", l1_ip_parameters.svdHitsR());
      nTuple_dileptons_->column("l1_svd_z", l1_ip_parameters.svdHitsZ());
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
  const char *run_variables = "exp_no "
                              "run_no "
                              "is_mc "
                              "is_cntnm "
                              "boost_x "
                              "boost_y "
                              "boost_z";

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
                               "good_mu "
                               "good_el "
                               "good_k "
                               "id_asn "
                               "id_tru "
                               "id_mom "
                               "pid_k_pi "
                               "pid_k_pr "
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

  const char *dilepton_variables = "exp_no "
                                   "run_no "
                                   "stream "
                                   "is_mc "
                                   "is_cntnm "
                                   "evt_no "
                                   "fw_r2 "
                                   "hadronb "
                                   "evt_type "
                                   "evt_sign "
                                   "n_kaons "
                                   "llcostha "
                                   "l0_chrge "
                                   "l0_mass "
                                   "l0_idasn "
                                   "l0_idtru "
                                   "l0_idmom "
                                   "l0_eidp "
                                   "l0_muidp "
                                   "l0_muidr "
                                   "l0_plab "
                                   "l0_pcm "
                                   "l0_e_cm "
                                   "l0_px_cm "
                                   "l0_py_cm "
                                   "l0_pz_cm "
                                   "l0_cospl "
                                   "l0_ip_dr "
                                   "l0_ip_dz "
                                   "l0_svd_r "
                                   "l0_svd_z"
                                   "l1_chrge "
                                   "l1_mass "
                                   "l1_idasn "
                                   "l1_idtru "
                                   "l1_idmom "
                                   "l1_eidp "
                                   "l1_muidp "
                                   "l1_muidr "
                                   "l1_plab "
                                   "l1_pcm "
                                   "l1_e_cm "
                                   "l1_px_cm "
                                   "l1_py_cm "
                                   "l1_pz_cm "
                                   "l1_cospl "
                                   "l1_ip_dr "
                                   "l1_ip_dz "
                                   "l1_svd_r "
                                   "l1_svd_z";

  nTuple_leptons_ = tm->ntuple("Leptons", lepton_variables, 1);
  nTuple_kaons_ = tm->ntuple("Kaons", kaon_variables, 2);
  nTuple_dileptons_ = tm->ntuple("Dilepton", dilepton_variables, 3);
  
  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
