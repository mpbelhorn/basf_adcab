//______________________________________________________________________________
// Filename: Adcab.cc
// Version: 2011.02.09.A
// Author: M.P. Belhorn
// Original Date: 2010-07-19
// Description: Analysis of the (A)nomalous (D)ilepton (C)harge (A)symmetry in
//   (B)s0 decays.
//______________________________________________________________________________

#include "Adcab.h"              // Adcab Analysis header.

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
  Module_descr *dscr = new Module_descr( "Adcab", module );
  
  // Set up module parameters.
  dscr->define_param ( "JPsi_Veto_OS_Only",
      "Apply J/Psi veto to opposite-sign pairs only",
      &module->flag_jpsi_veto_os_only);

  // Provide path to pass paramaters to BeamEnergy class.
  BeamEnergy::define_global( dscr );
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
  ptypeEMinus      = ( Ptype( "E-" ) );
  ptypeEPlus       = ( Ptype( "E+" ) );
  ptypeNuE         = ( Ptype( "NUE" ) );
  ptypeNuEBar      = ( Ptype( "NUEB" ) );
  ptypeMuMinus     = ( Ptype( "MU-" ) );
  ptypeMuPlus      = ( Ptype( "MU+" ) );
  ptypeNuMu        = ( Ptype( "NUM" ) );
  ptypeNuMuBar     = ( Ptype( "NUMB" ) );
  ptypeKMinus      = ( Ptype( "K-" ) );
  ptypeKPlus       = ( Ptype( "K+" ) );
  ptypeDsMinus     = ( Ptype( "DS-" ) );
  ptypeDsPlus      = ( Ptype( "DS+" ) );
  ptypeDsStarMinus = ( Ptype( "DS*-" ) );
  ptypeDsStarPlus  = ( Ptype( "DS*+" ) );

  // Initialize BASF parameters.
  flag_jpsi_veto_os_only = 1;
  
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
  std::cout
      << "\n\n"
      << "____________________________________________________________\n"
      << " Adcab Analysis Module terminated successfully \n\n"
      << std::endl;

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
  experimentNumber = 0;
  runNumber = 0;
  numberOfEvents = 0;
  numDileptonEvents = 0;
  
  // Set Diagnostic Variables to 0.
  num_bs_after_lepton_level.first = 0;
  num_bs_after_lepton_level.second = 0;
  num_bs_after_pair_removal.first = 0;
  num_bs_after_pair_removal.second = 0;
  num_bs_at_after_event_selection.first = 0;
  num_bs_at_after_event_selection.second = 0;

  // Set default flags.
  flagMC = false;
  flagERROR = false;
  flagVERBOSELOG = false;
  flagSELECTBESTCANDIDATE = false;
  
  // Set interaction point and error to default values.
  ip = HepPoint3D( 0, 0, 0 );
  ipErr = HepSymMatrix( 3, 0 );
  ipUsable = 0;

  // Get BELLE_EVENT runhead manager. This is the data stored in the 
  //   belletdf.tdf panther table.
  Belle_runhead_Manager &runhead_manager = Belle_runhead_Manager::get_manager();
  Belle_runhead &runhead = runhead_manager( ( Panther_ID ) 1 );
  experimentNumber = runhead.ExpNo();
  runNumber = runhead.RunNo();
  
  // Initialise BeamEnergy class.
  BeamEnergy::begin_run();
  
  // Set beam energy related class variables.
  beamEnergyCMFrame = BeamEnergy::E_beam_corr();  
  beamEnergyError   = BeamEnergy::E_beam_err();
  lerBeamEnergy     = BeamEnergy::E_LER();
  herBeamEnergy     = BeamEnergy::E_HER();
  kekbBeamEnergy    = BeamEnergy::E_beam_orig();
  kekbLerBeamEnergy = BeamEnergy::E_LER_orig();
  kekbHerBeamEnergy = BeamEnergy::E_HER_orig();
  beamCrossingAngle = BeamEnergy::Cross_angle();
  cmBoostVector     = -BeamEnergy::CMBoost();

  // Initialize PID functions.
  eid::init_data();
  
  // Get interaction point profile data from $BELLE_POSTGRES_SERVER. 
  IpProfile::begin_run();
  
  // Set interaction point and error to run values.
  if ( IpProfile::usable() ) {
    ip = IpProfile::e_position();
    ipErr = IpProfile::e_position_err_b_life_smeared();
    ipUsable = 1;
  } else {
    ip = HepPoint3D( 0, 0, 0 );
    ipErr = HepSymMatrix( 3, 0 );
  }
  
  // Print run information to the log.
  std::cout
      << "\n\n"
      << "____________________________________________________________\n"
      << " New Run: "
      << "Experiment " << experimentNumber 
      << ", Run " << runNumber 
      << "\n" 
      << std::endl;

  if ( runhead.ExpMC() == 1 ) {
    flagMC = false; // Set Data type flag to Real Data.
    std::cout << " Data is Real." << std::endl;
  } else {
    flagMC = true;  // Set Data type flag to Monte Carlo.
    std::cout << " Data is Monte Carlo." << std::endl;
  }
  
  std::cout
      << " Actual Beam Energy: " << beamEnergyCMFrame 
      << " +/- " << beamEnergyError << " GeV\n"
      << " Reported Beam Energy: " << kekbBeamEnergy << " GeV\n"
      << " BE Class cmBoostVector: " << cmBoostVector << "\n"
      << "\n"
      << " BASF Parameter Flags Settings:\n"
      << "   flag_jpsi_veto_os_only = " << flag_jpsi_veto_os_only << "\n"
      << "____________________________________________________________\n"
      << std::endl;

  return;
}

// This function is run once at the end of a data run. It does nothing but
//   write a message to somewhere. This message does not appear in log
//   although MPB believes it should.
void
Adcab::end_run(BelleEvent* evptr, int *status )
{ 
  (void)evptr;
  (void)status;
  
  std::cout
      << "\n\n"
      << "***************************************************\n"
      << "* WHERE IS THIS FUNCTION end_run() EXECUTED?!?!?! *\n"
      << "***************************************************\n"
      << std::endl;

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
  eventNumber = EvMgr[0].EvtNo() & ~(~0 << 28);

  // Check the event classification information for HadronB criteria.
  Evtcls_hadronic_flag_Manager &hadronFlagManager
      = Evtcls_hadronic_flag_Manager::get_manager();
  Evtcls_hadronic_flag &hadronFlags( hadronFlagManager( Panther_ID( 1 ) ) );
  float hadronBFlag = hadronFlags.hadronic_flag( 2 );

  // Print hadron flag to the log.
  if ( hadronBFlag < 10 ) {
    std::cout << "HADRONB CUT" << std::endl;
    // The following line prevents this event from being written to the ntuple
    //   This is not necessary and can be done at the Cern ROOT analysis
    //   level with the information contained in the ntuple.

    // continue;
  }
  
  // Write the value of the MC flag to the log for diagnostics.
  // std::cout << "(flagMC " << flagMC << ") ";

  // Get the Fox-wolfram R2 value for the event. Spherical events accepted (Low
  //   R2 values). If R2 is not calculated in the hadron info table, event will
  //   be accepted by default with an R2 of -1 to point out in the ntuple that
  //   the cut is "effectively off".
  double foxWolframR2 = -1;
  Evtcls_hadron_info_Manager& 
      hadMgr = Evtcls_hadron_info_Manager::get_manager();
  if ( hadMgr.count() ) {
    foxWolframR2 = hadMgr[0].R2();
  }

  // Write the event number to the log for diagnostics.
  // std::cout << experimentNumber << "," 
  //     << runNumber << ","
  //     << eventNumber << std::endl;
  
  // Instantiate constant data structures stored in external header files.
  PDGmasses masses;
  AdcabCuts cuts;
  
  // Add event to counter.
  ++numberOfEvents;
  
  // Set particle counters.
  int numCndtEPlus = 0;
  int numCndtEMinus = 0;
  int numCndtMuPlus = 0;
  int numCndtMuMinus = 0;

  // Set dilepton(++/--) counters.
  int numCndtEPlusPlus = 0;
  int numCndtMuPlusPlus = 0;
  int numCndtEMinusMinus = 0;
  int numCndtMuMinusMinus = 0;

  // Set correctly reconstructed MC particle counters.
  int numCrctEPlus = 0;
  int numCrctEMinus = 0;
  int numCrctMuPlus = 0;
  int numCrctMuMinus = 0;
  
  int numCrctEPlusPlus = 0;
  int numCrctMuPlusPlus = 0;
  int numCrctEMinusMinus = 0;
  int numCrctMuMinusMinus = 0;
    
  // Define lists (vector template) to store event particles.
  // Need a list for all mother and daughter particle species.
  // Reserve a large amount of memory ( 50 elements ) to avoid 
  //   relocating list to a new block of memory in case of large
  //   events.
  static std::vector< Particle > initialElectronList( 50 );
  static std::vector< Particle > electronList( 50 );
  static std::vector< Particle > initialMuonList( 50 );
  static std::vector< Particle > muonList( 50 );
  static std::vector< Particle > leptonList( 100 );
  static std::vector< Particle > kaonList( 50 );

  // Define list to store candidate dilepton events.
  static std::vector< DileptonEvent > dileptonEventList( 50 );
  
  // Ensure lists are empty so as not to consume too much memory.
  initialElectronList.clear();
  electronList.clear();
  initialMuonList.clear();
  muonList.clear();
  leptonList.clear();
  kaonList.clear();
  dileptonEventList.clear();
  
  // Alias the MDST charged manager, which contains
  // the measured charged tracks for each event.
  Mdst_charged_Manager &chg_mgr = Mdst_charged_Manager::get_manager();
  
  // Create an alias "null" for the null value of generated signal MC.
  // I.E. if something is not signal MC, it will be ascribed this value.
  const Gen_hepevt &null = Gen_hepevt_Manager::get_manager().get_NULL();
  
  // Populate the lepton candidate lists.
  for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin();
      i != chg_mgr.end(); ++i ) {

    // Alias the current particle as "chg".
    const Mdst_charged &chg = *i;
    
    // Get electron and muon liklihoods.
    eid chgEid( chg );
    double eidProb = chgEid.prob( 3, -1, 5 );
    Muid_mdst chgMuid( chg );
    double muidProb = chgMuid.Muon_likelihood();

    // Reject particle if below both electron and muon liklihood cuts.
    if ( eidProb < cuts.minEidProb && muidProb < cuts.minMuidProb ) {
      continue;
    }
    
    // Cut on IP dr and dz.
    // This is to make sure that particles were created near the IP.
    // TODO - 2010.08.11 - Is mass hypothesis = 3 (kaon) appropriate? Check with
    //                       authorities!
    IpDrDz ipDrDzParameters( chg, ip, 3 );
    if ( abs( ipDrDzParameters.dr() ) > cuts.maxIpDr ||
         abs( ipDrDzParameters.dz() ) > cuts.maxIpDz ) {
      continue;
    }
      
    // Mdst_trk member function mhyp( int hypID ) returns the fitted
    //   track parameters assuming certain particle mass hypotheses set
    //   by hypID. Possible mass hypotheses are:
    //   hypID = 0:e; 1:mu; 2:pi; 3:K; 4:p.

    // Mdst_trk_fit member function nhits( int detID) returns the
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
    Mdst_trk_fit &tfit = chg.trk().mhyp( 1 );
    if ( tfit.nhits( 3 ) < cuts.minSvdRHits ) continue;
    if ( tfit.nhits( 4 ) < cuts.minSvdZHits ) continue;

    // Reject particle if likelihoods are the same from eid and muid.
    if ( eidProb == muidProb ) continue;

    // Assume the particle species is that of the highest PID liklihood.
    if ( eidProb > muidProb ) {
      // Cut on minimum EID.
      if ( eidProb < cuts.minEidProb ) continue;
      
      // Assuming particle is an electron.
      Particle eCandidate( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );

      // Reject if the particle track has polar angle pointing outside the
      //   barrel (p is given closest to coord. origin - see mdst table).
      if ( eCandidate.p().cosTheta() < cuts.minLeptonCosTheta ||
           eCandidate.p().cosTheta() > cuts.maxLeptonCosTheta ) {
        continue;
      }

      // Cut on lepton CM momentum.
      // Get the lab-frame momentum and boost it to CM frame.
      HepLorentzVector candidatePCm( eCandidate.p() );
      candidatePCm.boost( cmBoostVector );
      double candidatePCmMag = candidatePCm.vect().mag();
      if ( candidatePCmMag < cuts.minLeptonMomentumCm ||
           candidatePCmMag > cuts.maxLeptonMomentumCm ) {
        continue;
      }

      // Store the MC truth to the lepton candidate.
      setMCtruth( eCandidate );

      // Add eCandidate to list of e+/- candidates.
      initialElectronList.push_back( eCandidate );

      // Collect diagnostic information.
      if ( eCandidate.relation().genHepevt() ) {
        double motherId = 0;
        if ( eCandidate.relation().genHepevt().mother() ) {
          motherId = eCandidate.relation().genHepevt().mother().idhep();
        }
        if ( motherId == -531 ) {
          num_bs_after_lepton_level.first++;
        }
        if ( motherId == 531 ) {
          num_bs_after_lepton_level.second++;
        }
      }
    } else {
      // Cut on minimum MUID.
      if ( muidProb < cuts.minMuidProb ) continue;
      
      // Assuming particle is a muon.
      Particle muCandidate( chg,
          chg.charge() > 0 ? ptypeMuPlus : ptypeMuMinus );

      // Reject if the particle track has polar angle pointing outside the
      //   barrel (p is given closest to coord. origin - see mdst table).
      if ( muCandidate.p().cosTheta() < cuts.minLeptonCosTheta ||
           muCandidate.p().cosTheta() > cuts.maxLeptonCosTheta ) {
        continue;
      }

      // Cut on lepton CM momentum.
      // Get the lab-frame momentum and boost it to CM frame.
      HepLorentzVector candidatePCm( muCandidate.p() );
      candidatePCm.boost( cmBoostVector );
      double candidatePCmMag = candidatePCm.vect().mag();
      if ( candidatePCmMag < cuts.minLeptonMomentumCm ||
           candidatePCmMag > cuts.maxLeptonMomentumCm ) {
        continue;
      }

      // Store the MC truth to the lepton candidate.
      setMCtruth( muCandidate );
      
      // Collect diagnostic information.
      if ( muCandidate.relation().genHepevt() ) {
        double motherId = 0;
        if ( muCandidate.relation().genHepevt().mother() ) {
          motherId = muCandidate.relation().genHepevt().mother().idhep();
        }
        if ( motherId == -531 ) {
          num_bs_after_lepton_level.first++;
        }
        if ( motherId == 531 ) {
          num_bs_after_lepton_level.second++;
        }
      }
      
      // Add eCandidate to list of e+/- candidates.
      muonList.push_back( muCandidate );
    }
  } // End for() loop populating lepton lists.

  // Write diagnostic information to the log.
  if ( flagVERBOSELOG ) {
  std::cout << num_bs_after_lepton_level.first << " "
      << num_bs_after_lepton_level.second << std::endl;
  }

  // Report number of electons in electron list for diagnostic purposes.
  // std::cout << "Number of electrons in event: "
  //           << initialElectronList.size() << std::endl;

  // Remove possible pair production electrons and J/Psi candidates.
  for ( std::vector< Particle >::iterator j = initialElectronList.begin();
      j != initialElectronList.end(); ++j ) {
    const Particle &eCndt = *j;
    
    // Flag for electrons that are not candidates for pair production daughters.
    // By default, assume all electrons are good.
    bool flagGoodElectron = true;

    // If the invariant mass of an electron candidate and every other opposite
    //   charged tracks is smaller than the cut value (nominally 100 MeV), the
    //   electron candidate is rejected.
    for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin();
        i != chg_mgr.end(); ++i ) {
      const Mdst_charged &chg = *i;
      Particle otherChg( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );
      
      // We need to know if the pair is same sign or opposite sign.
      bool ss_pair = ( eCndt.charge() == otherChg.charge() );

      // Calculate the invariant mass of the electron candidate and the other
      //   charged track.
      HepLorentzVector eCndtP = eCndt.p();
      HepLorentzVector otherChgP = otherChg.p();
      double electronChargedMass = ( eCndtP + otherChgP ).m();      
      
      // Report invariant mass for diagnostic purposes.
      // std::cout << "e p/m mass: " << electronChargedMass << std::endl;

      // If at any time eCndt proves to be likely from pair production,
      //   or a J/Psi, the flag is switched.
      if ( !( ss_pair ) && ( electronChargedMass < cuts.minEPlusEMinusMass ) ) {
        flagGoodElectron = false;
      }
      
      double deltaMass = electronChargedMass - cuts.massJPsi;
      // Same sign charge pairs may or may not be included in J/Psi veto
      //   depending on flag.
      if ( flag_jpsi_veto_os_only && ss_pair ) {
        // Nothing to do in this case.
        continue;
      } else if ( cuts.minElElJPsiCandidate < deltaMass ||
          deltaMass < cuts.maxElElJPsiCandidate ) {
        flagGoodElectron = false;
      }
    }

    // While eCndt is still a good candidate, add it to the electron list.
    if ( flagGoodElectron )
      electronList.push_back( eCndt );
  }
  
  // Remove possible J/Psi daughter muons.
  for ( std::vector< Particle >::iterator j = initialMuonList.begin();
      j != initialMuonList.end(); ++j ) {
    const Particle &muCndt = *j;
    
    // Flag for muons that are not candidates J/psi daughters.
    // By default, assume all muons are good.
    bool flagGoodMuon = true;

    // If the difference of the invariant mass of a muon candidate and every
    //   other opposite charged tracks is within the cut range, the candidate
    //   is rejected.
    for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin();
        i != chg_mgr.end(); ++i ) {
      const Mdst_charged &chg = *i;
      Particle otherChg( chg, chg.charge() > 0 ? ptypeMuPlus : ptypeMuMinus );

      // J/Psi veto only same sign charge muon pairs, depending on parameter.
      if ( flag_jpsi_veto_os_only &&
          ( muCndt.charge() == otherChg.charge() ) ) {
        continue;
      }
      
      // Calculate the invariant mass of the muon candidate and the other
      //   charged track.
      HepLorentzVector muCndtP = muCndt.p();
      HepLorentzVector otherChgP = otherChg.p();
      double muonChargedMass = ( muCndtP + otherChgP ).m();      
      
      // Report invariant mass for diagnostic purposes.
      // std::cout << "e p/m mass: " << muonChargedMass << std::endl;

      // If at any time muCndt proves to be likely from a J/Psi, the flag
      //   is switched.
      double deltaMass = muonChargedMass - cuts.massJPsi;
      if ( cuts.minMuMuJPsiCandidate < deltaMass ||
          deltaMass < cuts.maxMuMuJPsiCandidate ) {
        flagGoodMuon = false;
      }
    }

    // While muCndt is still a good candidate, add it to the electron list.
    if ( flagGoodMuon )
      muonList.push_back( muCndt );
  }
  
  // Report size of cleaned electron list for diagnotic purposes.
  // std::cout << "Number of electrons in event: " << electronList.size() 
  //     << std::endl;
  
  // Combine all leptons into a single list.
  leptonList.reserve( electronList.size() + muonList.size() );
  leptonList.insert( leptonList.end(), electronList.begin(),
                     electronList.end() );
  leptonList.insert( leptonList.end(), muonList.begin(), muonList.end() );
 
  
  // If wanting a verbose log, loop over the lepton list and write diagnostic
  //   selection information to the log.
  if ( flagVERBOSELOG ) {
    for ( std::vector< Particle >::iterator j = leptonList.begin();
        j != leptonList.end(); ++j ) {
      Particle &lepton0 = *j;
      if ( lepton0.relation().genHepevt() ) {
        double motherId = 0;
        if ( lepton0.relation().genHepevt().mother() ) {
          motherId = lepton0.relation().genHepevt().mother().idhep();
        }
        if ( motherId == -531 ) {
          num_bs_after_pair_removal.first++;
        }
        if ( motherId == 531 ) {
          num_bs_after_pair_removal.second++;
        }
      }
    }

  // Write diagnostic information to the log.
  std::cout << "           " 
      << num_bs_after_pair_removal.first
      << " "
      << num_bs_after_pair_removal.second << std::endl;
  }

  // Find dilepton event candidates.
  // Loop over the lepton list.
  for ( std::vector< Particle >::iterator j = leptonList.begin();
      j != leptonList.end(); ++j ) {
    
    // Loop over remaining leptons in the list to check all lepton pairs.
    for ( std::vector< Particle >::iterator i = j;
        i != leptonList.end(); ++i ) {
      
      // Exclude the case where both iterators point to the same particle.
      if ( i == j ) continue;
      
      // Alias the first lepton as "lepton0"
      //  and the second electron as "lepton1".
      Particle &lepton0 = *j;
      Particle &lepton1 = *i;
      
      // Cut on each lepton's CM momentum.
      // TODO - 2010.08.10 - Perhaps this should be done when populating the
      //                     lepton lists.
      // Get the lab-frame momentum for each lepton and boost it to CM frame.
      HepLorentzVector lepton0MomentumCm( lepton0.p() );
      HepLorentzVector lepton1MomentumCm( lepton1.p() );
      lepton0MomentumCm.boost( cmBoostVector );
      lepton1MomentumCm.boost( cmBoostVector );

      // Make the CM momentum cut.
      double lepton0MomentumCmMag = lepton0MomentumCm.vect().mag();
      double lepton1MomentumCmMag = lepton1MomentumCm.vect().mag();
      if ( lepton0MomentumCmMag < cuts.minLeptonMomentumCm ||
          lepton0MomentumCmMag > cuts.maxLeptonMomentumCm ) {
        continue;
      }

      if ( lepton1MomentumCmMag < cuts.minLeptonMomentumCm ||
          lepton1MomentumCmMag > cuts.maxLeptonMomentumCm ) {
        continue;
      }
      
      // Determine higher momentum lepton and add it to an event candidate.
      //   lepton0 is henceforth considered the higher momentum lepton.
      Particle &leptonHigherP = ( lepton0MomentumCmMag > lepton1MomentumCmMag ?
          lepton0 : lepton1 );
      Particle &leptonLowerP = ( lepton0MomentumCmMag > lepton1MomentumCmMag ?
          lepton1 : lepton0 );
      DileptonEvent eventCandidate( leptonHigherP, leptonLowerP,
          cmBoostVector );
      
      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosThetaLLCm = eventCandidate.cosThetaLL();
      if ( cosThetaLLCm < cuts.minCosThetaLLCm ||
          cosThetaLLCm > cuts.maxCosThetaLLCm ) {
        continue;
      }

      // Add event candidate to the list of dilepton events.
      dileptonEventList.push_back( eventCandidate );

      if ( flagVERBOSELOG ) {
        double lepton0MotherId = eventCandidate.l0().idMother();
        double lepton1MotherId = eventCandidate.l1().idMother();
        if ( lepton0MotherId == -531 ) num_bs_at_after_event_selection.first++;
        if ( lepton0MotherId ==  531 ) num_bs_at_after_event_selection.second++;
        if ( lepton1MotherId == -531 ) num_bs_at_after_event_selection.first++;
        if ( lepton1MotherId ==  531 ) num_bs_at_after_event_selection.second++;
      }
    }
  }
  
  // Write diagnostic information to the log.
  if ( flagVERBOSELOG ) {
    std::cout << "                     " 
        << num_bs_at_after_event_selection.first
        << " "
        << num_bs_at_after_event_selection.second << std::endl;
  }

  // Report the size of the number of dilepton events for diagnostic purposes.
  // std::cout << "numDileptonEvents = " << numDileptonEvents << "\n"
  //     << "   numberOfEvents = " << numberOfEvents << std::endl;

  // Choose single best event candidate if more than one.
  if ( flagSELECTBESTCANDIDATE ) {
    DileptonEvent bestEvent;
    double bestEventMomentum;
    for ( std::vector< DileptonEvent >::iterator i = dileptonEventList.begin();
        i != dileptonEventList.end(); ++i ) {
      DileptonEvent &currentEvent = *i;
      if ( i == dileptonEventList.begin() ) {
        bestEventMomentum = currentEvent.pSum();
        bestEvent = currentEvent;
      }
      double currentEventMomentum = currentEvent.pSum();
      if ( currentEventMomentum > bestEventMomentum ) {
        bestEventMomentum = currentEventMomentum;
        bestEvent = currentEvent;
      }
    }
    // If the list of dilepton events is empty, better not
    //   treat the past-end-of-vector object as an event!
    if ( dileptonEventList.begin() != dileptonEventList.end() ) {
      dileptonEventList.clear();
      dileptonEventList.push_back( bestEvent );
    }
  }

  // Write the best candidate information to n-tuple.
  for ( std::vector< DileptonEvent >::iterator i = dileptonEventList.begin();
      i != dileptonEventList.end(); ++i ) {
    DileptonEvent &eventCandidate = *i;

    // Rebuild the dr and dz track parameters for each lepton.
    int l0MHyp = ( abs( eventCandidate.l0().idAssigned() ) == 13 ) ? 1 : 0;
    int l1MHyp = ( abs( eventCandidate.l1().idAssigned() ) == 13 ) ? 1 : 0;
    IpDrDz l0IpDrDz( eventCandidate.l0().particle().relation().mdstCharged(),
        ip, l0MHyp );
    IpDrDz l1IpDrDz( eventCandidate.l1().particle().relation().mdstCharged(),
        ip, l1MHyp );

    // Write data to the n-tuple. As of 2010.08.25, each row in the n-tuple will
    //   consist of a single dilepton event candidate.
    // Column names can be no greater than eight (8) characters long.
    // Format:       "name    ", value )
    // Event singular information
    nTuple_->column( "exp_no"  , experimentNumber );
    nTuple_->column( "run_no"  , runNumber );
    nTuple_->column( "evt_no"  , eventNumber );
    nTuple_->column( "fw_r2"   , foxWolframR2 );
    nTuple_->column( "hadronb" , hadronBFlag );
    
    // Event candidate information.
    nTuple_->column( "evnttype", eventCandidate.eventType() );
    nTuple_->column( "llcostha", eventCandidate.cosThetaLL() );
    nTuple_->column( "pcm_sum" , eventCandidate.pSum() );
    nTuple_->column( "pcm_dif" , eventCandidate.pDifference() );

    // Lepton information.
    nTuple_->column( "l0_chg"  , eventCandidate.l0().particle().charge() );
    nTuple_->column( "l1_chg"  , eventCandidate.l1().particle().charge() );
    nTuple_->column( "l0_m"    , eventCandidate.l0().p().mag() );
    nTuple_->column( "l1_m"    , eventCandidate.l1().p().mag() );
    nTuple_->column( "l0_idasn", eventCandidate.l0().idAssigned() );
    nTuple_->column( "l1_idasn", eventCandidate.l1().idAssigned() );
    nTuple_->column( "l0_idtru", eventCandidate.l0().idTruth() );
    nTuple_->column( "l1_idtru", eventCandidate.l1().idTruth() );
    nTuple_->column( "l0_idmom", eventCandidate.l0().idMother() );
    nTuple_->column( "l1_idmom", eventCandidate.l1().idMother() );
    nTuple_->column( "l0_eidp" , eventCandidate.l0().likelihoodE() );
    nTuple_->column( "l1_eidp" , eventCandidate.l1().likelihoodE() );
    nTuple_->column( "l0_muidp", eventCandidate.l0().likelihoodMu() );
    nTuple_->column( "l1_muidp", eventCandidate.l1().likelihoodMu() );
    nTuple_->column( "l0_muidr", eventCandidate.l0().klmHitsChi2PerN() );
    nTuple_->column( "l1_muidr", eventCandidate.l1().klmHitsChi2PerN() );
    nTuple_->column( "l0_pcm"  , eventCandidate.l0().pCm().vect().mag() );
    nTuple_->column( "l1_pcm"  , eventCandidate.l1().pCm().vect().mag() );
    nTuple_->column( "l0_plab" , eventCandidate.l0().p().vect().mag() );
    nTuple_->column( "l1_plab" , eventCandidate.l1().p().vect().mag() );
    nTuple_->column( "l0_cme"  , eventCandidate.l0().pCm().e() );
    nTuple_->column( "l0_cme"  , eventCandidate.l0().pCm().e() );
    nTuple_->column( "l0_cmpx" , eventCandidate.l0().pCm().px() );
    nTuple_->column( "l0_cmpx" , eventCandidate.l0().pCm().px() );
    nTuple_->column( "l0_cmpy" , eventCandidate.l0().pCm().py() );
    nTuple_->column( "l0_cmpy" , eventCandidate.l0().pCm().py() );
    nTuple_->column( "l0_cmpz" , eventCandidate.l0().pCm().pz() );
    nTuple_->column( "l0_cmpz" , eventCandidate.l0().pCm().pz() );
    nTuple_->column( "l0_costh", eventCandidate.l0().p().cosTheta() );
    nTuple_->column( "l1_costh", eventCandidate.l1().p().cosTheta() );
    nTuple_->column( "l0_dr"   , l0IpDrDz.dr() );
    nTuple_->column( "l1_dr"   , l1IpDrDz.dr() );
    nTuple_->column( "l0_dz"   , l0IpDrDz.dz() );
    nTuple_->column( "l1_dz"   , l1IpDrDz.dz() );
    nTuple_->column( "l0_svdr" , eventCandidate.l0().svdHitsR() );
    nTuple_->column( "l1_svdr" , eventCandidate.l1().svdHitsR() );
    nTuple_->column( "l0_svdz" , eventCandidate.l0().svdHitsZ() );
    nTuple_->column( "l1_svdz" , eventCandidate.l1().svdHitsZ() );
    nTuple_->dumpData();
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

  nTuple_ = tm->ntuple( "Dilepton", ntList, 1 );
  
  return;
}

#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
