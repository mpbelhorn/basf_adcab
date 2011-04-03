//
//******************************************************************************
// Filename: Adcab.cc
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010-07-19
// Description: Analysis of the (A)nomalous (D)ilepton (C)harge (A)symmetry in
//   (B)s0 decays.
//******************************************************************************


//******************************************************************************
//  Preamble
//******************************************************************************

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

#include <panther/panther.h>    // Panther.
#include BELLETDF_H             // Panther.
#include HEPEVT_H               // Panther.
#include MDST_H                 // Panther.
#include EVTCLS_H               // Panther - Event Classification info.

#include "HEPconstants.h"       // PDG masses and constants.
#include "AnalysisTools.h"      // General analysis functions and utilities.
#include "AdcabCuts.h"          // Analysis specfic selection cut constants.
#include "geninfo.h"            // Custom analysis functions to use Zupanc's MC.
#include "userinfo.h"           // Custom analysis functions to use Zupanc's MC.

#if defined(BELLE_NAMESPACE)    // Namespace container for backwards
namespace Belle {               //  compatibility with older versions of
#endif                          //  BELLE Library (used for b200611xx onward).
                                //  Must be in all files.


//******************************************************************************
// Adcab Class Definition
//******************************************************************************

// Declare analysis class, inheriting Module from BASF.
class Adcab : public Module {
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
  int eventsInHeaderNumber;
  
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
  bool flagMC;         // Data type flag. 'true' = MC, 'false' = Real Data.
  bool flagERROR;      // Error flag. 'true' = fatal error or unrelaiable data.
  bool flagVERBOSELOG; // Turns on verbose diagnostic messages in the log.
  bool flagLOGBSASSYM; // Writes Bs - BsBar charge asym diagnostics to the log.

private:
  BelleTuple *nTuple_;  // Pointer for writing to the n-tuple. 
};


// Registers analysis module in BASF.
extern "C" Module_descr *mdcl_Adcab() 
{ 
  // Creates pointer to allocated Adcab class object.
  Adcab *module = new Adcab;
  
  // Creates pointer "dscr" to description of Adcab module.
  Module_descr *dscr = new Module_descr( "Adcab", module );
  
  // Provide path to pass paramaters to BeamEnergy class.
  BeamEnergy::define_global( dscr );
  
  return dscr;
}


// Adcab constructor definition.
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

  // Inform of succesful module load.
  std::cout 
    << "\n\n"
    << "**************************************************\n"
    << "* Adcab Analysis Module loaded successfully *\n"
    << "**************************************************\n"
    << std::endl;
  
  return;
}


//******************************************************************************
// Module Initialization Functions
//******************************************************************************

// init() definition.
void Adcab::init(int *)
{
  return;
}

// term() definition.
void Adcab::term()
{
  std::cout
    << "\n\n"
    << "******************************************************\n"
    << "* Adcab Analysis Module terminated successfully *\n"
    << "******************************************************\n"
    << std::endl;

  return;
}


//******************************************************************************
// Event Analysis Functions
//******************************************************************************

// This function is exectuted once per data run. Thus it resets the class flags
//   and event run counters, it initializes the runhead information (data type,
//   experiment number, run number, mdst access), it retrieves the run-dependent
//   beam and interaction point information. A basic run summery is dumped to
//   the log for the record.
void Adcab::begin_run(BelleEvent* evptr, int *status) 
{
  (void)evptr;
  (void)status;
  
  // Set run information to default values.
  experimentNumber  = 0;
  runNumber         = 0;
  numberOfEvents    = 0;
  numDileptonEvents = 0;
  eventsInHeaderNumber = 0;
  
  // Set Diagnostic Variables to 0.
  num_bs_after_lepton_level.first = 0;
  num_bs_after_lepton_level.second = 0;
  num_bs_after_pair_removal.first = 0;
  num_bs_after_pair_removal.second = 0;
  num_bs_at_after_event_selection.first = 0;
  num_bs_at_after_event_selection.second = 0;

  // Set default flags.
  flagMC    = false;
  flagERROR = false;
  flagVERBOSELOG = false;
  
  // Set interaction point and error to default values.
  ip       = HepPoint3D( 0, 0, 0 );
  ipErr    = HepSymMatrix( 3, 0 );
  ipUsable = 0;

  // cmFourMomentum = HepLorentzVector( 0, 0, 0, 0 );
  
  // Get BELLE runhead manager.
  Belle_runhead_Manager &rhd_mgr = Belle_runhead_Manager::get_manager();
  std::vector< Belle_runhead >::const_iterator rhd = rhd_mgr.begin();
  
  // Check access to BELLE runhead.
  if ( rhd == rhd_mgr.end() ) {
    // If BELLE runhead not available, flag the error, inform the log
    //   and return.
    fprintf( stderr, "Cannot access Belle_runhead\n" );
    flagERROR = true;
    
    return;
  } else {
    // Set runhead information tags.
    experimentNumber = rhd_mgr[0].ExpNo();
    runNumber = rhd_mgr[0].RunNo();
    eventsInHeaderNumber = rhd_mgr[0].NEvt();
    
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
      << "*******************\n"
      << "* Run Information *\n"
      << "*******************\n"
      << std::endl;

    if ( rhd->ExpMC() == 1 ) {
      flagMC = false; // Set Data type flag to Real Data.
      std::cout << "Data is Real." << std::endl;
    } else {
      flagMC = true;  // Set Data type flag to Monte Carlo.
      std::cout << "Data is Monte Carlo." << std::endl;
    }
    
    std::cout
      << "Experiment " << experimentNumber << ", Run " << runNumber 
      << ", Events " << eventsInHeaderNumber << "\n" 
      << "Actual Beam Energy: " << beamEnergyCMFrame 
      << " +/- " << beamEnergyError << " GeV\n"
      << "Reported Beam Energy: " << kekbBeamEnergy << " GeV\n"
      << "BE Class cmBoostVector: " << cmBoostVector << "\n"
      << std::endl;
  }
  
  return;
}

// This function is run once at the end of a data run. It does nothing but
//   write a message to somewhere. This message does not appear in log
//   although MPB believes it should.
void Adcab::end_run(BelleEvent* evptr, int *status )
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
void Adcab::event(BelleEvent* evptr, int* status)
{
  (void)evptr;
  (void)status;

  // Get the event number. 
  // The bitwise operation "& ~(0 << 28)" forces the event number to count from 
  //   zero again after EvtNo() reaches 2^28. Not sure why this is necessary?
  Belle_event_Manager& EvMgr = Belle_event_Manager::get_manager();
  eventNumber = EvMgr[0].EvtNo() & ~(~0 << 28);

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
  //           << runNumber << ","
  //           << eventNumber << std::endl;
  
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
  // Reserve a large amount of memory ( 100 elements ) to avoid 
  //   relocating list to a new block of memory in case of large
  //   events.
  static std::vector< Particle > uncleanElectronList( 100 );
  static std::vector< Particle > electronList( 100 );
  static std::vector< Particle > muonList( 100 );
  static std::vector< Particle > leptonList( 100 );
  static std::vector< Particle > kaonList( 100 );
  static std::vector< Particle > dsPlusList( 100 );
  static std::vector< Particle > dsMinusList( 100 );
  static std::vector< Particle > dsStarPlusList( 100 );
  static std::vector< Particle > dsStarMinusList( 100 );

  // Define list to store candidate dilepton events.
  static std::vector< DileptonEvent > dileptonEventList( 100 );
  
  // Ensure lists are empty so as not to consume too much memory.
  uncleanElectronList.clear();
  electronList.clear();
  muonList.clear();
  leptonList.clear();
  kaonList.clear();
  dsPlusList.clear();
  dsMinusList.clear();
  dsStarPlusList.clear();
  dsStarMinusList.clear();
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
    if ( eidProb < cuts.minEidProb && muidProb < cuts.minMuidProb )
      continue;
    
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

    // Reject particle if liklihoods are the same from eid and muid.
    if ( eidProb == muidProb ) continue;

    // Assume the particle species is that of the highest PID liklihood.
    if ( eidProb > muidProb ) {

      // Assuming particle is an electron.
      Particle eCndt( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );

      // Cut on track momentum.
      if ( eCndt.p().vect().mag() < cuts.minElectronMomentum ) continue;

      // Cut on minimum EID.
      if ( eidProb < cuts.minEidProb ) continue;
      
      // Store the MC truth to the lepton candidate.
      setMCtruth( eCndt );

      // Add eCandidate to list of e+/- candidates.
      uncleanElectronList.push_back( eCndt );

      // Collect diagnostic information.
      if ( eCndt.relation().genHepevt() ) {
        double motherId = 0;
        if ( eCndt.relation().genHepevt().mother() ) {
          motherId = eCndt.relation().genHepevt().mother().idhep();
        }
        if ( motherId == -531 ) {
          num_bs_after_lepton_level.first++;
        }
        if ( motherId == 531 ) {
          num_bs_after_lepton_level.second++;
        }
      }
    } else {
      // Assuming particle is a muon.
      Particle muCandidate( chg, 
                            chg.charge() > 0 ? ptypeMuPlus : ptypeMuMinus );

      // Cut on track momentum.
      if ( muCandidate.p().vect().mag() < cuts.minMuonMomentum ) continue;

      // Cut on minimum MUID.
      if ( muidProb < cuts.minMuidProb ) continue;
      
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
  //           << uncleanElectronList.size() << std::endl;

  // Remove possible pair production electrons.
  for ( std::vector< Particle >::iterator j = uncleanElectronList.begin();
        j != uncleanElectronList.end(); ++j ) {
    const Particle &eCndt = *j;
    
    // Flag for electrons that are not candidates for pair production daughters.
    // By default, assume all electrons are good.
    bool flagNonPairProduction = true;

    // If the invariant mass of an electron candidate and every other opposite
    //   charged tracks is smaller than the cut value (nominally 100 MeV), the
    //   electron candidate is rejected.
    for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin();
          i != chg_mgr.end(); ++i ) {
      const Mdst_charged &chg = *i;
      Particle otherChg( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );

      // Same sign charge pairs cannot be pair production electrons.
      if ( eCndt.charge() == otherChg.charge() ) continue;
      
      // Calculate the invariant mass of the electron candidate and the other
      //   charged track.
      HepLorentzVector eCndtP = eCndt.p();
      HepLorentzVector otherChgP = otherChg.p();
      double ePlusEMinusMass = ( eCndtP + otherChgP ).m();      
      
      // Report invariant mass for diagnostic purposes.
      // std::cout << "e p/m mass: " << ePlusEMinusMass << std::endl;

      // If at any time eCndt proves to be likely from pair production,
      //   the flag is switched.
      if ( ePlusEMinusMass < cuts.minEPlusEMinusMass )
        flagNonPairProduction = false;
    }

    // While eCndt is still a good candidate, add it to the electron list.
    if ( flagNonPairProduction )
      electronList.push_back( eCndt );
  } 
  
  // Report size of cleaned electron list for diagnotic purposes.
  // std::cout << "Number of electrons in event: " << electronList.size() 
  //           << std::endl;
  
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
                                  lepton0 :
                                  lepton1 );
      Particle &leptonLowerP = ( lepton0MomentumCmMag > lepton1MomentumCmMag ?
                                 lepton1 :
                                 lepton0 );
      DileptonEvent eventCandidate( leptonHigherP, leptonLowerP );
      
      // Cut on jet-like events where the included angle between the leptons
      //   in the CM frame is near 0 or Pi.
      double cosThetaLLCm = eventCandidate.cosThetaLLCm( cmBoostVector );
      if ( cosThetaLLCm < cuts.minCosThetaLLCm ||
           cosThetaLLCm > cuts.maxCosThetaLLCm ) {
        continue;
      }

      // TODO - 2010.07.27 - Hastings makes a cut requiring hits in the
      //                     barrel CsI. Is this needed?
      
      // Add event candidate to the list of dilepton events.
      dileptonEventList.push_back( eventCandidate );
      ++numDileptonEvents;
      if ( flagVERBOSELOG ) {
        double lepton0MotherId = eventCandidate.l0MotherId();
        double lepton1MotherId = eventCandidate.l1MotherId();
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
  //           << "   numberOfEvents = " << numberOfEvents << std::endl;

  // TODO - 2010.08.12 - Choose single best event candidate if more than one.
  //Loop over the dilepton event candidates and write information to n-tuple.
  for ( std::vector< DileptonEvent >::iterator i = dileptonEventList.begin();
        i != dileptonEventList.end(); ++i ) {
    DileptonEvent &eventCandidate = *i;

    // Unpack the event candidate.
    Particle lepton0 = eventCandidate.lepton0();
    Particle lepton1 = eventCandidate.lepton1();

    // Recall the event type, lepton species and lepton mother IDs.
    double eventType = eventCandidate.eventType();
    double lepton0Species = lepton0.pType().lund();
    double lepton1Species = lepton1.pType().lund();
    double lepton0Id = eventCandidate.l0Id();
    double lepton1Id = eventCandidate.l1Id();
    double lepton0MotherId = eventCandidate.l0MotherId();
    double lepton1MotherId = eventCandidate.l1MotherId();

    // If the data is MC, check that the event was identified correctly.
    int mcBsParents = -1;
    if ( flagMC ) {
      mcBsParents = eventCandidate.mcBsParents();
    }
    
    // Get the lab-frame momentum for each lepton and boost it to CM frame.
    HepLorentzVector lepton0MomentumCm( lepton0.p() );
    HepLorentzVector lepton1MomentumCm( lepton1.p() );
    lepton0MomentumCm.boost( cmBoostVector );
    lepton1MomentumCm.boost( cmBoostVector );

    // Recall the lepton muid and eid probability values.
    const Mdst_charged &lepton0MdstCharged = lepton0.relation().mdstCharged();
    const Mdst_charged &lepton1MdstCharged = lepton1.relation().mdstCharged();
    eid lepton0Eid( lepton0MdstCharged );
    eid lepton1Eid( lepton1MdstCharged );
    Muid_mdst lepton0Muid( lepton0MdstCharged );
    Muid_mdst lepton1Muid( lepton1MdstCharged );
    double lepton0EidProb = lepton0Eid.prob( 3, -1, 5 );
    double lepton1EidProb = lepton1Eid.prob( 3, -1, 5 );
    double lepton0MuidProb = lepton0Muid.Muon_likelihood();
    double lepton1MuidProb = lepton1Muid.Muon_likelihood();

    // Rebuild the dr and dz track parameters for each lepton.
    int lepton0MassHyp = 0;
    int lepton1MassHyp = 0;
    if ( abs( lepton0Species ) == 13 ) {
      lepton0MassHyp = 1;
    }
    if ( abs( lepton1Species ) == 13 ) {
      lepton1MassHyp = 1;
    }    
    IpDrDz l0IpDrDz( lepton0.relation().mdstCharged(), ip, lepton0MassHyp );
    IpDrDz l1IpDrDz( lepton1.relation().mdstCharged(), ip, lepton1MassHyp );
    double lepton0Dr = l0IpDrDz.dr();
    double lepton0Dz = l0IpDrDz.dz();
    double lepton1Dr = l1IpDrDz.dr();
    double lepton1Dz = l1IpDrDz.dz();

    // Rebuild the SVD hit numbers.
    Mdst_trk_fit &l0TrackFit = lepton0.relation().mdstCharged().trk().mhyp( 1 );
    double lepton0SvdRHits = l0TrackFit.nhits( 3 );
    double lepton0SvdZHits = l0TrackFit.nhits( 4 );
    Mdst_trk_fit &l1TrackFit = lepton1.relation().mdstCharged().trk().mhyp( 1 );
    double lepton1SvdRHits = l1TrackFit.nhits( 3 );
    double lepton1SvdZHits = l1TrackFit.nhits( 4 );

    // Recall the CM frame lepton opening angle.
    double cosThetaLLCm = eventCandidate.cosThetaLLCm( cmBoostVector );

    // Write data to the n-tuple. As of 2010.08.25, each row in the n-tuple will
    //   consist of a single dilepton event candidate.
    // Column names can be no greater than eight (8) characters long.
    // Format:       "name    ", value )
    nTuple_->column( "l0_spec" , lepton0Species );
    nTuple_->column( "l1_spec" , lepton1Species );
    nTuple_->column( "l0_momid", lepton0MotherId );
    nTuple_->column( "l1_momid", lepton1MotherId );
    nTuple_->column( "l0_id"   , lepton0Id );
    nTuple_->column( "l1_id"   , lepton1Id );
    nTuple_->column( "l0_eidp" , lepton0EidProb );
    nTuple_->column( "l1_eidp" , lepton1EidProb );
    nTuple_->column( "l0_muidp", lepton0MuidProb );
    nTuple_->column( "l1_muidp", lepton1MuidProb );
    nTuple_->column( "l0_p3"   , lepton0.p().vect().mag() );
    nTuple_->column( "l1_p3"   , lepton1.p().vect().mag() );
    nTuple_->column( "l0_pt"   , lepton0.p().vect().perp() );
    nTuple_->column( "l1_pt"   , lepton1.p().vect().perp() );
    nTuple_->column( "l0_p3cm" , lepton0MomentumCm.vect().mag() );
    nTuple_->column( "l1_p3cm" , lepton1MomentumCm.vect().mag() );
    nTuple_->column( "l0_m"    , lepton0.p().mag() );
    nTuple_->column( "l1_m"    , lepton1.p().mag() );
    nTuple_->column( "l0_chg"  , lepton0.charge() );
    nTuple_->column( "l1_chg"  , lepton1.charge() );
    nTuple_->column( "llcostha", cosThetaLLCm );
    nTuple_->column( "fw_r2"   , foxWolframR2 );
    nTuple_->column( "exp_no"  , experimentNumber );
    nTuple_->column( "run_no"  , runNumber );
    nTuple_->column( "evt_no"  , eventNumber );
    nTuple_->column( "mcbsprnt", mcBsParents );
    nTuple_->column( "evnttype", eventType );
    nTuple_->column( "l0_dr"   , lepton0Dr );
    nTuple_->column( "l0_dz"   , lepton0Dz );
    nTuple_->column( "l1_dr"   , lepton1Dr );
    nTuple_->column( "l1_dz"   , lepton1Dz );
    nTuple_->column( "l0_svdr" , lepton0SvdRHits );
    nTuple_->column( "l0_svdz" , lepton0SvdZHits );
    nTuple_->column( "l1_svdr" , lepton1SvdRHits );
    nTuple_->column( "l1_svdz" , lepton1SvdZHits );
    nTuple_->dumpData();
  }
  
  return;
}


// Specifies n-tuples to write and names for root histograms.
void Adcab::hist_def()
{
  extern BelleTupleManager *BASF_Histogram;   // Define a BASF Histogram
  BelleTupleManager *tm = BASF_Histogram;   
  const char *ntList = "l0_spec "
                       "l1_spec "
                       "l0_momid "
                       "l1_momid "
                       "l0_id "
                       "l1_id "
                       "l0_eidp "
                       "l1_eidp "
                       "l0_muidp "
                       "l1_muidp "
                       "l0_p3 "
                       "l1_p3 "
                       "l0_pt "
                       "l1_pt "
                       "l0_p3cm "
                       "l1_p3cm "
                       "l0_m "
                       "l1_m "
                       "l0_chg "
                       "l1_chg "
                       "llcostha "
                       "fw_r2 "
                       "exp_no "
                       "run_no "
                       "evt_no "
                       "mcbsprnt "
                       "evnttype "
                       "l0_dr "
                       "l0_dz "
                       "l1_dr "
                       "l1_dz "
                       "l0_svdr "
                       "l0_svdz "
                       "l1_svdr "
                       "l1_svdz";

  nTuple_ = tm->ntuple( "Dilepton", ntList, 1 );
  
  return;
}


#if defined(BELLE_NAMESPACE)  // Needed to close compiler namespace
} // namespace Belle          // container for backwards compatibility
#endif                        // with older BELLE Libraries.
