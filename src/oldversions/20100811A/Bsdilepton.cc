//
//*******************************************************************************
// Filename: Bsdilepton.cc
// Version: 2010.08.11.A
// Author: M.P. Belhorn
// Original Date: 2010-07-19
// Description: Analysis of the dilepton charge asymmetry in Bs0 decays.
//*******************************************************************************


//*******************************************************************************
//	Preamble
//*******************************************************************************

#include <cmath>						// Uses cmath functions.

#include "belle.h"						// BELLE Library.
#include "event/BelleEvent.h"			// For managing BELLE events.
#include "tuple/BelleTupleManager.h"	// For managing BELLE nTuples.
#include "benergy/BeamEnergy.h"			// For determining run beam energy.
#include "basf/module.h"				// For interfacing with BASF.
#include "basf/module_descr.h"			// ???
#include "particle/Particle.h"			// The BELLE Particle Class.
#include "particle/utility.h"			// Additional functions for use with the Particle class.
#include "kid/atc_pid.h"				// For particle species separation.
#include "mdst/Muid_mdst.h"				// For muon identification.
#include "eid/eid.h"					// For electron identification.
#include "mdst/mdst.h"					// For MDST files.
#include "ip/IpProfile.h"				// Beam Interaction Point (IP) analysis tools. Position unit = cm.

#include <panther/panther.h>			// Panther.
#include BELLETDF_H						// Panther.
#include HEPEVT_H						// Panther.
#include MDST_H							// Panther.

#include "HEPconstants.h"				// PDG masses and constants.
#include "AnalysisTools.h"				// General analysis functions and utilities.
#include "BsdileptonCuts.h"				// Analysis specfic selection cut constants.

#if defined(BELLE_NAMESPACE)			// Namespace container for backwards
namespace Belle {						//	compatibility with older versions
#endif									//	of BELLE Library (used for b200611xx onward).
										//	Must be in all files.


//*******************************************************************************
// Bsdilepton Class Definition
//*******************************************************************************

// Declare analysis class, inheriting Module from BASF.
class Bsdilepton : public Module {

public:
	// Class functions.
	Bsdilepton();								// Constructor.
	~Bsdilepton() {}							// Destructor. Not used.
	
	// BASF Module functions.
	void init( int * );							// Executed once when module is loaded in BASF.
	void disp_stat( const char* ) {}			// Not used.
	void hist_def();							// Executed on BASF "define histogram" call.
	void begin_run( BelleEvent*, int* );		// Executed at beginning of each run.
	void event( BelleEvent*, int* );			// Exectued for each event in run.
	void end_run( BelleEvent*, int* );			// Executed at end of each run.
	void other( int*, BelleEvent*, int* ) {}	// Not used.
	void term();								// Executed once when module is terminated in BASF.
	
	// Custom analysis class functions.
	void setBeamInformation();					// Calculates event beam energy.
	
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
	double kekbBeamEnergy;		// Uncalibrated energy reported by KEKB.
	double kekbLerBeamEnergy;	// Uncalibrated energy reported by KEKB.
	double kekbHerBeamEnergy;	// Uncalibrated energy reported by KEKB.
	double beamCrossingAngle;
				
private:
	BelleTuple *m_nt;	// Define pointer for writing to the n-tuple.
	
	// Flags
	bool flagMC;		// Data type flag. Values indicate 'true' = MC, 'false' = Real Data.
	bool flagERROR;		// Error flag. Value 'true' indicates fatal errors or unrelaiable data.
	
};


// Registers analysis module in BASF.
extern "C" Module_descr *mdcl_Bsdilepton() 
{	
	// Creates pointer to allocated Bsdilepton class object.
	Bsdilepton *module = new Bsdilepton;
	
	// Creates pointer "dscr" to description of Bsdilepton module.
	Module_descr *dscr = new Module_descr( "Bsdilepton", module );
	
	// Provide path to pass paramaters to BeamEnergy class.
	BeamEnergy::define_global( dscr );
	
	return dscr;
}


// Bsdilepton constructor definition.
Bsdilepton::Bsdilepton()
{	
	// Define particle types (Ptype) constants for particle class objects.
	// Valid names are those in the qq98 decay file located at
	// $BELLE_TOP_DIR/share/data-files/qq98/decay.dec
	ptypeEMinus			= ( Ptype( "E-" ) );
	ptypeEPlus			= ( Ptype( "E+" ) );
	ptypeNuE			= ( Ptype( "NUE" ) );
	ptypeNuEBar			= ( Ptype( "NUEB" ) );
	ptypeMuMinus		= ( Ptype( "MU-" ) );
	ptypeMuPlus			= ( Ptype( "MU+" ) );
	ptypeNuMu			= ( Ptype( "NUM" ) );
	ptypeNuMuBar		= ( Ptype( "NUMB" ) );
	ptypeKMinus			= ( Ptype( "K-" ) );
	ptypeKPlus			= ( Ptype( "K+" ) );
	ptypeDsMinus		= ( Ptype( "DS-" ) );
	ptypeDsPlus			= ( Ptype( "DS+" ) );
	ptypeDsStarMinus	= ( Ptype( "DS*-" ) );
	ptypeDsStarPlus		= ( Ptype( "DS*+" ) );

	// Inform of succesful module load.
	std::cout << "\n********* Bsdilepton Analysis Module loaded successfully *********\n" << std::endl;
	
	return;
}


//*******************************************************************************
// Module Initialization Functions
//*******************************************************************************

// init() definition.
void Bsdilepton::init(int *)
{
	return;
}

// term() definition.
void Bsdilepton::term()
{
	std::cout << "\n********* Bsdilepton Analysis Module terminated successfully *********\n"
			  << "Analysis summary information follows...\n" << std::endl;

	return;
}


//*******************************************************************************
// Event Analysis Functions
//*******************************************************************************

// Sets data flags and initializes some analysis functions.
void Bsdilepton::begin_run(BelleEvent* evptr, int *status) 
{
	(void)evptr;
	(void)status;
	
	// Set run information to default values.
	experimentNumber	= 0;
	runNumber 			= 0;
	numberOfEvents		= 0;
	numDileptonEvents	= 0;
	
	// Set default flags.
	flagMC				= false;
	flagERROR			= false;
	
	// Set interaction point and error to default values.
	ip					= HepPoint3D( 0, 0, 0 );
	ipErr				= HepSymMatrix( 3, 0 );
	ipUsable			= 0;

	// cmFourMomentum	= HepLorentzVector( 0, 0, 0, 0 );
	
	// Get BELLE runhead manager.
	Belle_runhead_Manager &rhd_mgr = Belle_runhead_Manager::get_manager();
	std::vector< Belle_runhead >::const_iterator rhd = rhd_mgr.begin();
	
	// Check access to BELLE runhead.
	if ( rhd == rhd_mgr.end() ) {
		fprintf( stderr, "Cannot access Belle_runhead\n" );
		flagERROR = true;	// If unable to read BELLE runhead, set flagERROR to true setting to indicate error.
		
		return;
	} 
	// If no error reading BELLE runhead, set run information and determine data type.
	else {
		// Set runhead information tags.
		experimentNumber = rhd_mgr[0].ExpNo();
		runNumber = rhd_mgr[0].RunNo();
		
		// Initialise BeamEnergy class.
		BeamEnergy::begin_run();
		
		beamEnergyCMFrame = BeamEnergy::E_beam_corr();	
		beamEnergyError = BeamEnergy::E_beam_err();
		lerBeamEnergy = BeamEnergy::E_LER();
		herBeamEnergy = BeamEnergy::E_HER();
		kekbBeamEnergy = BeamEnergy::E_beam_orig();
		kekbLerBeamEnergy = BeamEnergy::E_LER_orig();
		kekbHerBeamEnergy = BeamEnergy::E_HER_orig();
		beamCrossingAngle = BeamEnergy::Cross_angle();
		cmBoostVector = -BeamEnergy::CMBoost();

		// Initialize PID functions.
		eid::init_data();
		
		// Get interaction point profile data from $BELLE_POSTGRES_SERVER. 
		IpProfile::begin_run();
		
		// Set interaction point and error to run values.
		if ( IpProfile::usable() ) {
			ip = IpProfile::e_position();
			ipErr = IpProfile::e_position_err_b_life_smeared();
			ipUsable = 1;
	    } 
		else {
			ip = HepPoint3D( 0, 0, 0 );
			ipErr = HepSymMatrix( 3, 0 );
	 	}
		
		// Print run information to the log.
		std::cout << "\n\n*********************** Run Information ***********************" << std::endl;
		if ( rhd->ExpMC() == 1 ) {
			flagMC = false;	// Set Data type flag to Real Data.
			std::cout << "Data is Real." << std::endl;
		}
		else {
			flagMC = true;	// Set Data type flag to Monte Carlo.
			std::cout << "Data is Monte Carlo." << std::endl;
		}
		
		std::cout << "Experiment " << experimentNumber << ", Run " << runNumber << "\n" 
			  << "Actual Beam Energy: " << beamEnergyCMFrame << " +/- " << beamEnergyError << " GeV\n"
			  << "Reported Beam Energy: " << kekbBeamEnergy << " GeV\n"
			  << "BE Class cmBoostVector: " << cmBoostVector << "\n"
			  << std::endl;
	}
	
	return;
}


void Bsdilepton::end_run(BelleEvent* evptr, int *status )
{	
	(void)evptr;
	(void)status;
	
	std::cout << "*************WHERE IS THIS FUNCTION end_run() EXECUTED?!?!?!*********" << std::endl;

	return;
}

/* 2010.07.27 - This function is no longer used.
// Determines actual beam energy for run from runhead information.
void Bsdilepton::setBeamInformation()
{	
	// Define the beam crossing angle.
	const double theta( 0.022 );
	
	// WARNING - 2010.06.29 - Energies from the runhead may not be reliable.
	//						  Consider using Benergy() functions which 
	//						  retrieve information about the energies 
	//						  calculated from B meson data.
	//
	// Get LER and HER energies from BELLE runhead.
	Belle_runhead_Manager &HeadMgr = Belle_runhead_Manager::get_manager();
	Belle_runhead &brun = *HeadMgr.begin();
	double herEnergy = brun.EHER();
	double lerEnergy = brun.ELER();
	
	// Set total LER and HER 4 momentum class variable. The LER travels in the lab-frame negative z direction.
	// and the HER lies in the x-z plane, making a the crossing angle of 22 milliradians.
	cmFourMomentum.setE( herEnergy + lerEnergy );
	cmFourMomentum.setPx( herEnergy * sin( theta ) );
	cmFourMomentum.setPy( 0 );
	cmFourMomentum.setPz( ( herEnergy * cos( theta ) ) - lerEnergy );
	
	// Set class variable cmBoostVector as boost vector to CM frame. Since vectors above are with proper
	// signs for the lab frame, we must reverse the sign of the boost vector in order to
	// move to from the lab frame to the CM frame.
	Hep3Vector cmBoostVectorTEMP = -cmFourMomentum.boostVector();
	std::cout << "Runhead cmBoostVector =" << cmBoostVectorTEMP << std::endl;
	
	// Compute (effective) beam energy in the CM frame and store to class variable.
	// Beam energy in CM frame is defined as Total CM energy divided by two.
	// See BELLE note 418 - T. Matsumoto and BELLE Note 567 - T. Browder.
	cmFourMomentum.boost( cmBoostVectorTEMP );
	std::cout << "Runhead beamEnergyCMFrame =" << cmFourMomentum.mag() / 2 << std::endl;
	
}
*/

// This is the primary function of the module.
// Bsdilepton::event() finds analyses each event to find 
// signal dilepton event candidates and records information
// about those events to the n-tuple.
void Bsdilepton::event(BelleEvent* evptr, int* status)
{
	(void)evptr;
	(void)status;
	
	// Instantiate constant data structures storeed in external header files.
   	PDGmasses masses;
	BsdileptonCuts cuts;
	
	// Add event to counter.
	++numberOfEvents;
	
	// Set particle counters.
	int numCndtEPlus = 0;
	int numCndtEMinus = 0;
	int numCndtMuPlus = 0;
	int numCndtMuMinus = 0;

	// Set dilepton(++/--/total) counters.
	int numCndtEPlusPlus = 0;
	int numCndtMuPlusPlus = 0;
	int numCndtEMinusMinus = 0;
	int numCndtMuMinusMinus = 0;
	int numDileptonEvents = 0;

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
	//    relocating list to a new block of memory in case of large
	//    events.
	static std::vector< Particle > uncleanElectronList( 100 );
	static std::vector< Particle > electronList( 100 );
	static std::vector< Particle > muonList( 100 );
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
	kaonList.clear();
	dsPlusList.clear();
	dsMinusList.clear();
	dsStarPlusList.clear();
	dsStarMinusList.clear();
	dileptonEventList.clear();
	
	// Check that beam energy and CM boost vector are accessible.
	//std::cout << "\n Run beam energy (lab frame): " << TotalBeamEnergyLabFrame
	//		  << "\n Run beam energy (CM frame): " << cmFourMomentum.mag() 
	//		  << "\n Boost vector to CM: " << cmBoostVector << ".\n" << std::endl;
	
	// Alias the MDST charged manager, which contains
	// the measured charged tracks for each event.
	Mdst_charged_Manager &chg_mgr = Mdst_charged_Manager::get_manager();
	
	// Create an alias "null" for the null value of generated signal MC.
	// I.E. if something is not signal MC, it will be ascribed this value.
	const Gen_hepevt &null = Gen_hepevt_Manager::get_manager().get_NULL();
	
	// Populate the lepton candidate lists.
	for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin(); i != chg_mgr.end(); ++i ) {
		
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
		
     	// Cut on IP dr and dz - this is to make sure that particles were created near the IP.
		IPdrdz ipDrDzParameters;
		ipDrDzParameters.SetIPdrdz( chg, ip );
		if ( abs( ipDrDzParameters.Getdr() ) > cuts.maxIpDr || abs( ipDrDzParameters.Getdz() ) > cuts.maxIpDz )
			continue;
			
		// WARNING - 2010.28. - Not sure what mhyp( 3 ) indicates.
		// Reject particles with too few hits in the SVD.
		Mdst_trk_fit &tfit = chg.trk().mhyp( 3 );
		if ( tfit.nhits( 3 ) < cuts.minSvdRHits )       // Hits on the r - phi side.
			continue;
		if ( tfit.nhits( 4 ) < cuts.minSvdZHits )       // Hits on the z side.
			continue;

		// Reject particle if liklihoods are the same from eid and muid.
		if ( eidProb == muidProb ) {
			continue;
			// atc_pid selEMu( 3, 1, 5, 0, 2 );
			// double atcEMuProb = selEMu.prob( chg );
			
		}

		// Create an alias "isSigMC" for the MC information about chg to be used in if() statements.
		// Three cases exist for value of alias:
		// Data is MC, chg is signal MC  -> Alias is NOT null.
		// Data is MC, chg is background -> Alias is null.
		//    		    Data is Real -> Alias is null.
		const Gen_hepevt &isSigMC = ( flagMC ) ? get_hepevt( chg ) : null;
		
		// Assume the particle species is that of the highest PID liklihood.
		if ( eidProb > muidProb ) {

			// Assuming particle is an electron.
			Particle eCndt( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );

			// Cut on track momentum.
			if ( eCndt.p().vect().mag() < cuts.minElectronMomentum )
				continue;
			
			// If signal MC, add MC truth information to electron candidate.
			if ( isSigMC ) {
				eCndt.relation().genHepevt( isSigMC );
				chg.charge() > 0 ? ++numCrctEPlus : ++numCrctEMinus;
			}
			
			// Add eCandidate to list of e+/- candidates.
			uncleanElectronList.push_back( eCndt );
		}
		else {
		
			// Assuming particle is a muon.
			Particle muCandidate( chg, chg.charge() > 0 ? ptypeMuPlus : ptypeMuMinus );

			// Cut on track momentum.
			if ( muCandidate.p().vect().mag() < cuts.minMuonMomentum )
				continue;
			
			// If signal MC, add MC truth information to electron candidate.
			if ( isSigMC ) {
				muCandidate.relation().genHepevt( isSigMC );
				chg.charge() > 0 ? ++numCrctMuPlus : ++numCrctMuMinus;
			}
			
			// Add eCandidate to list of e+/- candidates.
			muonList.push_back( muCandidate );
		}
	
	} // End for() loop populating lepton lists.
	
	// Report number of electons in electron list for diagnostic purposes.
	// std::cout << "Number of electrons in event: " << uncleanElectronList.size() << std::endl;

	// Remove possible pair production electrons.
	for ( std::vector< Particle >::iterator j = uncleanElectronList.begin(); j != uncleanElectronList.end(); ++j ) {
		
		// Make a flag for electrons that are not candidates for pair production daughters.
		// By default, assume all electrons are good.
		bool flagNonPairProduction = true;

		const Particle &eCndt = *j;

		// If the invariant mass of an electron candidate and each other opposite charged tracks
		// is smaller than the cut value (nominally 100 MeV), the electron candidate is rejected.
		for ( std::vector< Mdst_charged >::const_iterator i = chg_mgr.begin(); i != chg_mgr.end(); ++i ) {

			const Mdst_charged &chg = *i;
			Particle otherChg( chg, chg.charge() > 0 ? ptypeEPlus : ptypeEMinus );

			// Reject same sign charge pairs.
			if ( eCndt.charge() == otherChg.charge() )
				continue;
			
			// double ePlusEMinusMass = ( eCndt.p() + otherChg.p() ).invariantMass();
			HepLorentzVector eCndtP = eCndt.p();
			HepLorentzVector otherChgP = otherChg.p();
			
			double ePlusEMinusMass = sqrt( 2 * ( ( masses.electron ) * ( masses.electron ) + ( ( eCndtP.e() * otherChg.e() ) - eCndtP.vect().dot( otherChgP.vect() ) ) ) );
			
			// Report ePlusEMinusMass for diagnostic purposes.
			// std::cout << "e p/m mass: " << ePlusEMinusMass << std::endl;

			if ( ePlusEMinusMass < cuts.minEPlusEMinusMass )
				flagNonPairProduction = false;

		}

		// If electron candidate is not a candidate for pair production,
		// add electron candidate to list of good electrons.
		if ( flagNonPairProduction )
			electronList.push_back( eCndt );

	}	
	
	// Report size of cleaned electron list for diagnotic purposes.
	// std::cout << "Number of electrons in event: " << electronList.size() << std::endl;
	
	// Find dielectron event candidates.
	// Loop over the electron list.
	for ( std::vector< Particle >::iterator j = electronList.begin(); j != electronList.end(); ++j ) {
		
		// Loop over remaining electrons in the list to check all electron pairs.
		for ( std::vector< Particle >::iterator i = j; i != electronList.end(); ++i ) {
			
			// Exclude the case where both iterators point to the same particle.
			if ( i == j )
				continue;
			
			// Alias the first electron as "electron0"
			//  and the second electron as "electron1".
			Particle &electron0 = *j;
			Particle &electron1 = *i;
				
			// Require electrons have same sign.
			// if ( electron0.charge() != electron1.charge() )
			//	continue;
			
			// Cut on each lepton's CM momentum.
			// TODO - 2010.08.10 - Perhaps this should be done when populating the lepton lists.
			// Get the lab-frame momentum for each lepton and boost it to CM frame.
			HepLorentzVector lepton0CmMomentum( electron0.p() );
			HepLorentzVector lepton1CmMomentum( electron1.p() );
			lepton0CmMomentum.boost( cmBoostVector );
			lepton1CmMomentum.boost( cmBoostVector );

			// Make the CM momentum cut.
			if ( lepton0CmMomentum.vect().mag() < cuts.minLeptonCmMomentum || lepton0CmMomentum.vect().mag() > cuts.maxLeptonCmMomentum )
				continue;

			if ( lepton1CmMomentum.vect().mag() < cuts.minLeptonCmMomentum || lepton1CmMomentum.vect().mag() > cuts.maxLeptonCmMomentum )
				continue;
			
			// TODO - 2010.07.27 - Cut on jets.
			// 					 - Hastings makes a cut requiring hits in the barrel CsI. Is this needed?
			
			// Make an event candidate.
			DileptonEvent eventCandidate( electron0, electron1 );

			// Add event candidate to the list of dilepton events.
			dileptonEventList.push_back( eventCandidate );
			++numDileptonEvents;

		}
	}

	// Report the size of the dilepton event list for diagnostic purposes.
	// std::cout << "dilepton event list length is " << dileptonEventList.size() << std::endl;

	// TODO - 2010.07.22 - Loop over same sign lepton candidates list and write information to n-tuple.

	return;
}


// Specifies n-tuples to write and names for root histograms.
void Bsdilepton::hist_def()
{
	extern BelleTupleManager *BASF_Histogram;		// Define a BASF Histogram
	BelleTupleManager *tm = BASF_Histogram;		
	const char *ntlist = "mbs";
	m_nt = tm->ntuple( "Bs Dilepton", ntlist, 1 );
	
	return;
}


#if defined(BELLE_NAMESPACE)	// Needed to close compiler namespace
} // namespace Belle			// container for backwards compatibility
#endif							// with older BELLE Libraries.
