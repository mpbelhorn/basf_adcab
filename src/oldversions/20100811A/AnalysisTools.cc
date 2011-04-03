//
//*******************************************************************************
// Filename: AnalysisTools.cc
// Version: 2010.06.24.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Definitions of custom analysis classes and functions. 
//*******************************************************************************

#include <cmath>						// Uses cmath functions.

#include "belle.h"						// BELLE Library.
#include "event/BelleEvent.h"			// For managing BELLE events.
#include "tuple/BelleTupleManager.h"	// For managing BELLE nTuples.
#include "benergy/BeamEnergy.h"			// For determining run beam energy.
#include "basf/module.h"				// For interfacing with BASF.
#include "basf/module_descr.h"			// ???
#include "particle/Particle.h"			// The BELLE Particle Class.
#include "kid/atc_pid.h"				// For particle species separation.
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

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//*******************************************************************************
// IPdrdz member function definitions.
//*******************************************************************************

// Constructor.
IPdrdz::IPdrdz()
{	
	set = false;
	dr = -44444;	// Units of cm. Initilized way outside detector.
	dz = -44444;	// Units of cm. Initilized way outside detector.
}


// Interaction Point pivot dr and dz parameters.
	// Obtains a charged particle's 5 fitted helix parameters (See "Track Parameterization" - Yukiyoshi Ohnishi).
	// Helix parameters are:
	//			 a = ( dr, phi_0, kappa, dz, Tan(lambda) ) where
	//          dr = displacement of helix from pivot point,
	//	     phi_0 = azimuthal angle to specify the pivot w.r.t the helix center,
	//       kappa = Inverse of transverse momentum where sign(kappa) indicates assumed charge,
	//          dz = displacement of helix from pivot point in z-direction.
	// tan(lambda) = slope of the track (tangent of dip angle).
	// By default, pivot of fit is assumed to be first hit wire in CDC. Here, we reparametrize the
	//		trajectory using the IP as the pivot in order to determine (using dr and dz) how close the
	//		charged tracks originate to the decay of the B meson, which should be close to the IP.
void IPdrdz::SetIPdrdz( const Mdst_charged& chg, HepPoint3D IP )
{	
	// Get particle's track assuming mass hypothesis 3 (???).
	Mdst_trk &trk = chg.trk();
	Mdst_trk_fit &trkfit = trk.mhyp( 3 );
	
	// Obtain the fitted CDC helix parameters and CDC pivot point.
	HepVector CDCHelixParameters( 5, 0 );
	for ( int i = 0; i < 5; i++ ) {
		CDCHelixParameters[ i ] = trkfit.helix( i );
	}
	HepPoint3D CDCPivot( trkfit.pivot( 0 ), trkfit.pivot( 1 ), trkfit.pivot( 2 ) );
	
	// Create a new set of helix parameters from the old ones.
	Helix IPHelixParameters( CDCPivot, CDCHelixParameters );
	
	// Transform the new parameters into a set using the IP as the pivot point.
	IPHelixParameters.pivot( IP );
	
	// Set the values of IPdrdz class dr and dz from the IP-pivot parameterization,
	// and indicate that they have been set correctly.
	dr = IPHelixParameters.dr();
	dz = IPHelixParameters.dz();
	set = true;

}

// Access IP-pivot dr helix parameter.
double IPdrdz::Getdr()
{
	return dr;
}

// Access IP-pivot dz helix parameter.
double IPdrdz::Getdz()
{
	return dz;
}



//*******************************************************************************
// HelicityCut member function definitions.
//*******************************************************************************

// Default Constructor.
RhoHelicityCut::RhoHelicityCut()
{	
	CosTheta12 = 0;
}

// Momentum Constructor.
RhoHelicityCut::RhoHelicityCut( HepLorentzVector pLabRho, HepLorentzVector pLabPi0, HepLorentzVector pLabB0)
{	
	CalcCosTheta12( pLabRho, pLabPi0, pLabB0 );
}

// Calculates and sets CosTheta12.
void RhoHelicityCut::CalcCosTheta12( HepLorentzVector pLabRho, HepLorentzVector pLabPi0, HepLorentzVector pLabB0)
{	
		RhoCMBoostVector = -( pLabRho ).boostVector();
		
		HepLorentzVector B0FourMomentumInRhoFrame = pLabB0;
		B0FourMomentumInRhoFrame.boost( RhoCMBoostVector );
		B0ThreeMomentumInRhoFrame = B0FourMomentumInRhoFrame.vect();
		
		HepLorentzVector Pi0FourMomentumInRhoFrame = pLabPi0;
		Pi0FourMomentumInRhoFrame.boost( RhoCMBoostVector );
		Pi0ThreeMomentumInRhoFrame = Pi0FourMomentumInRhoFrame.vect();
		CosTheta12 = ( ( -B0ThreeMomentumInRhoFrame ).dot( Pi0ThreeMomentumInRhoFrame ) ) / ( ( B0ThreeMomentumInRhoFrame.mag() ) * ( Pi0ThreeMomentumInRhoFrame.mag() ) );
}

double RhoHelicityCut::GetCosTheta12()
{	
	return CosTheta12;
}

//*******************************************************************************
// Definitions for use with Bs->Dilepton analysis.
//*******************************************************************************

// Default constructor.
DileptonEvent::DileptonEvent()
{
	// Function is blank.
}

// Dilepton constructor.
DileptonEvent::DileptonEvent( Particle lepton0, Particle lepton1 )
{
	lepton0_ = lepton0;
	lepton1_ = lepton1;
}

// Mutator for lepton0_.
void DileptonEvent::setLepton0( Particle lepton0 )
{
	lepton0_ = lepton0;
}

// Mutator for lepton1_.
void DileptonEvent::setLepton1( Particle lepton1 )
{
	lepton1_ = lepton1;
}

// Accessor for lepton0_.
Particle DileptonEvent::lepton0()
{
	return lepton0_;
}

// Accessor for lepton1_.
Particle DileptonEvent::lepton1()
{
	return lepton1_;
}

//*******************************************************************************
// General analysis function definitions.
//*******************************************************************************

// Determines if a photon was measured in CsI barrel or endcap.
int CsIHitLocation( HepLorentzVector GammaMomentum )
{
	// Start by gettting the photon's lab-frame polar (Theta) cosines.
	double GammacosTheta = GammaMomentum.cosTheta();
	
	// From polar cosine, determine which part of the calorimeter measured the photon.
	// For backward endcap face, Theta > 130.7 Degrees		 ->			   Cos(Theta) < -0.6521
	// For barrel face, 128.7 Degrees > Theta > 32.2 Degrees ->	 -0.6252 < Cos(Theta) <  0.8462
	// For forward endcap face, 31.4 Degrees > Theta		 ->	  0.8535 < Cos(Theta)
	// Set flag GammaCsIHitLocation = 3 (Not determined), = 2 (Backward Endcap), = 1 (Barrel), = 0 (Forward Endcap).
	int GammaCsIHitLocation = 3;
	if (GammacosTheta < -0.6521)
		GammaCsIHitLocation = 2;
	if ( ( -0.6252 < GammacosTheta ) && ( GammacosTheta < 0.8462 ) )
		GammaCsIHitLocation = 1;
	if ( 0.8535 < GammacosTheta )
		GammaCsIHitLocation = 0;
		
	return GammaCsIHitLocation;
}

// Determines whether to reject a neutral pion candidate due to photon hit energy and location in CsI calorimeter.
bool CsIPhotonCut( HepLorentzVector GammaMomentum, double MinCsIBarrelEnergy, double MinCsIEndcapEnergy )
{	
	// Find what part of the CsI the photon was measured in.
	int HitLocation = CsIHitLocation( GammaMomentum );
	
	// Cut all photons by default.
	bool CutPhoton = true;
	
	// Include photons based on thier energy and CsI hit location.
	switch ( HitLocation ) {
	
		// Case if measured in CsI Forward Endcap.
		case 0:
			// Reject iff E(gamma)_endcap < MinCsIEndcapEnergy.
			if ( GammaMomentum.e() > MinCsIEndcapEnergy )
				CutPhoton = false;
			break;
			
		// Case if measured in CsI Barrel.
		case 1:
			// Reject iff E(gamma)_barrel < MinCsIBarrelEnergy.
			if ( GammaMomentum.e() > MinCsIBarrelEnergy )
				CutPhoton = false;
			break;
			
		// Reject all other hit cases.
		case 2:	// Case if measured in CsI Backward Endcap.
		case 3:	// Case if CsI hit location is indeterminate.
		default:
			break;
	}
	
	return CutPhoton;
	
}



#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
