//
//*******************************************************************************
// Filename: AnalysisTools.h
// Version: 2010.06.24.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Declaration of custom analysis classes and functions. 
//*******************************************************************************

#ifndef ANALYSISTOOLS_H
#define ANALYSISTOOLS_H

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
#include "BsdileptonCuts.h"				// Analysis specfic selection cut constants.

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif


//*******************************************************************************
// IPdrdz class definition and prototypes.
//*******************************************************************************

// Class for Impact Parameters "dr" and "dz"
class IPdrdz {

public:
	IPdrdz();		// Default constructor.
	~IPdrdz() {}	// Destructor.
	
	void SetIPdrdz( const Mdst_charged&, HepPoint3D );	// Sets IP-pivot track helix dr and dz parameters.
	double Getdr();										// Access IP-pivot dr parameter.
	double Getdz();										// Access IP-pivot dz parameter.
	
	
private:
	bool set;
	double dr;
	double dz;
	
};


//*******************************************************************************
// Helicity Cut class definition and prototypes.
//*******************************************************************************

class RhoHelicityCut {

public:
	RhoHelicityCut();																// Default constructor.
	RhoHelicityCut( HepLorentzVector, HepLorentzVector, HepLorentzVector );			// Momentum Constructor.
	~RhoHelicityCut() {}															// Default destructor.
	
	void CalcCosTheta12( HepLorentzVector, HepLorentzVector, HepLorentzVector );
	double GetCosTheta12();
	
private:
	double CosTheta12;
	Hep3Vector RhoCMBoostVector;
	Hep3Vector B0ThreeMomentumInRhoFrame;
	Hep3Vector Pi0ThreeMomentumInRhoFrame;
	
};


//*******************************************************************************
// DileptonEvent class definition and prototypes.
//*******************************************************************************

class DileptonEvent {

public:
	DileptonEvent();						// Default constructor.
	DileptonEvent( Particle, Particle );	// Dilepton Particle constructor.
	~DileptonEvent() {}						// Default destructor.

	void setLepton0( Particle );			// Mutator for lepton0_.
	void setLepton1( Particle );			// Mutator for lepton1_.
	Particle lepton0();						// Accessor for lepton0_.
	Particle lepton1();						// Accessor for lepton1_.

private:
	// Private attributes.
	Particle lepton0_;
	Particle lepton1_;

};


//*******************************************************************************
// General function prototypes.
//*******************************************************************************

// Determines if a photon was measured in CsI barrel or endcap.
int CsIHitLocation( HepLorentzVector GammaMomentum );

// Determines whether to reject a neutral pion candidate due to photon hit energy and location in CsI calorimeter.
bool CsIPhotonCut( HepLorentzVector GammaMomentum, double MinCsIBarrelEnergy, double MinCsIEndcapEnergy );

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
#endif
