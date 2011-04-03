//
//*******************************************************************************
// Filename: BsdileptonCuts.h
// Version: 2010.07.20.A
// Author: M.P. Belhorn
// Original Date: 2010.07.20
// Description: Data structure containing Bs -> same sign dilepton selection cuts.
//*******************************************************************************

#ifndef CONSTANTS_H
#define CONSTANTS_H

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif
  
// Particle Selection Cuts.
struct BsdileptonCuts {
	
	// General charged track cuts.	
	static const double maxIpDr		= 0.2;	// Maximum helix displacement from IP in radial direction in cm.
	static const double maxIpDz		= 4.0;	// Maximum trajectory helix displacement from IP in z direction in cm.
	static const double minSvdRHits		= 1;	// Requires at least one hit in the SVD on the r-phi side.
	static const double minSvdZHits		= 2;	// Requires at least two hits in the SVD on the z side.
	static const double minLeptonCmMomentum = 1.1;	// GeV/c.
	static const double maxLeptonCmMomentum = 2.3;	// GeV/c.
	
	// e+/e- selection cuts.
	static const double minEidProb		= 0.2;	// Minimum acceptable eid.pid probability for electron selection.
	static const double minElectronMomentum	= 0.00;	// GeV/c
	static const double minEPlusEMinusMass	= 0.10;	// GeV/c
	
	// Muon selection cuts.
	static const double minMuidProb		= 0.8;	// Minimum acceptable Muid liklihood for muon selection.
	static const double minMuonMomentum	= 0.00; // GeV/c
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
