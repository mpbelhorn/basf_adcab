//
//*******************************************************************************
// Filename: HEPconstants.h
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Data structures containing often used HEP constants
//*******************************************************************************

#ifndef HEPCONSTANTS_H
#define HEPCONSTANTS_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// TODO - 2010.06.24 - Check Current Values.
// PDG Particle masses (all GeV/c^2). Not to significant figures - must check by hand.
struct PDGmasses {

  static const double piCharged	   = 0.139570;
  static const double pi0	         = 0.134976;
  static const double kaon         = 0.493677;
  static const double electron     = 0.000511;
  static const double muon         = 0.105658;
  static const double dStarCharged = 2.010000;
  static const double dStar0		   = 2.006700;
  static const double rho				   = 0.775800;
  static const double phi				   = 1.019460;
  static const double d0				   = 1.864500;
  static const double b0			     = 5.279000;
  static const double upsilon4S	   = 10.58000;
  static const double proton		   = 0.938272;
	
};

// Other constants.
struct HEPconstants {

  static const double rad2deg = 57.29577952;
	
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
