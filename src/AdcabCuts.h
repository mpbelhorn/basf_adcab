//
//******************************************************************************
// Filename: AdcabCuts.h
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010.07.20
// Description: Data structure containing Bs -> dilepton selection cuts.
//******************************************************************************

#ifndef CONSTANTS_H
#define CONSTANTS_H

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif
  
// Particle Selection Cuts.
struct AdcabCuts {
  
  // General charged track cuts.  
  static const int minSvdRHits            = 1;   // Nominally 1.
  static const int minSvdZHits            = 2;   // Nominally 2.   
  static const double maxIpDr             = 0.05; // cm.
  static const double maxIpDz             = 2.0; // cm.
  static const double minLeptonMomentumCm = 0.0; // = 1.1;  // GeV/c.
  static const double maxLeptonMomentumCm = 5.0; // = 2.3;  // GeV/c.
  static const double minLeptonCosTheta   = -1.0;
  static const double maxLeptonCosTheta   = 1.0;
  
  // Electron selection cuts.
  static const double minEidProb          = 0.80; // Low for diagnostics.
  static const double minElectronMomentum = 0.00; // GeV/c
  
  // Muon selection cuts.
  static const double minMuidProb         = 0.90; // Low for diagnostics.
  static const double minMuonMomentum     = 0.00; // GeV/c

  // Pair-production removal cuts.
  static const double minEPlusEMinusMass  = 0.10; // GeV/c

  // Dilepton event selection cuts.
  static const double minCosThetaLLCm     = -1; // = -0.80;
  static const double maxCosThetaLLCm     =  1; // =  0.95;

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
