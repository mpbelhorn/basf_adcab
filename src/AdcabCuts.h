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
  // Cut values here are low. Refinements are made at the ROOT level.
  // General charged track cuts.
  static const int minSvdRHits            = 1;     // Nominally 1.
  static const int minSvdZHits            = 2;     // Nominally 2.
  static const double maxIpDr             = 0.05;  // cm.
  static const double maxIpDz             = 2.00;  // cm.
  static const double minLeptonMomentumCm = 0.0;   // = 1.1;  // GeV/c.
  static const double maxLeptonMomentumCm = 5.0;   // = 2.3;  // GeV/c.
  static const double minLeptonCosTheta   = -1.0;  // Nominally in barrel PID.
  static const double maxLeptonCosTheta   =  1.0;  // Nominally in barrel PID.
  static const double massJPsi            = 3.096; // Gev

  // Electron selection cuts.
  static const double minEidProb           = 0.80;  // Low for diagnostics.
  static const double minElectronMomentum  = 0.00;  // GeV/c
  static const double min_dielectron_jpsi_mass_difference = -0.15; // Gev
  static const double max_dielectron_jpsi_mass_difference =  0.05; // Gev

  // Muon selection cuts.
  static const double minMuidProb          = 0.90;  // Low for diagnostics.
  static const double minMuonMomentum      = 0.00;  // GeV/c
  static const double min_dimuon_jpsi_mass_difference = -0.05; // Gev
  static const double max_dimuon_jpsi_mass_difference =  0.05; // Gev
  static const double maxKlmChi2PerHits    = 10.00; // Gev

  // Pair-production removal cuts.
  static const double minEPlusEMinusMass  = 0.10; // GeV/c

  // Kaon selection cuts.
  static const double minKaonToPionLikelihood   = 0.6;
  static const double minKaonToProtonLikelihood = 0.6;

  // Dilepton event selection cuts.
  static const double minCosThetaLLCm     = -1.0; // = -0.80;
  static const double maxCosThetaLLCm     =  1.0; // =  0.95;

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
