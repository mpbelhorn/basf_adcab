//
//******************************************************************************
// Filename: AnalysisTools.cc
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Definitions of custom analysis classes and functions. 
//******************************************************************************

#include <cmath>                     // Uses cmath functions.

#include "belle.h"                   // BELLE Library.
#include "event/BelleEvent.h"        // For managing BELLE events.
#include "tuple/BelleTupleManager.h" // For managing BELLE nTuples.
#include "benergy/BeamEnergy.h"      // For determining run beam energy.
#include "basf/module.h"             // For interfacing with BASF.
#include "basf/module_descr.h"       // ???
#include "particle/Particle.h"       // The BELLE Particle Class.
#include "kid/atc_pid.h"             // For particle species separation.
#include "eid/eid.h"                 // For electron identification.
#include "mdst/mdst.h"               // For MDST files.
#include "ip/IpProfile.h"            // Beam Interaction Point (IP) analysis
                                     //   tools. Position unit = cm.

#include <panther/panther.h>      // Panther.
#include BELLETDF_H               // Panther.
#include HEPEVT_H                 // Panther.
#include MDST_H                   // Panther.

#include "HEPconstants.h"         // PDG masses and constants.
#include "AnalysisTools.h"        // General analysis functions and utilities.
#include "AdcabCuts.h"       // Analysis specfic selection cut constants.

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//******************************************************************************
// IpDrDz member function definitions.
//******************************************************************************

// Default constructor.
IpDrDz::IpDrDz()
{ 
  dr_ = -44444;  // Units of cm. Initilized way outside detector.
  dz_ = -44444;  // Units of cm. Initilized way outside detector.
}

// Useful constructor.
IpDrDz::IpDrDz( const Mdst_charged& chg, HepPoint3D ip, int massHyp )
{
  // This constructor follows the same series of calculations as
  //   IpDrDz::setDrDz() but forgoes the use of keyword "this"
  //   and private mutators.
  Mdst_trk_fit &trkfit =chg.trk().mhyp( massHyp );
  HepVector cdcHelixParameters( 5, 0 );
  for ( int i = 0; i < 5; i++ ) {
    cdcHelixParameters[ i ] = trkfit.helix( i );
  }
  HepPoint3D cdcPivot( trkfit.pivot( 0 ),
                       trkfit.pivot( 1 ),
                       trkfit.pivot( 2 ) );
  Helix ipHelixParameters( cdcPivot, cdcHelixParameters );
  ipHelixParameters.pivot( ip );
  dr_  = ipHelixParameters.dr();
  dz_  = ipHelixParameters.dz();
}


// Interaction Point pivot dr and dz parameters.
// Obtains a charged particle's 5 fitted helix parameters
// See "Track Parameterization" - Yukiyoshi Ohnishi.
// Helix parameters are:
//   a           = ( dr, phi_0, kappa, dz, Tan(lambda) ) where
//   dr          = radial displacement of helix from pivot point,
//   phi_0       = azimuthal angle to specify the pivot w.r.t helix center,
//   kappa       = inverse of transverse momentum (sign gives assumed charge),
//   dz          = z-displacement of helix from pivot point.
//   tan(lambda) = slope of the track (tangent of dip angle).
// By default, pivot of fit is assumed to be first hit wire in CDC. Here, we
//   reparametrize the trajectory using the IP as the pivot in order to
//   determine (using dr and dz) how close the charged tracks originate to the
//   decay of the B meson, which should be close to the IP.
void IpDrDz::setDrDz( const Mdst_charged& chg, HepPoint3D ip, int massHyp )
{ 
  // Get the MSDT track information for chg.
  Mdst_trk &trk = chg.trk();

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
  //   only information from the CDC.
  // See mdst.tdf for information about what information is contained in the
  //   mdst files.
  Mdst_trk_fit &trkfit = trk.mhyp( massHyp );
  
  // Obtain the fitted CDC helix parameters and CDC pivot point.
  HepVector cdcHelixParameters( 5, 0 );
  for ( int i = 0; i < 5; i++ ) {
    cdcHelixParameters[ i ] = trkfit.helix( i );
  }
  HepPoint3D cdcPivot( trkfit.pivot( 0 ),
                       trkfit.pivot( 1 ),
                       trkfit.pivot( 2 ) );
  
  // Create a new set of helix parameters from the old ones.
  Helix ipHelixParameters( cdcPivot, cdcHelixParameters );
  
  // Transform the new parameters into a set using the IP as the pivot point.
  ipHelixParameters.pivot( ip );
  
  // Set the values of IPdrdz class dr and dz from the IP-pivot
  //   parameterization, and indicate that they have been set correctly.
  dr_  = ipHelixParameters.dr();
  dz_  = ipHelixParameters.dz();

}

// Accessor for dr_ helix parameter.
double IpDrDz::dr()
{
  return dr_;
}

// Accessor for dz_ helix parameter.
double IpDrDz::dz()
{
  return dz_;
}

//******************************************************************************
// LeptonCandidate Member functions.
//******************************************************************************

// Default Contructor.
LeptonCandidate::LeptonCandidate()
{
  // Function is blank.
}

// Lepton Candidate Constructor. 
LeptonCandidate::LeptonCandidate( Particle lepton, Hep3Vector cmBoostVector )
{
  lepton_ = lepton;
  cmBoostVector_ = cmBoostVector;
}

// Mutator for lepton_.
void LeptonCandidate::setLepton( Particle lepton )
{
  lepton_ = lepton;
}

// Mutator for cmBoostVector_.
void LeptonCandidate::setCmBoostVector( Hep3Vector cmBoostVector )
{
  cmBoostVector_ = cmBoostVector;
}

// Accessor for lepton_.
Particle LeptonCandidate::lepton()
{
  return lepton_;
}

// Accessor for cmBoostVector_.
Hep3Vector LeptonCandidate::cmBoostVector()
{
  return cmBoostVector_;
}

// Alias for for lepton().
Particle LeptonCandidate::particle()
{
  return lepton_;
}

// Returns the pythia particle ID code of the assignment given to particle
//   lepton_ at its creation.
double LeptonCandidate::idAssigned()
{
  return lepton_.pType().lund();
}

// Returns the pythia particle ID code of lepton_ as determined by the MC truth
//   table. Returns 0 if truth table is unavailable.
double LeptonCandidate::idTruth()
{
  double id = 0;
  if ( lepton_.relation().genHepevt() ) {
    if ( lepton_.relation().genHepevt().idhep() ) {
      id = lepton_.relation().genHepevt().idhep();
    }
  }
  return id;
}

// Returns the pythia code of the lepton mother as determined from the truth
//   table. Returns 0 if truth table is unavailable.
double LeptonCandidate::idMother()
{
  double id = 0;
  if ( lepton_.relation().genHepevt() ) {
    if ( lepton_.relation().genHepevt().mother() ) {
      id = lepton_.relation().genHepevt().mother().idhep();
    }
  }
  return id;
}

// Returns the CM frame 4 momentum of lepton_.
HepLorentzVector LeptonCandidate::pCm()
{
  HepLorentzVector leptonPCm( lepton_.p() );
  leptonPCm.boost( cmBoostVector_ );
  return leptonPCm;
}

// Returns the lab frame 4 momentum of lepton_.
HepLorentzVector LeptonCandidate::p()
{
  return lepton_.p();
}

// Returns the eid likelihood of lepton_.
double LeptonCandidate::likelihoodE()
{
  const Mdst_charged &leptonMdstCharged = lepton_.relation().mdstCharged();
  eid leptonEid( leptonMdstCharged );
  return leptonEid.prob( 3, -1, 5 );
}

// Returns the muid likelihood of lepton_.
double LeptonCandidate::likelihoodMu()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.Muon_likelihood();
}

// Returns the Chi^2 of the associated hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
double LeptonCandidate::klmHitsChi2()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.Chi_2();
}

// Returns the number of the hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
int LeptonCandidate::klmHitsN()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.N_layer_hit_brl() + leptonMuid.N_layer_hit_end();
}

// Returns the Chi^2 of the associated hits in the KLM divided by the number of
//   the hits in the KLM assuming the track is a muon. If the track has no
//   associated hits in the KLM layers, the function returns 0 (consistant with 
//   an electron).
double LeptonCandidate::klmHitsChi2PerN()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  int klmHitsN = leptonMuid.N_layer_hit_brl() + leptonMuid.N_layer_hit_end();
  if ( klmHitsN <= 0 ) {
    return 0;
  } else {
    return leptonMuid.Chi_2() / klmHitsN;
  }
}

// Returns the number of hits in the r-side of the SVD for lepton_.
double LeptonCandidate::svdHitsR()
{
  return lepton_.relation().mdstCharged().trk().mhyp( 1 ).nhits( 3 );
}

// Returns the number of hits in the z-side of the SVD for lepton_.
double LeptonCandidate::svdHitsZ()
{
  return lepton_.relation().mdstCharged().trk().mhyp( 1 ).nhits( 4 );
}

//******************************************************************************
// Definitions for use with Bs->Dilepton analysis.
//******************************************************************************

// Default constructor.
DileptonEvent::DileptonEvent()
{
  // Function is blank.
}

// Dilepton constructor.
DileptonEvent::DileptonEvent( Particle lepton0, Particle lepton1, Hep3Vector cmBoostVector )
{
  l0_ = LeptonCandidate( lepton0, cmBoostVector );
  l1_ = LeptonCandidate( lepton1, cmBoostVector );
  cmBoostVector_ = cmBoostVector;
}

// Mutator for l0_.
void DileptonEvent::setL0( Particle lepton0 )
{
  l0_ = LeptonCandidate( lepton0, cmBoostVector_ );
}

// Mutator for l1_.
void DileptonEvent::setL1( Particle lepton1 )
{
  l1_ = LeptonCandidate( lepton1, cmBoostVector_ );
}

// Mutator for cmBoostVector_.
void DileptonEvent::setCmBoostVector( Hep3Vector cmBoostVector )
{
  cmBoostVector_ = cmBoostVector;
}

// Accessor for l0_.
LeptonCandidate DileptonEvent::l0()
{
  return l0_;
}

// Accessor for l1_.
LeptonCandidate DileptonEvent::l1()
{
  return l1_;
}

// Determines the type of dilepton event and returns an integer for the
//   cases of 22:(e+e+)/(e-e-);     -22:e+e-;
//            26:(mu+e+)/(mu-e-);   -26:mu+mu-
//            24:(mu+mu+)/(mu-mu-); -24:(mu+e-)/(mu-e+).
// Negative codes imply opposite-sign whereas positive codes imply same-sign.
// The magnitude of the code is the sum of the magnitudes of the Lund ID values.
// The Lund values as of 2010/08/12 are e(+/-) = (+/-)11, mu(+/-) = (+/-)13.
double DileptonEvent::eventType()
{
  double sign = l0_.particle().charge() * l1_.particle().charge();
  double lundSum = abs( l0_.idAssigned() ) + abs( l1_.idAssigned() );
  return sign * lundSum;
}

// Calculates the cosine of the opening angle theta_ll between the two leptons
//   in the CM frame.
double DileptonEvent::cosThetaLL()
{
  Hep3Vector lepton0P3Cm = l0_.pCm().vect();
  Hep3Vector lepton1P3Cm = l1_.pCm().vect();
  
  double cosThetaLL;
  double pTotal2 = lepton0P3Cm.mag2() * lepton1P3Cm.mag2();
  if ( pTotal2 <= 0 ) {
    // This case can only happen if a lepton 3-momenta is not real or the 0
    // vector. Thus, I set the cosine to a value that will be cut.
    cosThetaLL = 1.0;
  } else {
    cosThetaLL = lepton0P3Cm.dot( lepton1P3Cm ) / sqrt( pTotal2 );
    if ( cosThetaLL >  1.0 ) cosThetaLL =  1.0;
    if ( cosThetaLL < -1.0 ) cosThetaLL = -1.0;
  }
  return cosThetaLL;
}

// Returns the scalar sum of the lepton CM 3-momenta.
double DileptonEvent::pSum()
{
  return l0_.pCm().vect().mag() + l1_.pCm().vect().mag();
}

// Returns the positive scalar difference of the lepton CM 3-momenta.
double DileptonEvent::pDifference()
{
  return abs( l0_.pCm().vect().mag() - l1_.pCm().vect().mag() );
}

//******************************************************************************
// General analysis function definitions.
//******************************************************************************

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
