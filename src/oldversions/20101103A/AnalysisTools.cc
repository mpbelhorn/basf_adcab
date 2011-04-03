//
//******************************************************************************
// Filename: AnalysisTools.cc
// Version: 2010.09.20.A
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
#include "BsdileptonCuts.h"       // Analysis specfic selection cut constants.

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
// Definitions for use with Bs->Dilepton analysis.
//******************************************************************************

// Default constructor.
DileptonEvent::DileptonEvent()
{
  // Function is blank.
}

// Dilepton constructor. Puts the greater charged lepton in position 0.
DileptonEvent::DileptonEvent( Particle lepton0, Particle lepton1 )
{
  if ( lepton0.charge() > lepton1.charge() ) {
    lepton0_ = lepton0;
    lepton1_ = lepton1;
  } else {
    lepton0_ = lepton1;
    lepton1_ = lepton0;
  }
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

// Calculates the cosine of the opening angle theta_ll between the two leptons
//   in the CM frame.
double DileptonEvent::cosThetaLLCm( Hep3Vector cmBoostVector )
{
  HepLorentzVector lepton0PCm = lepton0_.p();
  HepLorentzVector lepton1PCm = lepton1_.p();
  lepton0PCm.boost( cmBoostVector );
  lepton1PCm.boost( cmBoostVector );
  Hep3Vector lepton0P3Cm = lepton0PCm.vect();
  Hep3Vector lepton1P3Cm = lepton1PCm.vect();
  
  double cosTheta;
  double pTotal2 = lepton0P3Cm.mag2() * lepton1P3Cm.mag2();
  if ( pTotal2 <= 0 ) {
    // This case can only happen if a lepton 3-momenta is not real or the 0
    // vector. Thus, I set the cosine to a value that will be cut.
    cosTheta = 1.0;
  } else {
    cosTheta = lepton0P3Cm.dot( lepton1P3Cm ) / sqrt( pTotal2 );
    if ( cosTheta >  1.0 ) cosTheta =  1.0;
    if ( cosTheta < -1.0 ) cosTheta = -1.0;
  }
  return cosTheta;
}

// Determines the type of dilepton event and returns an integer for the
//   cases of 22:(e+e+)/(e-e-);     -22:e+e-;
//            24:(mu+e+)/(mu-e-);   -24:mu+mu-
//            26:(mu+mu+)/(mu-mu-); -26:(mu+e-)/(mu-e+).
// Negative codes imply opposite-sign whereas positive codes imply same-sign.
// The magnitude of the code is the sum of the magnitudes of the Lund ID values.
// The Lund values as of 2010/08/12 are e(+/-) = (+/-)11, mu(+/-) = (+/-)13.
double DileptonEvent::eventType()
{
  double sign = lepton0_.charge() * lepton1_.charge();
  double lundSum = abs( lepton0_.lund() ) + abs( lepton1_.lund() );
  return sign * lundSum;
}

// Determines if the mothers of both leptons are MC signal Bs's (Pythia ID 531).
// Returns 1 if true, 0 if not true, and -1 if not both MC info available. 
int DileptonEvent::mcBsParents()
{
  if ( lepton0_.relation().genHepevt() && lepton1_.relation().genHepevt() ) {
    if ( lepton0_.relation().genHepevt().mother() &&
         lepton1_.relation().genHepevt().mother() ) {
      double l0MotherId = lepton0_.relation().genHepevt().mother().idhep();
      double l1MotherId = lepton1_.relation().genHepevt().mother().idhep();
      if ( abs( l0MotherId ) == abs( l1MotherId ) &&
           abs( l0MotherId ) == 531 ) {
        return 1;
      }
    } else {
      return 0;
    }
  } else {
    return -1;
  }
}

// Returns the LUND ID for lepton0_. If no MC truthtable info
//  exists, then 0 is returned.
double DileptonEvent::l0Id()
{
  double mcId = 0;
  if ( lepton0_.relation().genHepevt() ) {
    if ( lepton0_.relation().genHepevt().idhep() ) {
      mcId = lepton0_.relation().genHepevt().idhep();
    }
  }
  return mcId;
}

// Returns the LUND ID for lepton1_. If no MC truthtable info
//  exists, then 0 is returned.
double DileptonEvent::l1Id()
{
  double mcId = 0;
  if ( lepton1_.relation().genHepevt() ) {
    if ( lepton1_.relation().genHepevt().idhep() ) {
      mcId = lepton1_.relation().genHepevt().idhep();
    }
  }
  return mcId;
}

// Returns the LUND ID for lepton0_'s MC parent. If no MC truthtable info
//  exists, then 0 is returned.
double DileptonEvent::l0MotherId()
{
  double motherId = 0;
  if ( lepton0_.relation().genHepevt() ) {
    if ( lepton0_.relation().genHepevt().mother() ) {
      motherId = lepton0_.relation().genHepevt().mother().idhep();
    }
  }
  return motherId;
}

// Returns the LUND ID for lepton1_'s MC parent. If no MC truthtable info
//  exists, then 0 is returned.
double DileptonEvent::l1MotherId()
{
  double motherId = 0;
  if ( lepton1_.relation().genHepevt() ) {
    if ( lepton1_.relation().genHepevt().mother() ) {
      motherId = lepton1_.relation().genHepevt().mother().idhep();
    }
  }
  return motherId;
}

//******************************************************************************
// General analysis function definitions.
//******************************************************************************

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
