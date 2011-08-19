//______________________________________________________________________________
// Filename: AnalysisTools.cc
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010-06-24
// Description: Definitions of custom analysis classes and functions. 
//______________________________________________________________________________

#include "IpParameters.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default constructor.
IpParameters::IpParameters()
{ 
  dr_ = -44444;  // Units of cm. Initilized way outside detector.
  dz_ = -44444;  // Units of cm. Initilized way outside detector.
}

// Useful constructor.
IpParameters::IpParameters(const Mdst_charged& chg, HepPoint3D ip, int massHyp)
{
  init(chg, ip, massHyp);
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
void
IpParameters::init( const Mdst_charged& chg, HepPoint3D ip, int massHyp )
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
  //   only information from the CDC. The mass hypothesis does not affect the
  //   number of svd hits.
  // See mdst.tdf for information about what information is contained in the
  //   mdst files.
  Mdst_trk_fit &trkfit = trk.mhyp(massHyp);
  
  // Obtain the fitted CDC helix parameters and CDC pivot point.
  HepVector cdcHelixParameters(5, 0);
  for (int i = 0; i < 5; i++) {
    cdcHelixParameters[i] = trkfit.helix(i);
  }
  HepPoint3D cdcPivot(trkfit.pivot(0), trkfit.pivot(1), trkfit.pivot(2));
  
  // Create a new set of helix parameters from the old ones.
  Helix ipHelixParameters(cdcPivot, cdcHelixParameters);
  
  // Transform the new parameters into a set using the IP as the pivot point.
  ipHelixParameters.pivot(ip);
  
  // Set the values of IPdrdz class dr and dz from the IP-pivot
  //   parameterization, and indicate that they have been set correctly.
  dr_ = ipHelixParameters.dr();
  dz_ = ipHelixParameters.dz();
  svd_hits_r_ = trkfit.nhits(3);
  svd_hits_z_ = trkfit.nhits(4);
}

// Accessor for dr_ helix parameter.
double
IpParameters::dr()
{
  return dr_;
}

// Accessor for dz_ helix parameter.
double
IpParameters::dz()
{
  return dz_;
}

// Returns the number of hits in the SVD on the r-phi side.
int IpParameters::svdHitsR()
{
  return svd_hits_r_;
}

// Returns the number of hits in the SVD on the r-phi side.
int IpParameters::svdHitsZ()
{
  return svd_hits_z_;
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
