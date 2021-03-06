//______________________________________________________________________________
// Filename: TrackParameters.h
// Version: 2011.08.07.A
// Author: M.P. Belhorn
// Original Date: 2010.06.24
// Description: A class for managing track parameterization.
//______________________________________________________________________________

#include "TrackParameters.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default constructor.
TrackParameters::TrackParameters()
    : dr_(new double(-3000)), // Units of cm. Initilized way outside detector.
      dz_(new double(-3000)), // Value of dr & dz can be used track errors.
      mass_hypothesis_(new int(2)),
      helix_(new Helix(HepPoint3D(0, 0, 0), Hep3Vector(0, 0, 1), 1))
{ 
  // Intentionally blank.
}

// Useful constructor.
TrackParameters::TrackParameters(const Particle& particle, HepPoint3D new_pivot)
    : dr_(new double(-3000)), // Units of cm. Initilized way outside detector.
      dz_(new double(-3000)), // Value of dr & dz can be used track errors.
      mass_hypothesis_(new int(2)),
      helix_(new Helix(HepPoint3D(0, 0, 0), Hep3Vector(0, 0, 1), 1))
{
  // Get Mdst_trk from particle object if it exists.
  Mdst_trk &mdst_track = particle.mdstCharged().trk();
  if (!mdst_track) {
    *dr_ = -2000;
    *dz_ = -2000;
  } else {
    // Mdst_trk member function mhyp( int hypID ) returns the fitted
    //   track parameters assuming certain particle mass hypotheses set
    //   by hypID. Possible mass hypotheses are:
    //   hypID = 0:e; 1:mu; 2:pi; 3:K; 4:p. 
    if (abs(particle.lund()) == 11) *mass_hypothesis_ = 0;
    else if (abs(particle.lund()) == 13) *mass_hypothesis_ = 1;
    else if (abs(particle.lund()) == 321) *mass_hypothesis_ = 3;
    else if (abs(particle.lund()) == 2212) *mass_hypothesis_ = 4;

    Mdst_trk_fit &track_fit = mdst_track.mhyp(*mass_hypothesis_);
    if (!track_fit) {
      *dr_ = -1000;
      *dz_ = -1000;
    } else {
      // By default, pivot of fit is assumed to be first hit wire in CDC.
      HepPoint3D pivot(track_fit.pivot_x(),
          track_fit.pivot_y(), track_fit.pivot_z());

      // See Belle Note "Track Parameterization" - Yukiyoshi Ohnishi.
      HepVector helix_parameters(5, 0);
      helix_parameters[0] = track_fit.helix(0); // dr
      helix_parameters[1] = track_fit.helix(1); // phi_0
      helix_parameters[2] = track_fit.helix(2); // kappa
      helix_parameters[3] = track_fit.helix(3); // dz
      helix_parameters[4] = track_fit.helix(4); // tan(lambda)

      HepSymMatrix helix_errors(5, 0);
      helix_errors[0][0] = track_fit.error(0);
      helix_errors[1][0] = track_fit.error(1);
      helix_errors[1][1] = track_fit.error(2);
      helix_errors[2][0] = track_fit.error(3);
      helix_errors[2][1] = track_fit.error(4);
      helix_errors[2][2] = track_fit.error(5);
      helix_errors[3][0] = track_fit.error(6);
      helix_errors[3][1] = track_fit.error(7);
      helix_errors[3][2] = track_fit.error(8);
      helix_errors[3][3] = track_fit.error(9);
      helix_errors[4][0] = track_fit.error(10);
      helix_errors[4][1] = track_fit.error(11);
      helix_errors[4][2] = track_fit.error(12);
      helix_errors[4][3] = track_fit.error(13);
      helix_errors[4][4] = track_fit.error(14);

      Helix helix(pivot, helix_parameters, helix_errors);
      helix.pivot(new_pivot);

      *dr_ = helix.dr();
      *dz_ = helix.dz();
    }
  }
}

// Copy constructor.
TrackParameters::TrackParameters(const TrackParameters &that)
    : dr_(new double(*that.dr_)),
      dz_(new double(*that.dz_)),
      mass_hypothesis_(new int(*that.mass_hypothesis_)),
      helix_(new Helix(*that.helix_))
{
  // Intentionally blank.
}

// Assignment operator.
TrackParameters &TrackParameters::operator= (const TrackParameters &that)
{
  if (this != &that) {
    *dr_ = *(that.dr_);
    *dz_ = *(that.dz_);
    *mass_hypothesis_ = *(that.mass_hypothesis_);
    *helix_ = *(that.helix_);
  }
  return *this;
}

// Destructor.
TrackParameters::~TrackParameters()
{
  delete helix_;
  delete mass_hypothesis_;
  delete dz_;
  delete dr_;
}

// Accessor for dr_ helix parameter.
double &TrackParameters::dr()
{
  return *dr_;
}

// Accessor for dz_ helix parameter.
double &TrackParameters::dz()
{
  return *dz_;
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
