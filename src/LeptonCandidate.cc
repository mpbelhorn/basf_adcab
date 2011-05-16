//______________________________________________________________________________
// Filename: LeptonCandidate.cc
// Version: 2011.05.15.A
// Author: M.P. Belhorn
// Original Date: 2011.05.15
// Description: Definitions of LeptonCandidate class members.
//______________________________________________________________________________

#include "LeptonCandidate.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Default Contructor.
LeptonCandidate::LeptonCandidate()
{
  // Function is blank.
}

// Lepton Candidate Constructor. 
LeptonCandidate::LeptonCandidate( Particle lepton, Hep3Vector cm_boost )
{
  lepton_ = lepton;
  cm_boost_ = cm_boost;
}

// Mutator for lepton_.
void
LeptonCandidate::set_lepton( Particle lepton )
{
  lepton_ = lepton;
}

// Mutator for cm_boost_.
void
LeptonCandidate::set_cm_boost_vector( Hep3Vector cm_boost )
{
  cm_boost_ = cm_boost;
}

// Accessor for lepton_.
Particle LeptonCandidate::lepton()
{
  return lepton_;
}

// Accessor for cm_boost_.
Hep3Vector
LeptonCandidate::cm_boost()
{
  return cm_boost_;
}

// Returns the pythia particle ID code of the assignment given to particle
//   lepton_ at its creation.
double
LeptonCandidate::idAssigned()
{
  return lepton_.pType().lund();
}

// Returns the mass hypothesis needed for the interaction point parameters.
// mass_hyp = 1 for muons, mass_hyp = 0 for electrons.
int
LeptonCandidate::massHypothesis()
{
  int mass_hyp = 0;
  if (abs(idAssigned()) == 13) {
     mass_hyp = 1;
  }
  return mass_hyp;
}

// Returns the pythia particle ID code of lepton_ as determined by the MC truth
//   table. Returns 0 if truth table is unavailable.
double
LeptonCandidate::idTrue()
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
double
LeptonCandidate::idMom()
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
HepLorentzVector
LeptonCandidate::pCm()
{
  HepLorentzVector leptonPCm( lepton_.p() );
  leptonPCm.boost( cm_boost_ );
  return leptonPCm;
}

// Returns the lab frame 4 momentum of lepton_.
HepLorentzVector
LeptonCandidate::p()
{
  return lepton_.p();
}

// Returns the eid likelihood of lepton_.
double
LeptonCandidate::electronProbability()
{
  const Mdst_charged &leptonMdstCharged = lepton_.relation().mdstCharged();
  eid leptonEid( leptonMdstCharged );
  return leptonEid.prob( 3, -1, 5 );
}

// Returns the muid likelihood of lepton_.
double
LeptonCandidate::muonProbability()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.Muon_likelihood();
}

// Returns the Chi^2 of the associated hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
double
LeptonCandidate::klmHitsChi2()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.Chi_2();
}

// Returns the number of the hits in the KLM assuming the track is a
//   muon. If the track is not in Muid_mdst, the function returns 0.
int
LeptonCandidate::klmHits()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  return leptonMuid.N_layer_hit_brl() + leptonMuid.N_layer_hit_end();
}

// Returns the Chi^2 of the associated hits in the KLM divided by the number of
//   the hits in the KLM assuming the track is a muon. If the track has no
//   associated hits in the KLM layers, the function returns 0 (consistant with 
//   an electron).
double
LeptonCandidate::klmChi2PerHits()
{
  Muid_mdst leptonMuid( lepton_.relation().mdstCharged() );
  if ( klmHits() <= 0 ) {
    return 0;
  } else {
    return klmHitsChi2() / klmHits();
  }
}

// Returns the number of hits in the r-side of the SVD for lepton_.
double
LeptonCandidate::svdRadialHits()
{
  return lepton_.relation().mdstCharged().trk().mhyp( 1 ).nhits( 3 );
}

// Returns the number of hits in the z-side of the SVD for lepton_.
double
LeptonCandidate::svdAxialHits()
{
  return lepton_.relation().mdstCharged().trk().mhyp( 1 ).nhits( 4 );
}

#if defined(BELLE_NAMESPACE)
} //namespace Belle
#endif
