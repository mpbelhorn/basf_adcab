#ifndef USERINFO_H
#define USERINFO_H

//______________________________________________________________________________
// Filename: Userinfo.h
// Author: Unknown.
// Modified: M.P.Belhorn (matt.belhorn@gmail.com).
// Description: Class for custom particle information. Used with A. Zupanc's
//     MC truthtable functions.
// _____________________________________________________________________________

#include "belle.h"
#include "particle/Particle.h"
#include "particle/Momentum.h"
#include "particle/ParticleUserInfo.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "CLHEP/Matrix/Matrix.h"
#include <panther/panther.h>
#include BELLETDF_H
#include HEPEVT_H
#include MDST_H
#include EVTCLS_H

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// For Interface to Set UserInfo Class
void setUserInfo(Particle &p, unsigned ch=0);
void setUserInfo(std::vector<Particle>& p, unsigned ch=0);

// UserInfo Class
class
UserInfo : public ParticleUserInfo
{
 public:
  UserInfo();                 // Default Constructor.
  UserInfo(unsigned);         // Useful Constructor.
  UserInfo(const UserInfo &); // Copy Constructor.
  virtual ~UserInfo();

  // Constructs self object.
  UserInfo * clone(void) const;
  // Copy operator.
  UserInfo & operator = (const UserInfo &);

 public:
  const int & decayMode(void) const { return m_decayMode; }
  void decayMode(const int &v) { m_decayMode = v; }

  const int & genHepevtLink(void) const { return m_genHepevtLink; }
  void genHepevtLink(const int &v) { m_genHepevtLink = v; }

  const int & genHepevtChecked(void) const { return m_genHepevtChecked; }
  void genHepevtChecked(const int &v) { m_genHepevtChecked = v; }

  // Vertex fit Chi^2.
  const double & chisq(void) const { return m_chisq; }
  void chisq(const double &v) { m_chisq = v; }

  /// Vertex fit confidence level.
  const double & cl(void) const { return m_cl; }
  void cl(const double &v) { m_cl = v; }

  /// Vertex fit degrees of freedom.
  const unsigned & ndf(void) const { return m_ndf; }
  void ndf(const unsigned &v) { m_ndf = v; }

  // Muon ID
  void muonId(const Muid_mdst &muid) {
    muon_likelihood_ = muid.Muon_likelihood();
    klm_hits_ = muid.N_layer_hit_brl() + muid.N_layer_hit_end();
    klm_chi2_ = muid.Chi_2();
  }
  const double & muonLikelihood(void) { return muon_likelihood_; }
  const double klmSignature(void) {
    return (klm_hits_ > 0 ?  klm_chi2_ / klm_hits_ : 0);
  }

  // Electron ID
  void electronLikelihood(const double &electron_likelihood) {
    electron_likelihood_ = electron_likelihood;
  }
  const double & electronLikelihood(void) { return electron_likelihood_; }

  // Kaon ID
  void kaonToPionLikelihood(const double &kaon_to_pi_likelihood) {
    kaon_to_pion_likelihood_ = kaon_to_pi_likelihood;
  }
  void kaonToProtonLikelihood(const double &kaon_to_proton_likelihood) {
    kaon_to_proton_likelihood_ = kaon_to_proton_likelihood;
  }
  const double & kaonToPionLikelihood(void) { return kaon_to_pion_likelihood_; }
  const double & kaonToProtonLikelihood(void) { return kaon_to_proton_likelihood_; }

  void svdHits(const Mdst_charged &mdst) {
    svd_r_electron_hits_ = mdst.trk().mhyp(0).nhits(3);
    svd_z_electron_hits_ = mdst.trk().mhyp(0).nhits(4);
    svd_r_muon_hits_ = mdst.trk().mhyp(1).nhits(3);
    svd_z_muon_hits_ = mdst.trk().mhyp(1).nhits(4);
    svd_r_kaon_hits_ = mdst.trk().mhyp(3).nhits(3);
    svd_z_kaon_hits_ = mdst.trk().mhyp(3).nhits(4);
  }
  void svdHits(const double &r_hits, const double &z_hits) {
    svd_r_hits_ = r_hits;
    svd_z_hits_ = z_hits;
  }
  const double svdRHits(const int &mass_hyp) {
    if (mass_hyp == 0) {
      return svd_r_electron_hits_;
    } else if (mass_hyp == 1) {
      return svd_r_muon_hits_;
    } else if (mass_hyp == 3) {
      return svd_r_kaon_hits_;
    } else {
      return 0;
    }
  }
  const double svdZHits(const int &mass_hyp) {
    if (mass_hyp == 0) {
      return svd_z_electron_hits_;
    } else if (mass_hyp == 1) {
      return svd_z_muon_hits_;
    } else if (mass_hyp == 3) {
      return svd_z_kaon_hits_;
    } else {
      return 0;
    }
  }
  const double svdRHits(void) { return svd_r_hits_; }
  const double svdZHits(void) { return svd_z_hits_; }

  // Track Parameters.
  void ipDeltaR(const double &ip_dr) { ip_dr_ = ip_dr; }
  void ipDeltaZ(const double &ip_dz) { ip_dr_ = ip_dz; }
  const double & ipDeltaR(void) { return ip_dr_; }
  const double & ipDeltaZ(void) { return ip_dz_; }

  // CM-Frame momentum.
  void pCm(const HepLorentzVector &lab_momentum, const Hep3Vector &boost) {
    p_cm_ = lab_momentum;
    p_cm_.boost(boost);
  }
  const HepLorentzVector & pCm(void) { return p_cm_; }

 private:
  int m_genHepevtLink;
  int m_genHepevtChecked;
  int m_decayMode;
  double   m_chisq;
  double   m_cl;
  unsigned m_ndf;

  HepLorentzVector p_cm_;
  double muon_likelihood_;
  double electron_likelihood_;
  double kaon_to_pion_likelihood_;
  double kaon_to_proton_likelihood_;
  double klm_chi2_;
  unsigned klm_hits_;
  double svd_r_hits_;
  double svd_z_hits_;
  double svd_r_electron_hits_;
  double svd_z_electron_hits_;
  double svd_r_muon_hits_;
  double svd_z_muon_hits_;
  double svd_r_kaon_hits_;
  double svd_z_kaon_hits_;
  double ip_dr_;
  double ip_dz_;
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
