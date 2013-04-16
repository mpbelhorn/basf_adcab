//______________________________________________________________________________
// Filename: Userinfo.cc
// Author: Unknown.
// Modified: M.P.Belhorn (matt.belhorn@gmail.com).
// Description: Class for custom particle information. Used with A. Zupanc's
//     MC truthtable fuctions.
// _____________________________________________________________________________

#include "belle.h"
#include "userinfo.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// For Interface to Set UserInfo Class
void
setUserInfo(Particle &p, unsigned ch)
{
  if(! &p.userInfo()) {
    p.userInfo(UserInfo(ch));
  }
}

void
setUserInfo(std::vector<Particle> &p, unsigned ch)
{
  for (unsigned i = 0; i < p.size(); ++i) {
    setUserInfo(p[i], ch);
  }
}

// UserInfo Class
UserInfo::UserInfo()
    : m_genHepevtLink(0),
      m_genHepevtChecked(0),
      m_decayMode(0),
      m_chisq(-1.),
      m_cl(-1.),
      m_ndf(0),
      p_cm_(HepLorentzVector(0,0,0,0)),
      muon_likelihood_(0),
      electron_likelihood_(0),
      kaon_to_pion_likelihood_(-1),
      kaon_to_proton_likelihood_(-1),
      klm_chi2_(-1),
      klm_hits_(0),
      svd_r_hits_(-1),
      svd_z_hits_(-1),
      svd_r_electron_hits_(-1),
      svd_z_electron_hits_(-1),
      svd_r_muon_hits_(-1),
      svd_z_muon_hits_(-1),
      svd_r_kaon_hits_(-1),
      svd_z_kaon_hits_(-1),
      ip_dr_(-9999),
      ip_dz_(-9999)
{
  // Empty function.
}

UserInfo::UserInfo(unsigned ch)
    : m_genHepevtLink(0),
      m_genHepevtChecked(0),
      m_decayMode(0),
      m_chisq(-1.),
      m_cl(-1.),
      m_ndf(0),
      p_cm_(HepLorentzVector(0,0,0,0)),
      muon_likelihood_(0),
      electron_likelihood_(0),
      kaon_to_pion_likelihood_(-1),
      kaon_to_proton_likelihood_(-1),
      klm_chi2_(-1),
      klm_hits_(0),
      svd_r_hits_(-1),
      svd_z_hits_(-1),
      svd_r_electron_hits_(-1),
      svd_z_electron_hits_(-1),
      svd_r_muon_hits_(-1),
      svd_z_muon_hits_(-1),
      svd_r_kaon_hits_(-1),
      svd_z_kaon_hits_(-1),
      ip_dr_(-9999),
      ip_dz_(-9999)
{
  // Empty function.
}

UserInfo::UserInfo(const UserInfo &x)
    : m_genHepevtLink(x.m_genHepevtLink),
      m_genHepevtChecked(x.m_genHepevtChecked),
      m_decayMode(x.m_decayMode),
      m_chisq(x.m_chisq),
      m_cl(x.m_chisq),
      m_ndf(x.m_ndf),
      p_cm_(x.p_cm_),
      muon_likelihood_(x.muon_likelihood_),
      electron_likelihood_(x.electron_likelihood_),
      kaon_to_pion_likelihood_(x.kaon_to_pion_likelihood_),
      kaon_to_proton_likelihood_(x.kaon_to_proton_likelihood_),
      klm_chi2_(x.klm_chi2_),
      klm_hits_(x.klm_hits_),
      svd_r_hits_(x.svd_r_hits_),
      svd_z_hits_(x.svd_z_hits_),
      svd_r_electron_hits_(x.svd_r_electron_hits_),
      svd_z_electron_hits_(x.svd_z_electron_hits_),
      svd_r_muon_hits_(x.svd_r_muon_hits_),
      svd_z_muon_hits_(x.svd_z_muon_hits_),
      svd_r_kaon_hits_(x.svd_r_kaon_hits_),
      svd_z_kaon_hits_(x.svd_z_kaon_hits_),
      ip_dr_(x.ip_dr_),
      ip_dz_(x.ip_dz_)
{
  // Empty function.
}

UserInfo::~UserInfo()
{
  // Empty function.
}

UserInfo*
UserInfo::clone(void) const
{
  UserInfo *x = new UserInfo(*this);
  return x;
}

UserInfo &
UserInfo::operator = (const UserInfo &x)
{
  if(this == &x) return *this;
  m_genHepevtLink = x.m_genHepevtLink;
  m_genHepevtChecked = x.m_genHepevtChecked;
  m_decayMode = x.m_decayMode;
  m_chisq = x.m_chisq;
  m_cl = x.m_cl;
  m_ndf = x.m_ndf;
  p_cm_ = x.p_cm_;
  muon_likelihood_ = x.muon_likelihood_;
  electron_likelihood_ = x.electron_likelihood_;
  kaon_to_pion_likelihood_ = x.kaon_to_pion_likelihood_;
  kaon_to_proton_likelihood_ = x.kaon_to_proton_likelihood_;
  klm_chi2_ = x.klm_chi2_;
  klm_hits_ = x.klm_hits_;
  svd_r_hits_ = x.svd_r_hits_;
  svd_z_hits_ = x.svd_z_hits_;
  svd_r_electron_hits_ = x.svd_r_electron_hits_;
  svd_z_electron_hits_ = x.svd_z_electron_hits_;
  svd_r_muon_hits_ = x.svd_r_muon_hits_;
  svd_z_muon_hits_ = x.svd_z_muon_hits_;
  svd_r_kaon_hits_ = x.svd_r_kaon_hits_;
  svd_z_kaon_hits_ = x.svd_z_kaon_hits_;
  ip_dr_ = x.ip_dr_;
  ip_dz_ = x.ip_dz_;
  return *this;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
