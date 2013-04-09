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
      m_ndf(0)
{
  // Empty function.
}

UserInfo::UserInfo(unsigned ch)
    : m_genHepevtLink(0),
      m_genHepevtChecked(0),
      m_decayMode(0),
      m_chisq(-1.),
      m_cl(-1.),
      m_ndf(0)
{
  // Empty function.
}

UserInfo::UserInfo(const UserInfo &x)
    : m_genHepevtLink(x.m_genHepevtLink),
      m_genHepevtChecked(x.m_genHepevtChecked),
      m_decayMode(x.m_decayMode),
      m_chisq(x.m_chisq),
      m_cl(x.m_chisq),
      m_ndf(x.m_ndf)
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
  return *this;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
