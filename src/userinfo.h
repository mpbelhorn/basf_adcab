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
#include "CLHEP/Matrix/Matrix.h"

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

  /// Input vertex Chi^2.
  void chisq(const double &v) { m_chisq = v; }

  /// Vertex fit confidence level.
  void cl(const double &v)    { m_cl = v; }

  /// Vertex fit degrees of freedom.
  void ndf(const unsigned &v) { m_ndf = v; }

  /// Output.
  const double   & chisq(void) const { return m_chisq; }
  const double   & cl(void)    const { return m_cl; }
  const unsigned & ndf(void)   const { return m_ndf; }

 private:
  int m_genHepevtLink;
  int m_genHepevtChecked;
  int m_decayMode;
  double   m_chisq;
  double   m_cl;
  unsigned m_ndf;
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
