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

class UserInfo : public ParticleUserInfo
{
public:
	/// Default constructor
	UserInfo();

	UserInfo(unsigned);

	/// Copy constructor
	UserInfo(const UserInfo &);

	/// Destructor
	virtual ~UserInfo();

	/// constructs self object.
	UserInfo * clone(void) const;

	/// Copy operator
	UserInfo & operator = (const UserInfo &);

public:
	const int & decayMode(void) const { return m_decayMode; }
	void decayMode(const int &v) { m_decayMode = v; }

	const int & genHepevtLink(void) const { return m_genHepevtLink; }
	void genHepevtLink(const int &v) { m_genHepevtLink = v; }

	const int & genHepevtChecked(void) const { return m_genHepevtChecked; }
	void genHepevtChecked(const int &v) { m_genHepevtChecked = v; }

private:
	int m_genHepevtLink;
	int m_genHepevtChecked;
	int m_decayMode;


};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
