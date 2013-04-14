//______________________________________________________________________________
// Filename: VertexFit.cpp
// Version 0.1
// Date: 2013.04.09
// Original Date: 2013.04.09
// Author: M.P. Belhorn
// Description: Functions to perform kinematic vertex fits.
//______________________________________________________________________________

#include "VertexFit.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

int
fitPhiVertex(Particle &p)
{
  setUserInfo(p);
  kvertexfitter kvf;
  // Check the pointer goes somewhere.
  addTrack2fit(kvf, p.relation().child(0));
  addTrack2fit(kvf, p.relation().child(1));
  unsigned err = kvf.fit();
  UserInfo &info = dynamic_cast<UserInfo&>(p.userInfo());
  if (err == 0) {
    info.cl(kvf.cl());
    info.chisq(kvf.chisq());
    info.ndf(kvf.dgf());
    makeMother(kvf, p); // calculate "D0" using the fitting result.
    return 0;
  } else {
    info.cl(-1.0);
    HepPoint3D vtx(999.,999.,999.);
    HepSymMatrix errVtx(3,0);
    p.momentum().decayVertex(vtx, errVtx);
    return 1;
  }
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
