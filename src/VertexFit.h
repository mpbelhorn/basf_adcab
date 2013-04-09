#ifndef VERTEXFIT_H
#define VERTEXFIT_H

#include "belle.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "userinfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

void setMCtruth(std::vector<Particle> &plist);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
