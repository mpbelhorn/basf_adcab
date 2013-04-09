#ifndef VERTEXFIT_H
#define VERTEXFIT_H

#include "belle.h"
#include "particle/Particle.h"
#include "particle/utility.h"
#include "userinfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

int fitPhiVertex(Particle &p);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
