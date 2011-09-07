//
//******************************************************************************
// Filename: EntryTypes.h
// Version: 2010.11.03.A
// Author: M.P. Belhorn
// Original Date: 2010.07.20
// Description: Data structure containing Adcab ntuple row identifiers.
//******************************************************************************

#ifndef ENTRYTYPES_H
#define ENTRYTYPES_H

#if defined( BELLE_NAMESPACE )
namespace Belle {
#endif
  
// Particle Selection Cuts.
struct EntryTypes {
  static const int event  = 1;
  static const int lepton = 2;
  static const int kaon   = 3;
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif
