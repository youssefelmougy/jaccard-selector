//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SLATERTYPEORBITAL_BUILDER_H
#define QMCPLUSPLUS_SLATERTYPEORBITAL_BUILDER_H

#include "Configuration.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"
#include "Numerics/SlaterBasisSet.h"

namespace qmcplusplus
{

class STOBuilder: public QMCTraits
{

public:

  typedef SlaterCombo<RealType>                      RadialOrbitalType;
  typedef SphericalBasisSet<RadialOrbitalType>     CenteredOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///the species
  std::string m_species;
  ///constructor
  STOBuilder(xmlNodePtr cur=NULL);

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
  {
    m_orbitals = oset;
    m_species = acenter;
  }

  bool addGrid(xmlNodePtr cur)
  {
    return true;
  }

  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

  bool putCommon(xmlNodePtr cur);

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: STOBuilder.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
