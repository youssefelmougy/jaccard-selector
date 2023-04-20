//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_GAUSSIANORBITAL_MOLECULARORBITALS_H
#define QMCPLUSPLUS_GAUSSIANORBITAL_MOLECULARORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBasis.h"
#include "Numerics/GaussianBasisSet.h"

namespace qmcplusplus
{

/**Class to add a set of Gaussian Type atomic orbital basis functions
 *to the collection of basis functions.
 *
 @brief Example of a ROT (Radial Orbital Type)
*/
class GTOMolecularOrbitals: public OrbitalBuilderBase
{
public:

  typedef GaussianCombo<RealType>                    RadialOrbitalType;
  typedef SphericalOrbitalSet<RadialOrbitalType>     CenteredOrbitalType;
  typedef MolecularOrbitalBasis<CenteredOrbitalType> BasisSetType;

  ///constructor
  GTOMolecularOrbitals(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions);

  ///implement vritual function
  bool put(xmlNodePtr cur);

  ///returns a BasisSet
  BasisSetType* addBasisSet(xmlNodePtr cur);

private:

  enum {DONOT_EXPAND=0, GAUSSIAN_EXPAND=1, NATURAL_EXPAND};

  bool Normalized;
  ParticleSet& IonSys;
  BasisSetType*      BasisSet;
  DistanceTableData* d_table;
  map<string,int>    RnlID;
  map<string,int>    CenterID;
  ///map for (n,l,m,s) to its quantum number index
  map<string,int> nlms_id;

  int expandYlm(const string& rnl, const QuantumNumberType& nlms,
                int num, CenteredOrbitalType* aos, xmlNodePtr cur1,
                int expandlm);

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: GTOMolecularOrbitals.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
