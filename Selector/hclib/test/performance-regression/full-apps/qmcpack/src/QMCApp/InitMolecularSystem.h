//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file InitMolecularSystem.h
 * @brief Declaration of InitMolecularSystem
 */
#ifndef QMCPLUSPLUS_INITMOLECULARSYSTEM_H
#define QMCPLUSPLUS_INITMOLECULARSYSTEM_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>

namespace qmcplusplus
{

class ParticleSet;
class ParticleSetPool;

/* Engine to initialize the initial electronic structure for a molecular system
 */
class InitMolecularSystem : public OhmmsElementBase
{

public:

  InitMolecularSystem(ParticleSetPool* pset, const char* aname = "mosystem");

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();

  /** initialize els for an atom
   */
  void initAtom(ParticleSet* ions, ParticleSet* els);
  /** initialize els position for a molecule
   *
   * Use the valence of each ionic species on a sphere
   */
  void initMolecule(ParticleSet* ions, ParticleSet* els);
  /** initialize els for the systems with a mixed boundary
   *
   * Use the bound of the ionic systems and uniform random positions within a reduced box
   */
  void initWithVolume(ParticleSet* ions, ParticleSet* els);

private:

  /** pointer to ParticleSetPool
   *
   * QMCHamiltonian needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   * Any number of ParticleSet can be used to describe
   * a QMCHamiltonian.
   */
  ParticleSetPool* ptclPool;

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 5848 $   $Date: 2013-05-15 09:33:24 -0400 (Wed, 15 May 2013) $
 * $Id: InitMolecularSystem.h 5848 2013-05-15 13:33:24Z jnkim $
 ***************************************************************************/
