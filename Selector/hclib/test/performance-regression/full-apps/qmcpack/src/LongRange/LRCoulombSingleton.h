//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file LRCoulombSingleton.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_LRCOULOMBSINGLETON_H
#define QMCPLUSPLUS_LRCOULOMBSINGLETON_H

#include <config.h>
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRHandlerSRCoulomb.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimLinearSpline.h"

namespace qmcplusplus
{

struct LRCoulombSingleton
{

  typedef OHMMS_PRECISION                                    RealType;
  typedef LRHandlerBase                                      LRHandlerType;
  typedef LinearGrid<RealType>                               GridType;
  //    typedef OneDimLinearSpline<RealType>                 RadFunctorType;
  typedef OneDimCubicSpline<RealType>                       RadFunctorType;

  static LRHandlerType* CoulombHandler;
  static LRHandlerType* CoulombDerivHandler;
  static LRHandlerType* getHandler(ParticleSet& ref);
  static LRHandlerType* getDerivHandler(ParticleSet& ref);
  /** create a linear spline function
   * @param aLR LRHandler
   * @param rcut cutoff radius
   * @param agrid pointer to a grid
   * @return a RadFunctorType
   *
   * The spline function is the short-range term after breaking up
   * \f$r V_{S} = r \times (V(r)-V_{L})\f$
   */
  static RadFunctorType* createSpline4RbyVs(LRHandlerType* aLR, RealType rcut,
      GridType* agrid=0);
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: j.k.rofling@gmail.com $
 * $Revision: 6229 $   $Date: 2014-02-17 13:34:17 -0500 (Mon, 17 Feb 2014) $
 * $Id: LRCoulombSingleton.h 6229 2014-02-17 18:34:17Z j.k.rofling@gmail.com $
 ***************************************************************************/
