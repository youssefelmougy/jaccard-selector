//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H
#define QMCPLUSPLUS_NUMERICALGRIDORBITALBUILDER_H

#include "Configuration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"

namespace qmcplusplus
{

struct NGOrbital: public OptimizableFunctorBase
{
  typedef real_type                    value_type;
  typedef real_type                    point_type;
  typedef OneDimGridBase<real_type>    grid_type;
#if QMC_BUILD_LEVEL>2
  typedef OneDimQuinticSpline<real_type> functor_type;
#else
  typedef OneDimCubicSpline<real_type> functor_type;
#endif
  functor_type myFunc;
  real_type Y, dY, d2Y, d3Y;

  NGOrbital(grid_type* agrid):myFunc(agrid) { }

  template<typename VV>
  NGOrbital(grid_type* agrid, const VV& nv):myFunc(agrid,nv) { }

  void checkInVariables(opt_variables_type& active) {}
  void checkOutVariables(const opt_variables_type& active) {}
  void resetParameters(const opt_variables_type& active) {}
  void reset() {}
  inline real_type f(real_type r)
  {
    return myFunc.f(r);
  }
  inline real_type df(real_type r)
  {
    return myFunc.df(r);
  }
  bool put(xmlNodePtr cur)
  {
    return true;
  }
  OptimizableFunctorBase* makeClone() const;

  inline real_type evaluate(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r);
  }
  inline value_type evaluateAll(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r,dY,d2Y);
  }

  inline value_type evaluateWithThirdDeriv(real_type r, real_type rinv)
  {
    return Y=myFunc.splint(r,dY,d2Y,d3Y);
  }

  inline value_type operator()(int i) const
  {
    return myFunc(i);
  }
  inline value_type& operator()(int i)
  {
    return myFunc(i);
  }
  inline grid_type& grid()
  {
    return myFunc.grid();
  }
  inline void setGridManager(bool willmanage)
  {
    myFunc.setGridManager(willmanage);
  }

  inline void spline(int imin, value_type yp1, int imax, value_type ypn)
  {
    myFunc.spline(imin,yp1,imax,ypn);
  }

  inline void resize(int n)
  {
    myFunc.resize(n);
  }
};

/**Class to convert SlaterTypeOrbital to a radial orbital on a log grid.
 *
 * For a center,
 *   - only one grid is used
 *   - any number of radial orbitals
 */
class NGOBuilder: public QMCTraits
{

public:
  //typedef OneDimGridBase<RealType>                        GridType;
  //typedef OneDimGridFunctor<RealType>                     RadialOrbitalType;
  typedef NGOrbital                                     RadialOrbitalType;
  typedef NGOrbital::grid_type                          GridType;
  typedef SphericalBasisSet<RadialOrbitalType,GridType> CenteredOrbitalType;

  ///true, if the RadialOrbitalType is normalized
  bool Normalized;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;
  ///input grid in case transform is needed
  GridType *input_grid;
  ///maximum cutoff
  RealType m_rcut;
  ///the quantum number of this node
  QuantumNumberType m_nlms;
  ///the species
  std::string m_species;
  ///type of input function
  std::string m_infunctype;

  ///constructor
  NGOBuilder(xmlNodePtr cur=NULL);
  ///destructor
  ~NGOBuilder();

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter);

  ///add a grid
  bool addGrid(xmlNodePtr cur);

  /** add a radial functor
   * @param cur xml element
   * @param nlms quantum number
   */
  bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms);

  /** put common element
   * @param cur xml element
   */
  bool putCommon(xmlNodePtr cur);

private:
  void addGaussian(xmlNodePtr cur);
  void addSlater(xmlNodePtr cur);
  void addNumerical(xmlNodePtr cur, const string& dsname);
  void addPade(xmlNodePtr cur);
  hid_t m_fileid;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: NGOBuilder.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
