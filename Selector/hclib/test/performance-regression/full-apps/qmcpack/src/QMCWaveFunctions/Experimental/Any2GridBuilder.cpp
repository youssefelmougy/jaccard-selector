//////////////////////////////////////////////////////////////////
// (c) Copyright 2003 by Jeongnim Kim and Jordan Vincent
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
#include "Utilities/OhmmsInfo.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Numerics/GaussianBasisSet.h"
#include "Numerics/SlaterBasisSet.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.h"
using namespace std;
namespace qmcplusplus
{

Any2GridBuilder::Any2GridBuilder(xmlNodePtr cur):
  Normalized(true),m_rcut(-1.0)
{
  if(cur != NULL)
  {
    putCommon(cur);
  }
}

bool Any2GridBuilder::putCommon(xmlNodePtr cur)
{
  const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
  if(a)
  {
    if(xmlStrEqual(a,(const xmlChar*)"no"))
      Normalized=false;
  }
  return true;
}

/** Add a new Slater Type Orbital with quantum numbers \f$(n,l,m,s)\f$
 * \param cur  the current xmlNode to be processed
 * \param nlms a vector containing the quantum numbers \f$(n,l,m,s)\f$
 * \return true is succeeds
 *
 This function puts the STO on a logarithmic grid and calculates the boundary
 conditions for the 1D Cubic Spline.  The derivates at the endpoint
 are assumed to be all zero.  Note: for the radial orbital we use
 \f[ f(r) = \frac{R(r)}{r^l}, \f] where \f$ R(r) \f$ is the usual
 radial orbital and \f$ l \f$ is the angular momentum.
*/
bool
Any2GridBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  const xmlChar *tptr = xmlGetProp(cur,(const xmlChar*)"type");
  string radtype("Gaussian");
  if(tptr)
    radtype = (const char*)tptr;
  tptr = xmlGetProp(cur,(const xmlChar*)"rmax");
  if(tptr)
    m_rcut = atof((const char*)tptr);
  int lastRnl = m_orbitals->Rnl.size();
  app_log() << "    basisGroup " << radtype << endl;
  m_nlms = nlms;
  if(radtype == "Gaussian")
  {
    addGaussian(cur);
  }
  else
    if(radtype == "Slater")
    {
      addSlater(cur);
    }
    else
      if(radtype == "Pade")
      {
        app_error() << "  Any2GridBuilder::addPade is disabled." << endl;
        abort();
        //addPade(cur);
      }
  if(lastRnl && m_orbitals->Rnl.size()> lastRnl)
  {
    //LOGMSG("\tSetting GridManager of " << lastRnl << " radial orbital to false")
    m_orbitals->Rnl[lastRnl]->setGridManager(false);
  }
  return true;
}

void Any2GridBuilder::addGaussian(xmlNodePtr cur)
{
  int L= m_nlms[1];
  GaussianCombo<RealType> gset(L,Normalized);
  gset.putBasisGroup(cur);
  GridType* agrid = m_orbitals->Grids[0];
  RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);
  if(m_rcut<0)
    m_rcut = agrid->rmax();
  Transform2GridFunctor<GaussianCombo<RealType>,RadialOrbitalType> transform(gset, *radorb);
  transform.generate(agrid->rmin(),m_rcut,agrid->size());
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(m_nlms);
}

void Any2GridBuilder::addSlater(xmlNodePtr cur)
{
  ////pointer to the grid
  GridType* agrid = m_orbitals->Grids[0];
  RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);
  SlaterCombo<RealType> sto(m_nlms[1],Normalized);
  sto.putBasisGroup(cur);
  //spline the slater type orbital
  Transform2GridFunctor<SlaterCombo<RealType>,RadialOrbitalType> transform(sto, *radorb);
  transform.generate(agrid->rmin(), agrid->rmax(),agrid->size());
  //transform.generate(agrid->rmax());
  //add the radial orbital to the list
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(m_nlms);
}

void Any2GridBuilder::addPade(xmlNodePtr cur)
{
  //Who is using this????
  //GridType* agrid = m_orbitals->Grids[0];
  //RadialOrbitalType *radorb = new OneDimCubicSpline<RealType>(agrid);
  //PadeOrbital<RealType> pade;
  //pade.putBasisGroup(cur);
  ////spline the slater type orbital
  //Transform2GridFunctor<PadeOrbital<RealType>,RadialOrbitalType> transform(pade, *radorb);
  //if(pade.rcut>0)
  //  transform.generate(pade.rcut);
  //else
  //  transform.generate(agrid->rmax());
  ////add the radial orbital to the list
  //m_orbitals->Rnl.push_back(radorb);
  //m_orbitals->RnlID.push_back(m_nlms);
}
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: Any2GridBuilder.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
