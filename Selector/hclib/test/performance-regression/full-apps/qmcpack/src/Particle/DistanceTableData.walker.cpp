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
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
using namespace qmcplusplus;

void SymmetricDTD::reset(int m, int nactive)
{
  if(m != N[SourceIndex] || nactive != N[WalkerIndex])
  {
    N[SourceIndex]=m;
    N[VisitorIndex]=m;
    resize(m*(m-1)/2,nactive);
    M.resize(m+1);
    J.resize(m*(m-1));
    M[0] = 0;
    int ij = 0;
    for(int i=0; i<m; i++)
    {
      for(int j=i+1; j<m; j++, ij++)
      {
        J[ij] = j;
      }
      M[i+1] = ij;
    }
  }
}

///evaluate the Distance Table using a set of Particle Positions
void SymmetricDTD::evaluate(const PosVector_t& a, int visitors, int copies)
{
  ///number of columns
  reset(visitors,copies);
  int ij=0;
  for(int i=0; i<visitors; i++)
  {
    for(int j=i+1; j<visitors; j++)
    {
      int i0 = i*copies;
      int j0 = j*copies;
      for(int iw=0; iw<copies; iw++, ij++, i0++, j0++)
      {
        SPPosition_t drij = a(j0)-a(i0);
        value_type sep = sqrt(dot(drij,drij));
        r(ij) = sep;
        rinv(ij) = 1.0/sep;
        dr(ij) = drij;
      }
    }
  }
}

void AsymmetricDTD::reset(int n1, int n2, int nactive)
{
  if( n1!=N[SourceIndex] || n2 != N[VisitorIndex] || nactive != N[WalkerIndex])
  {
    N[SourceIndex] = n1;
    N[VisitorIndex] = n2;
    int m = n1*n2;
    if(m)
    {
      resize(m,nactive);
      M.resize(n1+1);
      J.resize(m);
      M[0] = 0;
      int ij = 0;
      for(int i=0; i<n1; i++)
      {
        for(int j=0; j<n2; j++, ij++)
        {
          J[ij] = j;
        }
        M[i+1] = M[i]+n2;
      }
    }
  }
}

void AsymmetricDTD::evaluate(const PosVector_t& a, int visitors, int copies)
{
  ///number of columns
  int ns = Origin.getTotalNum();
  reset(ns,visitors,copies);
  int ij = 0;
  for(int i=0; i<ns; i++)
  {
    SPPosition_t r0 = Origin.R(i);
    int j0=0;
    for(int j=0; j<visitors; j++)
    {
      for(int ia=0; ia<copies; ia++, ij++, j0++)
      {
        SPPosition_t drij = a(j0)-r0;
        value_type sep = sqrt(dot(drij,drij));
        r(ij) = sep;
        rinv(ij) = 1.0/sep;
        dr(ij) = drij;
      }
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: DistanceTableData.walker.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
