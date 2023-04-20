//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
#ifndef OHMMS_DUPLICATE_SUPERCELL_H
#define OHMMS_DUPLICATE_SUPERCELL_H

// forward declaration
template<class T, unsigned D> class CrystalLattice;

// generic class declaration
template<class CL>
struct DuplicateSuperCell { };

// specialized class for CrystalLattice<T,D>
template<class T, unsigned D>
struct DuplicateSuperCell<CrystalLattice<T,D> >
{
  typedef typename CrystalLattice<T,D>::SingleParticlePos_t SingleParticlePos_t;
  inline static void apply(const CrystalLattice<T,D>& , vector<SingleParticlePos_t>& nc, T rmax)
  {
  }
};

// specialized class for CrystalLattice<T,1>
template<class T>
struct DuplicateSuperCell<CrystalLattice<T,1> >
{

  typedef typename CrystalLattice<T,1>::SingleParticlePos_t SingleParticlePos_t;

  inline static void apply(
    const CrystalLattice<T,1>& lat, vector<SingleParticlePos_t>& nc, T rmax)
  {
    if(nc.empty())
    {
      T xmin = lat.R(0,0); // only need the size
      T x = xmin;
      int nx=0;
      while(x < rmax)
      {
        x+= xmin;
        nx++;
      }
      // to ensure that the first cell is always included
      if(nx == 0)
        nx = 1;
      // box in itself is always added to nextcell
      // nextcell is in a supercell unit
      nc.push_back(SingleParticlePos_t(0.0e0));
      for(int ix =1; ix <= nx; ix++)
      {
        nc.push_back(SingleParticlePos_t(-static_cast<T>(ix)));
        nc.push_back(SingleParticlePos_t( static_cast<T>(ix)));
      }
    }
  }
};
// specialized class for CrystalLattice<T,2>
template<class T>
struct DuplicateSuperCell<CrystalLattice<T,2> >
{

  typedef typename CrystalLattice<T,2>::SingleParticlePos_t SingleParticlePos_t;
  inline static void apply(const CrystalLattice<T,2>& lat, vector<SingleParticlePos_t>& nc, T rmax)
  {
    if(nc.empty())
    {
      // unit vectors
      SingleParticlePos_t u0(1.0,0.0);
      SingleParticlePos_t u1(0.0,1.0);
      T xmin = sqrt(lat.Dot(u0,u0));
      T ymin = sqrt(lat.Dot(u1,u1));
      T x= xmin;
      int nx=0;
      int ny =0;
      int nz = 0;
      T Rmax = rmax*sqrt(2.0e0);
      while(x < Rmax)
      {
        x+= xmin;
        nx++;
      }
      x = ymin;
      while(x < Rmax)
      {
        x+= ymin;
        ny++;
      }
      // to ensure that the first cell is always included
      if(nx == 0)
        nx = 1;
      if(ny == 0)
        ny = 1;
      // box in itself is always added to nextcell
      // nextcell is in a supercell unit
      nc.push_back(SingleParticlePos_t(0.0e0, 0.0e0));
      for(int ix = -nx; ix <= nx; ix++)
      {
        T dx = static_cast<T>(ix);
        for(int iy = -ny; iy <= ny; iy++)
        {
          T dy = static_cast<T>(iy);
          if(ix != 0 || iy != 0)
          {
            nc.push_back(dx*u0+dy*u1);
          }
        }
      }
    }
  }
};

// specialized class for CrystalLattice<T,3>
template<class T>
struct DuplicateSuperCell<CrystalLattice<T,3> >
{

  typedef typename CrystalLattice<T,3>::SingleParticlePos_t SingleParticlePos_t;
  inline static void apply(const CrystalLattice<T,3>& lat, vector<SingleParticlePos_t>& nc, T rmax)
  {
    if(nc.empty())
    {
      // unit vectors
      SingleParticlePos_t u0(1.0,0.0,0.0);
      SingleParticlePos_t u1(0.0,1.0,0.0);
      SingleParticlePos_t u2(0.0,0.0,1.0);
      T xmin = sqrt(lat.Dot(u0,u0));
      T ymin = sqrt(lat.Dot(u1,u1));
      T zmin = sqrt(lat.Dot(u2,u2));
      T x= xmin;
      int nx=0;
      int ny =0;
      int nz = 0;
      T Rmax = rmax*sqrt(2.0e0);
      while(x < Rmax)
      {
        x+= xmin;
        nx++;
      }
      x = ymin;
      while(x < Rmax)
      {
        x+= ymin;
        ny++;
      }
      x = zmin;
      while(x < Rmax)
      {
        x+= zmin;
        nz++;
      }
      //nx++; ny++; nz++;
      // to ensure that the first cell is always included
      if(nx == 0)
        nx = 1;
      if(ny == 0)
        ny = 1;
      if(nz == 0)
        nz = 1;
      if(!lat.BoxBConds[0])
        nx = 0;
      if(!lat.BoxBConds[1])
        ny = 0;
      if(!lat.BoxBConds[2])
        nz = 0;
      // box in itself is always added to nextcell
      // nextcell is in a supercell unit
      nc.push_back(SingleParticlePos_t(0.0e0, 0.0e0, 0.0e0));
      for(int ix = -nx; ix <= nx; ix++)
      {
        T dx = static_cast<T>(ix);
        for(int iy = -ny; iy <= ny; iy++)
        {
          T dy = static_cast<T>(iy);
          for(int iz = -nz; iz <= nz; iz++)
          {
            T dz = static_cast<T>(iz);
            if(ix != 0 || iy != 0 || iz != 0)
            {
              nc.push_back(dx*u0+dy*u1+dz*u2);
            }
          }
        }
      }
    }
  }
};

#endif // OHMMS_DUPLICATE_SUPERCELL_H

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: DuplicateSuperCell.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
