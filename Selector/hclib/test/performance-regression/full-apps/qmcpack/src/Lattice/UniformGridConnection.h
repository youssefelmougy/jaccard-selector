//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_GRID_CONNECTION_H
#define OHMMS_GRID_CONNECTION_H

#include <vector>
#include <map>
#include <iostream>
using namespace std;
#include "Lattice/CrystalLattice.h"

/* \class GridConnection
   \brief templated over scalar type, dimension and a key to
   make connection list, similar to the nearest neighbor lists
*/

template<class T, unsigned D>
struct UniformGridConnection {};

// specialization
template<class T>
struct UniformGridConnection<T,3>
{

  typedef CrystalLattice<T,3> ParticleLayout_t;
  typedef typename ParticleLayout_t::SingleParticlePos_t   SingleParticlePos_t;

  vector<int> M;
  vector<int> ID;
  vector<SingleParticlePos_t> BC;

  UniformGridConnection() { }

  //!< return the first mesh connected to the ig grid
  inline int first(int ig) const
  {
    return M[ig];
  }

  //!< return the last mesh connected to the ig grid
  inline int last(int ig) const
  {
    return M[ig+1];
  }

  //!< return the mesh index
  inline int id(int j) const
  {
    return ID[j];
  }

  //!< return the correction vector according to the boundary conditions
  inline SingleParticlePos_t bc(int j) const
  {
    return BC[j];
  }

  //!< return the maximum number of connected cells
  inline int connect(const CrystalLattice<T,3>& Lattice, T rmax, int glevel)
  {
    const typename ParticleLayout_t::PtclGrid_t&  basegrid = *(Lattice.dGrid[glevel]);
    SingleParticlePos_t u0(basegrid.Delta[0],0.0,0.0);
    SingleParticlePos_t u1(0.0,basegrid.Delta[1],0.0);
    SingleParticlePos_t u2(0.0,0.0,basegrid.Delta[2]);
    T RmaxSq = rmax*rmax;
    ///calculate the extend of linked cells
    int nx = static_cast<int>(sqrt(RmaxSq/Lattice.Dot(u0,u0)))+1;
    int ny = static_cast<int>(sqrt(RmaxSq/Lattice.Dot(u1,u1)))+1;
    int nz = static_cast<int>(sqrt(RmaxSq/Lattice.Dot(u2,u2)))+1;
    M.resize(basegrid.getTotalNum()+1);
    int ntot = basegrid.getTotalNum()*nx*ny*nz;
    if(ID.size() < ntot)
      ID.reserve(ntot);
    if(BC.size() < ntot)
      BC.reserve(ntot);
    M[0] = 0;
    int maxnc = 0;
    SingleParticlePos_t dx(basegrid.Delta[0], basegrid.Delta[1], basegrid.Delta[2]);
    ///for each grid, search for connected grids
    for(int ig=0; ig<basegrid.getTotalNum(); ig++)
    {
      SingleParticlePos_t org(basegrid.Node[ig].Ri[0]+0.5*dx[0],
                              basegrid.Node[ig].Ri[1]+0.5*dx[1],
                              basegrid.Node[ig].Ri[2]+0.5*dx[2]),d;
      T x,y,z;
      int nconnected = 0;
      for(int ix=-nx; ix<=nx; ix++)
      {
        d[0] = x = org[0]+double(ix)*dx[0];
        if(Lattice.BoxBConds[0])
        {
          x = fmod(d[0],1.0);
          if(x<0)
            x += 1.0;
        }
        if(x<0 || x>=1)
          continue;
        d[0] -= x;
        for(int jx=-ny; jx<=ny; jx++)
        {
          d[1] = y = org[1]+double(jx)*dx[1];
          if(Lattice.BoxBConds[1])
          {
            y = fmod(d[1],1.0);
            if(y<0)
              y += 1.0;
          }
          if(y<0 || y>=1)
            continue;
          d[1] -= y;
          for(int kx=-nz; kx<=nz; kx++)
          {
            ///exclude itself
            if(ix == 0 && jx == 0 && kx == 0)
              continue;
            d[2] = z = org[2]+double(kx)*dx[2];
            if(Lattice.BoxBConds[2])
            {
              z = fmod(d[2],1.0);
              if(z<0)
                z += 1.0;
            }
            if(z<0 || z>=1)
              continue;
            d[2] -= z;
            ID.push_back(basegrid.loc(x,y,z));
            BC.push_back(Lattice.toCart(d));
            nconnected++;
          }
        }
      }
      M[ig+1] = M[ig]+nconnected;
      maxnc = max(maxnc,nconnected);
    }
//     //for each grid, search for connected grid
//     for(int ig=0; ig<basegrid.getTotalNum(); ig++) {
//       cout << ig << " has " << M[ig+1]-M[ig] << endl;
//       for(int j=M[ig]; j<M[ig+1]; j++) {
// 	cout << ID[j] << " " << BC[j] << endl;
//       }
//     }
    return maxnc; // return the maxmimum number of connected cells
  }
  /*** Extremely inefficient: using the new method above
  inline int connect(const CrystalLattice<T,3>& Lattice, T rmax, int glevel) {

    T RmaxSq = rmax*rmax;

    SingleParticlePos_t u0(1.0,0.0,0.0);
    SingleParticlePos_t u1(0.0,1.0,0.0);
    SingleParticlePos_t u2(0.0,0.0,1.0);

    int maxx1 = static_cast<int>(sqrt( RmaxSq/Lattice.Dot(u0,u0))) + 1;
    int maxx2 = static_cast<int>(sqrt( RmaxSq/Lattice.Dot(u1,u1))) + 1;
    int maxx3 = static_cast<int>(sqrt( RmaxSq/Lattice.Dot(u2,u2))) + 1;

    if(!Lattice.BoxBConds[0]) maxx1 = 0;
    if(!Lattice.BoxBConds[1]) maxx2 = 0;
    if(!Lattice.BoxBConds[2]) maxx3 = 0;

    // using multimap to sort image cells
    multimap<int,SingleParticlePos_t> N;

    // add the origin
    N.insert(pair<int,SingleParticlePos_t>(0,SingleParticlePos_t(0.0)));

    //summation over the direct lattice
    for(int ix1=-maxx1; ix1<=maxx1; ix1++) {
      for(int ix2=-maxx2; ix2<=maxx2; ix2++) {
  for(int ix3=-maxx3; ix3<=maxx3; ix3++) {
    if(ix1 == 0 && ix2 == 0 && ix3 == 0) continue;
    SingleParticlePos_t xlp(static_cast<T>(ix1),
  			  static_cast<T>(ix2),
  			  static_cast<T>(ix3));

    // using ix1^2+ix2^2+ix3^3 as key to sort the image cells
    N.insert(pair<int,SingleParticlePos_t>(ix1*ix1+ix2*ix2+ix3*ix3, xlp));
  }
      }
    }

    const typename ParticleLayout_t::PtclGrid_t&
      basegrid = *(Lattice.dGrid[glevel]);
    SingleParticlePos_t dx(1.0/basegrid.Delta[0],
  		   1.0/basegrid.Delta[1],
  		   1.0/basegrid.Delta[2]);

    int maxnc = 0;

    M.resize(basegrid.getTotalNum()+1);
    int ntot = basegrid.getTotalNum()*N.size();
    if(ID.size() < ntot) ID.reserve(ntot);
    if(BC.size() < ntot) BC.reserve(ntot);

    M[0] = 0;

    //for each grid, search for connected grid
    for(int ig=0; ig<basegrid.getTotalNum(); ig++) {

      SingleParticlePos_t org(basegrid.Node[ig].Ri[0],
  		      basegrid.Node[ig].Ri[1],
  		      basegrid.Node[ig].Ri[2]);

      typename multimap<int,SingleParticlePos_t>::iterator it = N.begin();
      int ibox = 0;
      int nconnected = 0;

      while(it != N.end()) {

  SingleParticlePos_t displ = (*it).second;

  for(int jg=0; jg<basegrid.getTotalNum(); jg++) {

    if(ibox == 0 && ig == jg) continue; // exclude itself
    SingleParticlePos_t nncell(basegrid.Node[jg].Ri[0]+displ[0],
  			     basegrid.Node[jg].Ri[1]+displ[1],
  			     basegrid.Node[jg].Ri[2]+displ[2]);
    SingleParticlePos_t dcell = nncell-org;
    int icx = int(fabs(dcell[0]*dx[0]));
    int icy = int(fabs(dcell[1]*dx[1]));
    int icz = int(fabs(dcell[2]*dx[2]));
    bool connected = false;
    if(icx <= 1 && icy <= 1 && icz <= 1) {
      connected = true; // next cells are always included
    } else {
      if(Lattice.Dot(dcell,dcell) <= 4*RmaxSq) connected = true;
    }
    if(connected) {
      ID.push_back(jg);
      BC.push_back(Lattice.toCart(displ));
      nconnected++; // number of connected cells for a mesh
    }

  }
  it++; ibox++;
      }
      maxnc = max(nconnected,maxnc);
      M[ig+1] = M[ig] + nconnected;
    }
    return maxnc; // return the maxmimum number of connected cells
  }
  */
  void print(ostream& os)
  {
    for(int ig=0; ig<M.size()-1; ig++)
    {
      cout << ig << " has neighboring cell "  << M[ig+1]-M[ig]<< endl;
      for(int ii=M[ig]; ii<M[ig+1]; ii++)
      {
        cout << ID[ii] << " " << BC[ii] << endl;
      }
    }
  }
};


#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: UniformGridConnection.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
