//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ASYMMETRICDISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_ASYMMETRICDISTANCETABLEDATAIMPL_H

namespace qmcplusplus
{

/**@ingroup nnlist
 * @brief A derived classe from DistacneTableData, specialized for dense-asymmetric case
 *
 * AsymmetricDTD stands for Asymmetric Distance Table Data with
 * distinct source and target sets.
 * The template parameter BC provides BC::apply member function
 * to evaluate Cartesian distances.
 */
template<typename T, unsigned D, int SC>
struct AsymmetricDTD
    : public DTD_BConds<T,D,SC>, public DistanceTableData
{
  const ParticleSet& Target;
  AsymmetricDTD(const ParticleSet& source,
                const ParticleSet& target)
    : DTD_BConds<T,D,SC>(source.Lattice), DistanceTableData(source,target)
    , Target(target)
  {
    create(1);
  }

  void create(int walkers)
  {
    int nw = (walkers>0)? walkers:1;
    reset(Origin->getTotalNum(),Target.getTotalNum(),nw);
  }

  /*!\fn  void reset(int n1, int n2, int nactive){
   *\param n1 the number of sources
   *\param n2 the number of targets
   *\param nactive the number of copies of the targets
   *\brief Resize the internal data and assign indices
   */
  inline void reset(int n1, int n2, int nactive)
  {
    if( n1!=N[SourceIndex] || n2 != N[VisitorIndex]
        || nactive != N[WalkerIndex])
    {
      N[SourceIndex] = n1;
      N[VisitorIndex] = n2;
      int m = n1*n2;
      if(m)
      {
        M.resize(n1+1);
        J.resize(m);
        PairID.resize(m);
        resize(m,nactive);
        M[0] = 0;
        int ij = 0;
        for(int i=0; i<n1; i++)
        {
          for(int j=0; j<n2; j++, ij++)
          {
            J[ij] = j;
            PairID[ij] = Origin->GroupID[i];
          }
          M[i+1] = M[i]+n2;
        }
        npairs_m = n1*n2;
      }
    }
  }
       
  inline void setTranspose()
  {
    const int ns = N[SourceIndex];
    const int nv = N[VisitorIndex];
    for(int i=0; i<nv; i++)
      for(int j=0; j<ns; j++)
        trans_r(i,j)=r_m[j*nv+i];
  
    if(NeedDisplacement)
      for(int i=0; i<nv; i++)
        for(int j=0; j<ns; j++)
          trans_dr(i,j)=dr_m[j*nv+i];
  }

  inline virtual void nearest_neighbor(vector<ripair>& ri,bool transposed=false) const
  {
    if(transposed)
    {
      for(int n=0; n<ri.size(); ++n)
        ri[n].first = 1e99;
      const int m = N[SourceIndex];
      const int nv = N[VisitorIndex];
      int shift = 0;
      for(int i=0; i<m; ++i,shift+=nv)
        for(int n=0; n<ri.size(); ++n)
        {
          ripair& rin = ri[n];
          RealType rp = r_m[shift+n];
          if(rp<rin.first)
          {
            rin.first  = rp;
            rin.second = i;
          }
        }
    }
    else
    {
      const int m = N[VisitorIndex];
      for(int n=0; n<ri.size(); ++n)
      {
        const int shift = M[n];
        ripair& rin = ri[n];
        rin.first = 1e99;
        for(int i=0; i<m; ++i)
        {
          RealType rp = r_m[shift+i];
          if(rp<rin.first)
          {
            rin.first  = rp;
            rin.second = i;
          }
        }
      }
    }
  }

  inline virtual void nearest_neighbors(int n,int neighbors,vector<ripair>& ri,bool transposed=false)
  {
    if(transposed)
    {
      const int m = N[SourceIndex];
      const int nv = N[VisitorIndex];
      int shift = 0;
      for(int i=0; i<m; ++i,shift+=nv)
      {
        ri[i].first  = r_m[shift+n];
        ri[i].second = i;
      }
    }
    else
    {
      const int m = N[VisitorIndex];
      const int shift = M[n];
      for(int i=0; i<m; ++i)
      {
        ri[i].first  = r_m[shift+i];
        ri[i].second = i;
      }
    }
    partial_sort(ri.begin(),ri.begin()+neighbors,ri.end());
  }

  virtual void nearest_neighbors_by_spec(int n,int neighbors,int spec_start,vector<ripair>& ri,bool transposed=false)
  {
    if(transposed)
    {
      const int nv = N[VisitorIndex];
      int shift = spec_start*nv;
      for(int i=0; i<ri.size(); ++i,shift+=nv)
      {
        ri[i].first  = r_m[shift+n];
        ri[i].second = i;
      }
    }
    else
    {
      const int shift = M[n]+spec_start;
      for(int i=0; i<ri.size(); ++i)
      {
        ri[i].first  = r_m[shift+i];
        ri[i].second = i;
      }
    }
    partial_sort(ri.begin(),ri.begin()+neighbors,ri.end());
  }

  ///evaluate the Distance Table using a set of Particle Positions
  //inline void evaluate(const WalkerSetRef& W) {
  //  int copies = W.walkers();
  //  int visitors = W.particles();
  //  int ns = Origin.getTotalNum();

  //  reset(ns,visitors,copies);
  //  for(int iw=0; iw<copies; iw++) {
  //    int nn=0;
  //    for(int i=0; i<ns; i++) {
  //      PosType r0 = Origin.R(i);
  //      for(int j=0; j<visitors; j++,nn++) {
  //        PosType drij = W.R(iw,j)-r0;
  //        RealType sep = sqrt(BC::apply(Origin.Lattice,drij));
  //#ifdef USE_FASTWALKER
  //        r2_m(nn,iw) = sep;
  //        rinv2_m(nn,iw) = 1.0/sep;
  //        dr2_m(nn,iw) = drij;
  //#else
  //        r2_m(iw,nn) = sep;
  //        rinv2_m(iw,nn) = 1.0/sep;
  //        dr2_m(iw,nn) = drij;
  //#endif
  //      }
  //    }
  //  }
  //}

  ///not so useful inline but who knows
  inline void evaluate(const ParticleSet& P)
  {
    for(int i=0,ij=0; i<N[SourceIndex]; i++)
      for(int j=0; j<N[VisitorIndex]; j++,ij++)
        dr_m[ij]=P.R[j]-Origin->R[i];
    //BC::apply(Origin.Lattice,dr_m,r_m,rinv_m);
    DTD_BConds<T,D,SC>::apply_bc(dr_m,r_m,rinv_m);
    ////reset(Origin.getTotalNum(),P.getTotalNum(),1);
    //int nn=0;
    //for(int i=0; i<N[SourceIndex]; i++) {
    //  PosType r0(Origin.R[i]);
    //  for(int j=0; j<N[VisitorIndex]; j++,nn++) {
    //    PosType drij(P.R[j]-r0);
    //    RealType sep2(BC::apply(Origin.Lattice,drij));
    //    RealType sep(sqrt(sep2));
    //    r_m[nn]    = sep;
    //    //rr_m[nn]   = sep2;
    //    rinv_m[nn] = 1.0/sep;
    //    dr_m[nn]   = drij;
    //  }
    //}
  }

  ///evaluate the temporary pair relations
  inline void move(const ParticleSet& P, const PosType& rnew, IndexType jat)
  {
    activePtcl=jat;
    for(int iat=0, loc=jat; iat<N[SourceIndex]; iat++,loc+=N[VisitorIndex])
    {
      PosType drij(rnew-Origin->R[iat]);
      //RealType sep2(BC::apply(Origin.Lattice,drij));
      RealType sep(std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij)));
      Temp[iat].r1=sep;
      Temp[iat].rinv1=1.0/sep;
      Temp[iat].dr1=drij;
      //Temp[iat].r0=r_m[loc];
      //Temp[iat].rinv0=rinv_m[loc];
      //Temp[iat].dr0=dr_m[loc];
    }
  }

  ///evaluate the temporary pair relations
  inline void moveby(const ParticleSet& P, const PosType& displ, IndexType jat)
  {
    activePtcl=jat;
    for(int ic=0, loc=jat; ic<N[SourceIndex]; ic++,loc+=N[VisitorIndex])
      temp_dr[ic]=displ+dr_m[loc];
    DTD_BConds<T,D,SC>::apply_bc(temp_dr,temp_r);
    //DTD_BConds<T,D,SC>::get_min_distanceX(P.Lattice,Temp);
  }


  ///evaluate the temporary pair relations
  inline void moveOnSphere(const ParticleSet& P, const PosType& rnew, IndexType jat)
  {
    activePtcl=jat;
    for(int iat=0, loc=jat; iat<N[SourceIndex]; iat++,loc+=N[VisitorIndex])
    {
      PosType drij(rnew-Origin->R[iat]);
      Temp[iat].r1=std::sqrt(DTD_BConds<T,D,SC>::apply_bc(drij));
      Temp[iat].dr1=drij;
    }
  }

  inline void update(IndexType jat)
  {
    for(int iat=0,loc=jat; iat<N[SourceIndex]; iat++, loc += N[VisitorIndex])
    {
      r_m[loc]=Temp[iat].r1;
      rinv_m[loc]=1/Temp[iat].r1;
      //rinv_m[loc]=Temp[iat].rinv1;
      dr_m[loc]=Temp[iat].dr1;
    }
  }
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 6361 $   $Date: 2014-09-24 15:18:11 -0400 (Wed, 24 Sep 2014) $
 * $Id: AsymmetricDistanceTableData.h 6361 2014-09-24 19:18:11Z jtkrogel $
 ***************************************************************************/
