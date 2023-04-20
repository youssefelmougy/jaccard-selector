//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Kris Delaney and Jeongnim Kim
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
/** @file PWBasis.h
 * @brief Declaration of Plane-wave basis set
 */
#ifndef QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H
#define QMCPLUSPLUS_PLANEWAVEBASIS_BLAS_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Numerics/HDFSTLAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Message/Communicate.h"
#include <complex>

/** If defined, use recursive method to build the basis set for each position
 *
 * performance improvement is questionable: load vs sin/cos
 */
//#define PWBASIS_USE_RECURSIVE

namespace qmcplusplus
{

/** Plane-wave basis set
 *
 * Rewrite of PlaneWaveBasis to utilize blas II or III
 * Support more general input tags
 */
class PWBasis: public QMCTraits
{

public:
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;
  typedef TinyVector<IndexType,3>       GIndex_t;

private:
  ///max of maxg[i]
  int maxmaxg;
  //Need to store the maximum translation in each dimension to use recursive PW generation.
  GIndex_t maxg;
  //The PlaneWave data - keep all of these strictly private to prevent inconsistencies.
  RealType ecut;
  ///twist angle in reduced
  PosType twist;
  ///twist angle in cartesian
  PosType twist_cart; //Twist angle in reduced and Cartesian.

  ///gvecs in reduced coordiates
  vector<GIndex_t> gvecs;
  ///Reduced coordinates with offset gvecs_shifted[][idim]=gvecs[][idim]+maxg[idim]
  vector<GIndex_t> gvecs_shifted;

  vector<RealType> minusModKplusG2;
  vector<PosType>  kplusgvecs_cart; //Cartesian.

  Matrix<ComplexType> C;
  //Real wavefunctions here. Now the basis states are cos(Gr) or sin(Gr), not exp(iGr)
  //We need a way of switching between them for G -> -G, otherwise the
  //determinant will have multiple rows that are equal (to within a constant factor)
  //of others, giving a zero determinant. For this, we build a vector (negative) which
  //stores whether a vector is "+" or "-" (with some criterion, to be defined). We
  //the switch from cos() to sin() based on the value of this input.
  vector<int> negative;
public:
  //enumeration for the value, laplacian, gradients and size
  enum {PW_VALUE, PW_LAP, PW_GRADX, PW_GRADY, PW_GRADZ, PW_MAXINDEX};

  Matrix<ComplexType> Z;

  Vector<ComplexType> Zv;
  /* inputmap is used for a memory efficient way of
   *
   * importing the basis-set and coefficients when the desired energy cutoff may be
   * lower than that represented by all data in the wavefunction input file.
   * The steps taken are:
   *  - Read all basis data.
   *  - Create map. inputmap[i] = j; j is correct PW index, i is input coef index.
   *    For basis elements outside cutoff, inputmap[i] = gvecs.size();
   *  - Coefficients are in same order as PWs in inputfile => simply file into
   *    storage matrix using the map as the input. All excess coefficients are
   *    put into [gvecs.size()] and not used. i.e. coefs need to be allocated 1 higher.
   * Such an approach is not needed for Gamma-point only calculations because the
   * basis is spherically ordered. However, when a twist-angle is used, the "sphere"
   * of allowed planewaves is shifted.
   */
  vector<int> inputmap;

  ///total number of basis functions
  int NumPlaneWaves;

  ///local copy of Lattice
  ParticleLayout_t Lattice;

  ///default constructor
  PWBasis():maxmaxg(0), NumPlaneWaves(0)
  {
  }

  ///constructor
  PWBasis(const PosType& twistangle):
    maxmaxg(0), twist(twistangle), NumPlaneWaves(0)
  {
  }

  ~PWBasis()
  {
  }

  ///basis size
  inline IndexType getBasisSetSize() const
  {
    return NumPlaneWaves;
  }

  ///set the twist angle
  void setTwistAngle(const PosType& tang);

  ///reset
  void reset();

  /** Read basisset from hdf5 file. Apply ecut.
   * @param h5basisgroup h5 node where basis is located
   * @param ecutoff cutoff energy
   * @param lat CrystalLattice
   * @param resizeContainer if true, resize internal storage.
   * @return the number of plane waves
   */
  int
  readbasis(hid_t h5basisgroup,RealType ecutoff,
            ParticleLayout_t &lat,
            const string& pwname="planewaves",
            const string& pwmultname="multipliers",
            bool resizeContainer=true);

  /** Remove basis elements if kinetic energy > ecut.
   *
   * Keep and indexmap so we know how to match coefficients on read.
   */
  void trimforecut();

#if defined(PWBASIS_USE_RECURSIVE)
  /** Fill the recursion coefficients matrix.
   *
   * @todo Generalize to non-orthorohmbic cells
   */
  inline void BuildRecursionCoefs(const PosType& pos)
  {
    PosType tau_red(Lattice.toUnit(pos));
//      RealType phi=TWOPI*tau_red[0];
//      RealType nphi=maxg0*phi;
//      ComplexType ct0(std::cos(phi),std::sin(phi));
//      ComplexType t(std::cos(nphi),-std::sin(nphi));
//      C0[0]=t;
//      for(int n=1; n<=2*maxg0; n++) C0[n] = (t *= ct0);
//
//      phi=TWOPI*tau_red[1];
//      nphi=maxg1*phi;
//      ct0=ComplexType(std::cos(phi),std::sin(phi));
//      t=ComplexType(std::cos(nphi),-std::sin(nphi));
//      C1[0]=t;
//      for(int n=1; n<=2*maxg1; n++) C1[n] = (t *= ct0);
//
//      phi=TWOPI*tau_red[2];
//      nphi=maxg2*phi;
//      ct0=ComplexType(std::cos(phi),std::sin(phi));
//      t=ComplexType(std::cos(nphi),-std::sin(nphi));
//      C2[0]=t;
//      for(int n=1; n<=2*maxg2; n++) C2[n] = (t *= ct0);
#pragma ivdep
    for(int idim=0; idim<3; idim++)
    {
      int ng=maxg[idim];
      RealType phi=TWOPI*tau_red[idim];
      RealType nphi=ng*phi;
      ComplexType Ctemp(std::cos(phi),std::sin(phi));
      ComplexType t(std::cos(nphi),-std::sin(nphi));
      ComplexType* restrict cp_ptr=C[idim];
      *cp_ptr++ = t;
      for(int n=1; n<=2*ng; n++)
      {
        *cp_ptr++ = (t *= Ctemp);
      }
    }
    //Base version
//#pragma ivdep
//      for(int idim=0; idim<3; idim++){
//        RealType phi=TWOPI*tau_red[idim];
//        ComplexType Ctemp(std::cos(phi),std::sin(phi));
//        register int ng=maxg[idim];
//        ComplexType* restrict cp_ptr=C[idim]+ng;
//        ComplexType* restrict cn_ptr=C[idim]+ng-1;
//        *cp_ptr=1.0;
//        for(int n=1; n<=ng; n++,cn_ptr--){
//          ComplexType t(Ctemp*(*cp_ptr++));
//          *cp_ptr = t;
//          *cn_ptr = conj(t);
//        }
//      }
    //Not valid for general supercell
    //      // Cartesian of twist for 1,1,1 (reduced coordinates)
    //      PosType G111(1.0,1.0,1.0);
    //      G111 = Lattice.k_cart(G111);
    //
    //      //Precompute a small number of complex factors (PWs along b1,b2,b3 lines)
    //      //using a fast recursion algorithm
    //#pragma ivdep
    //      for(int idim=0; idim<3; idim++){
    //        //start the recursion with the 111 vector.
    //        RealType phi = pos[idim] * G111[idim];
    //        register ComplexType Ctemp(std::cos(phi), std::sin(phi));
    //        register int ng=maxg[idim];
    //        ComplexType* restrict cp_ptr=C[idim]+ng;
    //        ComplexType* restrict cn_ptr=C[idim]+ng-1;
    //        *cp_ptr=1.0;
    //        for(int n=1; n<=ng; n++,cn_ptr--){
    //          ComplexType t(Ctemp*(*cp_ptr++));
    //          *cp_ptr = t;
    //          *cn_ptr = conj(t);
    //        }
    //      }
  }

  inline void
  evaluate(const PosType& pos)
  {
    BuildRecursionCoefs(pos);
    RealType twistdotr = dot(twist_cart,pos);
    ComplexType pw0(std::cos(twistdotr),std::sin(twistdotr));
    //Evaluate the planewaves for particle iat.
    for(int ig=0; ig<NumPlaneWaves; ig++)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations are for (twist+G).r
      ComplexType pw(pw0); //std::cos(twistdotr),std::sin(twistdotr));
      for(int idim=0; idim<3; idim++)
        pw *= C(idim,gvecs_shifted[ig][idim]);
      //pw *= C0[gvecs_shifted[ig][0]];
      //pw *= C1[gvecs_shifted[ig][1]];
      //pw *= C2[gvecs_shifted[ig][2]];
      Zv[ig]=pw;
    }
  }
  /** Evaluate all planewaves and derivatives for the iat-th particle
   *
   * The basis functions are evaluated for particles iat: first <= iat < last
   * Evaluate the plane-waves at current particle coordinates using a fast
   * recursion algorithm. Order of Y,dY and d2Y is kept correct.
   * These can be "dotted" with coefficients later to complete orbital evaluations.
   */
  inline void
  evaluateAll(const ParticleSet& P, int iat)
  {
    BuildRecursionCoefs(P.R[iat]);
    RealType twistdotr = dot(twist_cart,P.R[iat]);
    ComplexType pw0(std::cos(twistdotr),std::sin(twistdotr));
    //Evaluate the planewaves and derivatives.
    ComplexType* restrict zptr=Z.data();
    for(int ig=0; ig<NumPlaneWaves; ig++,zptr+=5)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations
      //are for (twist+G).r
      ComplexType pw(pw0);
      // THE INDEX ORDER OF C DOESN'T LOOK TOO GOOD: this could be fixed
      for(int idim=0; idim<3; idim++)
        pw *= C(idim,gvecs_shifted[ig][idim]);
      //pw *= C0[gvecs_shifted[ig][0]];
      //pw *= C1[gvecs_shifted[ig][1]];
      //pw *= C2[gvecs_shifted[ig][2]];
      zptr[0]= pw;
      zptr[1]= minusModKplusG2[ig]*pw;
      zptr[2]= kplusgvecs_cart[ig][0]*ComplexType(-pw.imag(),pw.real());
      zptr[3]= kplusgvecs_cart[ig][1]*ComplexType(-pw.imag(),pw.real());
      zptr[4]= kplusgvecs_cart[ig][2]*ComplexType(-pw.imag(),pw.real());
    }
  }
#else
  inline void
  evaluate(const PosType& pos)
  {
    //Evaluate the planewaves for particle iat.
    for(int ig=0; ig<NumPlaneWaves; ig++)
    {
      RealType phi(dot(kplusgvecs_cart[ig],pos));
      Zv[ig]=ComplexType(std::cos(phi),std::sin(phi));
    }
  }
  inline void
  evaluateAll(const ParticleSet& P, int iat)
  {
    ComplexType* restrict zptr=Z.data();
    for(int ig=0; ig<NumPlaneWaves; ig++,zptr+=5)
    {
      //PW is initialized as exp(i*twist.r) so that the final basis evaluations
      //are for (twist+G).r
      RealType phi(dot(kplusgvecs_cart[ig],P.R[iat]));
      ComplexType pw(std::cos(phi),std::sin(phi));
      zptr[0]= pw;
      zptr[1]= minusModKplusG2[ig]*pw;
      zptr[2]= kplusgvecs_cart[ig][0]*ComplexType(-pw.imag(),pw.real());
      zptr[3]= kplusgvecs_cart[ig][1]*ComplexType(-pw.imag(),pw.real());
      zptr[4]= kplusgvecs_cart[ig][2]*ComplexType(-pw.imag(),pw.real());
    }
  }
#endif
//    /** Fill the recursion coefficients matrix.
//     *
//     * @todo Generalize to non-orthorohmbic cells
//     */
//    void BuildRecursionCoefsByAdd(const PosType& pos)
//    {
//      // Cartesian of twist for 1,1,1 (reduced coordinates)
//      PosType G111(1.0,1.0,1.0);
//      G111 = Lattice.k_cart(G111);
//      //PosType redP=P.Lattice.toUnit(P.R[iat]);
//      //Precompute a small number of complex factors (PWs along b1,b2,b3 lines)
//      for(int idim=0; idim<3; idim++){
//        //start the recursion with the 111 vector.
//        RealType phi = pos[idim] * G111[idim];
//        int ng(maxg[idim]);
//        RealType* restrict cp_ptr=logC[idim]+ng;
//        RealType* restrict cn_ptr=logC[idim]+ng-1;
//        *cp_ptr=0.0;
//        //add INTEL vectorization
//        for(int n=1; n<=ng; n++,cn_ptr--){
//          RealType t(phi+*cp_ptr++);
//          *cp_ptr = t;
//          *cn_ptr = -t;
//        }
//      }
//    }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: PWBasis.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
