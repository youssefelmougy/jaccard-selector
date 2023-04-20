//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H
#define QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H

#include "OhmmsPETE/OhmmsArray.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

namespace qmcplusplus
{

/** base class for Single-particle orbital sets
 *
 * SPOSetBase stands for S(ingle)P(article)O(rbital)SetBase which contains
 * a number of single-particle orbitals with capabilities of evaluating \f$ \psi_j({\bf r}_i)\f$
 */
class SPOSetBase: public QMCTraits
{
public:
  typedef OrbitalSetTraits<ValueType>::IndexVector_t IndexVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
  typedef OrbitalSetTraits<ValueType>::ValueMatrix_t ValueMatrix_t;
  typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
  typedef OrbitalSetTraits<ValueType>::GradMatrix_t  GradMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessVector_t  HessVector_t;
  typedef OrbitalSetTraits<ValueType>::HessMatrix_t  HessMatrix_t;
  typedef OrbitalSetTraits<ValueType>::HessType      HessType;
  typedef Array<HessType,OHMMS_DIM>                  HessArray_t;
  typedef OrbitalSetTraits<ValueType>::GradHessType  GGGType;
  typedef OrbitalSetTraits<ValueType>::GradHessVector_t GGGVector_t;
  typedef OrbitalSetTraits<ValueType>::GradHessMatrix_t GGGMatrix_t;
  typedef ParticleSet::Walker_t                      Walker_t;
  typedef std::map<string,SPOSetBase*> SPOPool_t;

  ///index in the builder list of sposets
  int builder_index;
  ///true if C is an identity matrix
  bool Identity;
  ///true if C is an identity matrix
  bool Optimizable;
  ///flag to calculate ionic derivatives
  bool ionDerivs;
  ///total number of orbitals
  IndexType TotalOrbitalSize;
  ///number of Single-particle orbitals
  IndexType OrbitalSetSize;
  ///number of Single-particle orbitals
  IndexType BasisSetSize;
  ///index of the particle
  IndexType ActivePtcl;
  ///matrix to store temporary value before transpose
  ValueMatrix_t t_logpsi;
  ///matrix containing the coefficients
  ValueMatrix_t C;
  ///occupation number
  Vector<RealType> Occ;
  /// Optimizable variables
  opt_variables_type myVars;
  ///name of the basis set
  string className;
  /** name of the object
   *
   * Several user classes can own SPOSetBase and use objectName as counter
   */
  string objectName;

  /** constructor
   */
  SPOSetBase()
    :Identity(false),TotalOrbitalSize(0),OrbitalSetSize(0),BasisSetSize(0),
    ActivePtcl(-1),Optimizable(false),ionDerivs(false),builder_index(-1)
  {
    className="invalid";
  }

  /** destructor
   */
  virtual ~SPOSetBase() {}

  /** return the size of the orbital set
   */
  inline int size() const
  {
    return OrbitalSetSize;
  }

  /** print basic SPOSet information
   */
  void basic_report(const string& pad="");
  
  /** print SPOSet information
   */
  virtual void report(const string& pad="")
  {
    basic_report(pad);
  }


  /** return the size of the orbital set
   */
  inline int getOrbitalSetSize() const
  {
    return OrbitalSetSize;
  }

  inline int getBasisSetSize() const
  {
    return BasisSetSize;
  }


  bool setIdentity(bool useIdentity)
  {
    Identity=useIdentity;
    C.resize(OrbitalSetSize,BasisSetSize);
    for(int i=0; i<OrbitalSetSize; i++)
      C(i,i)=1.0;
    return true;
  }

  void checkObject();

  ///get C and Occ
  bool put(xmlNodePtr cur);

  virtual bool put(xmlNodePtr cur, SPOPool_t &spo_pool)
  {
    return put(cur);
  }

  ///reset
  virtual void resetParameters(const opt_variables_type& optVariables)=0;

  virtual void checkInVariables(opt_variables_type& active) {}
  virtual void checkOutVariables(const opt_variables_type& active) {}

  // Evaluate the derivative of the optimized orbitals with
  // respect to the parameters
  virtual void evaluateDerivatives
  (ParticleSet& P, int iat, const opt_variables_type& active,
   ValueMatrix_t& d_phi, ValueMatrix_t& d_lapl_phi) {}


  ///reset the target particleset
  virtual void resetTargetParticleSet(ParticleSet& P)=0;
  /** set the OrbitalSetSize
   * @param norbs number of single-particle orbitals
   */
  virtual void setOrbitalSetSize(int norbs)=0;

  virtual void
  evaluate (const ParticleSet& P, PosType &r, ValueVector_t &psi)
  {
    app_error() << "Need specialization for SPOSetBase::evaluate "
                << "(const ParticleSet& P, PosType &r).\n";
    abort();
  }

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)=0;

  /** evaluate values for the virtual moves, e.g., sphere move for nonlocalPP
   * @param VP virtual particle set
   * @param psiM single-particle orbitals psiM(i,j) for the i-th particle and the j-th orbital
   */
  virtual void
  evaluateValues(const ParticleSet& VP, ValueMatrix_t& psiM);

  /** evaluate the values, gradients and laplacians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)=0;

  /** evaluate the values, gradients and hessians of this single-particle orbital set
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)=0;

  /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last)particles
   * @param P current ParticleSet
   * @param first starting index of the particles
   * @param last ending index of the particles
   * @param logdet determinant matrix to be inverted
   * @param dlogdet gradients
   * @param d2logdet laplacians
   *
   * Call evaluate_notranspose to build logdet
   */
  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  virtual void
  evaluate(const ParticleSet& P, int first, int last
           , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  virtual void
  evaluateThirdDeriv(const ParticleSet& P, int first, int last
                     , GGGMatrix_t& grad_grad_grad_logdet);

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)=0;

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet);

  virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
                                    , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src, GradMatrix_t &gradphi);

  virtual void evaluateGradSource (const ParticleSet &P, int first, int last
                                   , const ParticleSet &source, int iat_src
                                   , GradMatrix_t &grad_phi, HessMatrix_t &grad_grad_phi, GradMatrix_t &grad_lapl_phi);

  virtual void evaluateBasis (const ParticleSet &P, int first, int last
                              , ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl);

  virtual void evaluateForDeriv (const ParticleSet &P, int first, int last
                                 , ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl);

  virtual inline void setpm(int x) {};

  virtual void copyParamsFromMatrix (const opt_variables_type& active
                                     , const ValueMatrix_t &mat, vector<RealType> &destVec);

  virtual PosType get_k(int orb)
  {
    return PosType();
  }

  /** make a clone of itself
   */
  virtual SPOSetBase* makeClone() const;

  virtual bool transformSPOSet()
  {
    return true;
  }


#ifdef QMC_CUDA

  /** evaluate the values of this single-particle orbital set
   * @param P current ParticleSet
   * @param r is the position of the particle
   * @param psi values of the SPO
   */
  virtual void
  evaluate (const ParticleSet& P, const PosType& r, vector<RealType> &psi);

  virtual void initGPU() {  }

  //////////////////////////////////////////
  // Walker-parallel vectorized functions //
  //////////////////////////////////////////
  virtual void
  reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool) { }

  virtual void
  evaluate (vector<Walker_t*> &walkers, int iat, gpu::device_vector<CudaValueType*> &phi);

  virtual void evaluate (vector<Walker_t*> &walkers, vector<PosType> &new_pos
                         , gpu::device_vector<CudaValueType*> &phi);

  virtual void
  evaluate (vector<Walker_t*> &walkers,
            vector<PosType> &new_pos,
            gpu::device_vector<CudaValueType*> &phi,
            gpu::device_vector<CudaValueType*> &grad_lapl_list,
            int row_stride);

  virtual void
  evaluate (vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi);
  virtual void
  evaluate (vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi);
#endif

protected:
  bool putOccupation(xmlNodePtr occ_ptr);
  bool putFromXML(xmlNodePtr coeff_ptr);
  bool putFromH5(const char* fname, xmlNodePtr coeff_ptr);
};

#if defined(ENABLE_SMARTPOINTER)
typedef boost::shared_ptr<SPOSetBase> SPOSetBasePtr;
#else
typedef SPOSetBase*                   SPOSetBasePtr;
#endif

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 6255 $   $Date: 2014-02-25 11:00:26 -0500 (Tue, 25 Feb 2014) $
 * $Id: SPOSetBase.h 6255 2014-02-25 16:00:26Z jnkim $
 ***************************************************************************/

