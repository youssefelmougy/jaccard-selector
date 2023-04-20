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
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/Communicate.h"
#include <io/hdf_archive.h>
#include <limits>

namespace qmcplusplus
{

template<typename T>
inline void transpose(const T* restrict in, T* restrict out, int m)
{
  for(int i=0,ii=0; i<m; ++i)
    for(int j=0,jj=i; j<m; ++j,jj+=m)
      out[ii++]=in[jj];
}

void SPOSetBase::evaluate(const ParticleSet& P, int first, int last,
                          ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  evaluate_notranspose(P,first,last,t_logpsi,dlogdet,d2logdet);
  transpose(t_logpsi.data(),logdet.data(),OrbitalSetSize);
  //evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet);
  //MatrixOperators::transpose(logdet);
}

void SPOSetBase::evaluate(const ParticleSet& P, int first, int last,
                          ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  evaluate_notranspose(P,first,last,t_logpsi,dlogdet,grad_grad_logdet);
  transpose(t_logpsi.data(),logdet.data(),OrbitalSetSize);
  //evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet);
  //MatrixOperators::transpose(logdet);
}

void SPOSetBase::evaluate(const ParticleSet& P, int first, int last,
                          ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
{
  logdet=0;
  evaluate_notranspose(P,first,last,t_logpsi,dlogdet,grad_grad_logdet,grad_grad_grad_logdet);
  transpose(t_logpsi.data(),logdet.data(),OrbitalSetSize);
  //evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet);
  //MatrixOperators::transpose(logdet);
}

void SPOSetBase::evaluateValues(const ParticleSet& P, ValueMatrix_t& psiM)
{
  APP_ABORT("SPOSetBase::evaluate(P,psiM) not implemented.");
}

void SPOSetBase::evaluateThirdDeriv(const ParticleSet& P, int first, int last,
                                    GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateThirdDeriv(). \n");
}

void SPOSetBase::evaluate_notranspose(const ParticleSet& P, int first, int last
                                      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluate_notranspose() for grad_grad_logdet. \n");
}

void SPOSetBase::evaluate_notranspose(const ParticleSet& P, int first, int last,
                                      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluate_notranspose() for grad_grad_grad_logdet. \n");
}


SPOSetBase* SPOSetBase::makeClone() const
{
  APP_ABORT("Missing  SPOSetBase::makeClone for "+className);
  return 0;
}


/** Parse the xml file for information on the Dirac determinants.
 *@param cur the current xmlNode
 */
bool SPOSetBase::put(xmlNodePtr cur)
{
  //initialize the number of orbital by the basis set size
  int norb= BasisSetSize;
  string debugc("no");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(norb,"size");
  aAttrib.add(debugc,"debug");
  aAttrib.put(cur);
  setOrbitalSetSize(norb);
  TotalOrbitalSize=norb;
  //allocate temporary t_logpsi
  t_logpsi.resize(TotalOrbitalSize,OrbitalSetSize);
  const xmlChar* h=xmlGetProp(cur, (const xmlChar*)"href");
  xmlNodePtr occ_ptr=NULL;
  xmlNodePtr coeff_ptr=NULL;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      occ_ptr=cur;
    }
    else if(cname.find("coeff") < cname.size() || cname == "parameter" || cname == "Var")
    {
      coeff_ptr=cur;
    }
    cur=cur->next;
  }
  if(coeff_ptr == NULL)
  {
    app_log() << "   Using Identity for the LCOrbitalSet " << endl;
    return setIdentity(true);
  }
  bool success=putOccupation(occ_ptr);
  if(h == NULL)
    success = putFromXML(coeff_ptr);
  else
    success = putFromH5((const char*)h, coeff_ptr);
  bool success2 = transformSPOSet();
  if(debugc=="yes")
  {
    app_log() << "   Single-particle orbital coefficients dims=" << C.rows() << " x " << C.cols() << endl;
    app_log() << C << endl;
  }
  return success && success2;
}

void SPOSetBase::checkObject()
{
  if(!(OrbitalSetSize == C.rows() && BasisSetSize == C.cols()))
  {
    app_error() << "   SPOSetBase::checkObject Linear coeffient for SPOSet is not consistent with the input." << endl;
    OHMMS::Controller->abort();
  }
}

bool SPOSetBase::putOccupation(xmlNodePtr occ_ptr)
{
  //die??
  if(BasisSetSize ==0)
  {
    APP_ABORT("SPOSetBase::putOccupation detected ZERO BasisSetSize");
    return false;
  }
  Occ.resize(max(BasisSetSize,OrbitalSetSize));
  Occ=0.0;
  for(int i=0; i<OrbitalSetSize; i++)
    Occ[i]=1.0;
  vector<int> occ_in;
  string occ_mode("table");
  if(occ_ptr == NULL)
  {
    occ_mode="ground";
  }
  else
  {
    const xmlChar* o=xmlGetProp(occ_ptr,(const xmlChar*)"mode");
    if(o)
      occ_mode = (const char*)o;
  }
  //Do nothing if mode == ground
  if(occ_mode == "excited")
  {
    putContent(occ_in,occ_ptr);
    for(int k=0; k<occ_in.size(); k++)
    {
      if(occ_in[k]<0) //remove this, -1 is to adjust the base
        Occ[-occ_in[k]-1]=0.0;
      else
        Occ[occ_in[k]-1]=1.0;
    }
  }
  else
    if(occ_mode == "table")
    {
      putContent(Occ,occ_ptr);
    }
  return true;
}

bool SPOSetBase::putFromXML(xmlNodePtr coeff_ptr)
{
  Identity=true;
  int norbs=0;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norbs,"size");
  aAttrib.add(norbs,"orbitals");
  aAttrib.put(coeff_ptr);
  if(norbs < OrbitalSetSize)
  {
    return false;
    APP_ABORT("SPOSetBase::putFromXML missing or incorrect size");
  }
  if(norbs)
  {
    Identity=false;
    vector<ValueType> Ctemp;
    Ctemp.resize(norbs*BasisSetSize);
    setIdentity(Identity);
    putContent(Ctemp,coeff_ptr);
    int n=0,i=0;
    vector<ValueType>::iterator cit(Ctemp.begin());
    while(i<OrbitalSetSize)
    {
      if(Occ[n]>numeric_limits<RealType>::epsilon())
      {
        std::copy(cit,cit+BasisSetSize,C[i]);
        i++;
      }
      n++;
      cit+=BasisSetSize;
    }
  }
  return true;
}

/** read data from a hdf5 file
 * @param norb number of orbitals to be initialized
 * @param fname hdf5 file name
 * @param occ_ptr xmlnode for occupation
 * @param coeff_ptr xmlnode for coefficients
 */
bool SPOSetBase::putFromH5(const char* fname, xmlNodePtr coeff_ptr)
{
#if defined(HAVE_LIBHDF5)
  int norbs=OrbitalSetSize;
  int neigs=BasisSetSize;
  string setname;
  OhmmsAttributeSet aAttrib;
  aAttrib.add(setname,"dataset");
  aAttrib.add(neigs,"size");
  aAttrib.add(neigs,"orbitals");
  aAttrib.put(coeff_ptr);
  setIdentity(false);
  if(setname.empty())
  {
    APP_ABORT("SPOSetBase::putFromH5 missing dataset attribute");
  }
  Matrix<RealType> Ctemp(BasisSetSize,BasisSetSize);
  hdf_archive hin(0);
  hin.open(fname);
  hin.read(Ctemp,setname);
  int n=0,i=0;
  while(i<norbs)
  {
    if(Occ[n]>0.0)
    {
      std::copy(Ctemp[n],Ctemp[n+1],C[i]);
      i++;
    }
    n++;
  }
#else
  APP_ABORT("SPOSetBase::putFromH5 HDF5 is disabled.")
#endif
  return true;
}


void SPOSetBase::basic_report(const string& pad)
{
  app_log()<<pad<<"size = "<<size()<<endl;
  app_log()<<pad<<"state info:"<<endl;
  //states.report(pad+"  ");
  app_log().flush();
}


void SPOSetBase::evaluateBasis (const ParticleSet &P, int first, int last,
                                ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad,
                                ValueMatrix_t &basis_lapl)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateBasis.\n");
}

void SPOSetBase::evaluateForDeriv (const ParticleSet &P, int first, int last,
                                   ValueMatrix_t &basis_val,  GradMatrix_t  &basis_grad,
                                   ValueMatrix_t &basis_lapl)
{
  APP_ABORT("Need specialization of SPOSetBase::evaluateBasis.\n");
}

void SPOSetBase::copyParamsFromMatrix (const opt_variables_type& active,
                                       const ValueMatrix_t &mat, vector<RealType> &destVec)
{
  APP_ABORT("Need specialization of SPOSetBase::copyParamsFromMatrix.");
}

void SPOSetBase::evaluateGradSource (const ParticleSet &P
                                     , int first, int last, const ParticleSet &source
                                     , int iat_src, GradMatrix_t &gradphi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

void SPOSetBase::evaluateGradSource (const ParticleSet &P, int first, int last,
                                     const ParticleSet &source, int iat_src,
                                     GradMatrix_t &grad_phi,
                                     HessMatrix_t &grad_grad_phi,
                                     GradMatrix_t &grad_lapl_phi)
{
  APP_ABORT("SPOSetlBase::evalGradSource is not implemented");
}

#ifdef QMC_CUDA

void SPOSetBase::evaluate(const ParticleSet& P, const PosType& r, vector<RealType> &psi)
{
  APP_ABORT("Not implemented.\n");
}


void SPOSetBase::evaluate (vector<Walker_t*> &walkers, int iat,
                           gpu::device_vector<CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSetBase.\n";
  abort();
}

void SPOSetBase::evaluate (vector<Walker_t*> &walkers, vector<PosType> &new_pos,
                           gpu::device_vector<CudaValueType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate in SPOSetBase.\n";
  abort();
}

void SPOSetBase::evaluate (vector<Walker_t*> &walkers,
                           vector<PosType> &new_pos,
                           gpu::device_vector<CudaValueType*> &phi,
                           gpu::device_vector<CudaValueType*> &grad_lapl_list,
                           int row_stride)
{
  app_error() << "Need specialization of vectorized eval_grad_lapl in SPOSetBase.\n";
  abort();
}

void SPOSetBase::evaluate (vector<PosType> &pos, gpu::device_vector<CudaRealType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSetBase.\n";
  abort();
}

void SPOSetBase::evaluate (vector<PosType> &pos, gpu::device_vector<CudaComplexType*> &phi)
{
  app_error() << "Need specialization of vectorized evaluate "
              << "in SPOSetBase.\n";
  abort();
}
#endif
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 6255 $   $Date: 2014-02-25 11:00:26 -0500 (Tue, 25 Feb 2014) $
 * $Id: SPOSetBase.cpp 6255 2014-02-25 16:00:26Z jnkim $
 ***************************************************************************/

