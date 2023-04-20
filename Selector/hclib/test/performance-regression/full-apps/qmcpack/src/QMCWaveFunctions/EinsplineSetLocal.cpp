//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/EinsplineSetLocal.h"

namespace qmcplusplus
{

void EinsplineSetLocal::evaluate (const ParticleSet& P, int iat,
                                  ValueVector_t& psi)
{
  PosType r (P.R[iat]);
  PosType ru(PrimLattice.toUnit(P.R[iat]));
  ru[0] -= std::floor (ru[0]);
  ru[1] -= std::floor (ru[1]);
  ru[2] -= std::floor (ru[2]);
  for(int j=0; j<OrbitalSetSize; j++)
  {
    complex<double> val;
    Orbitals[j]->evaluate(ru, val);
    double phase = -dot(r, Orbitals[j]->kVec);
    double s,c;
    sincos (phase, &s, &c);
    complex<double> e_mikr (c,s);
    val *= e_mikr;
    convert(val,psi[j]);
  }
}

void EinsplineSetLocal::evaluate (const ParticleSet& P, int iat,
                                  ValueVector_t& psi, GradVector_t& dpsi,
                                  ValueVector_t& d2psi)
{
  PosType r (P.R[iat]);
  PosType ru(PrimLattice.toUnit(P.R[iat]));
  ru[0] -= std::floor (ru[0]);
  ru[1] -= std::floor (ru[1]);
  ru[2] -= std::floor (ru[2]);
  complex<double> val;
  TinyVector<complex<double>,3> gu;
  Tensor<complex<double>,3> hess;
  complex<double> eye (0.0, 1.0);
  for(int j=0; j<OrbitalSetSize; j++)
  {
    complex<double> u;
    TinyVector<complex<double>,3> gradu;
    complex<double> laplu;
    Orbitals[j]->evaluate(ru, val, gu, hess);
    u  = val;
    // Compute gradient in cartesian coordinates
    gradu = dot(PrimLattice.G, gu);
    laplu = trace(hess, GGt);
    PosType k = Orbitals[j]->kVec;
    TinyVector<complex<double>,3> ck;
    ck[0]=k[0];
    ck[1]=k[1];
    ck[2]=k[2];
    double s,c;
    double phase = -dot(P.R[iat], k);
    sincos (phase, &s, &c);
    complex<double> e_mikr (c,s);
    convert(e_mikr*u,psi[j]);
    convert(e_mikr*(-eye * ck * u + gradu),dpsi[j]);
    convert(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu),d2psi[j]);
//#ifdef QMC_COMPLEX
//      psi[j]   = e_mikr * u;
//      dpsi[j]  = e_mikr*(-eye * ck * u + gradu);
//      d2psi[j] = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
//#else
//      psi[j]   = real(e_mikr * u);
//      dpsi[j]  = real(e_mikr*(-eye * ck * u + gradu));
//      d2psi[j] = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
//#endif
  }
}

void
EinsplineSetLocal::evaluate_notranspose (const ParticleSet& P, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads,
    ValueMatrix_t& lapls)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    complex<double> val;
    TinyVector<complex<double>,3> gu;
    Tensor<complex<double>,3> hess;
    complex<double> eye (0.0, 1.0);
    for(int j=0; j<OrbitalSetSize; j++)
    {
      complex<double> u;
      TinyVector<complex<double>,3> gradu;
      complex<double> laplu;
      Orbitals[j]->evaluate(ru, val, gu, hess);
      u  = val;
      gradu = dot(PrimLattice.G, gu);
      laplu = trace(hess, GGt);
      PosType k = Orbitals[j]->kVec;
      TinyVector<complex<double>,3> ck;
      ck[0]=k[0];
      ck[1]=k[1];
      ck[2]=k[2];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      convert(e_mikr * u, vals(i,j));
      convert(e_mikr*(-eye*u*ck + gradu),grads(i,j));
      convert(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu),lapls(i,j));
//#ifdef QMC_COMPLEX
//	vals(i,j)  = e_mikr * u;
//	grads(i,j) = e_mikr*(-eye*u*ck + gradu);
//	lapls(i,j) = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
//#else
//	vals(i,j)  = real(e_mikr * u);
//	grads(i,j) = real(e_mikr*(-eye*u*ck + gradu));
//	lapls(i,j) = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
//#endif
    }
  }
}

void
EinsplineSetLocal::evaluate_notranspose (const ParticleSet& P, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads,
    HessMatrix_t& grad_grad_psi)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(r));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    complex<double> val;
    TinyVector<complex<double>,3> gu;
    Tensor<complex<double>,3> hess;
    complex<double> eye (0.0, 1.0);
    for(int j=0; j<OrbitalSetSize; j++)
    {
      complex<double> u;
      TinyVector<complex<double>,3> gradu;
      Tensor<complex<double>,3> tmphs,hs;
      Orbitals[j]->evaluate(ru, val, gu, hess);
      u  = val;
      gradu = dot(PrimLattice.G, gu);
// FIX FIX FIX: store transpose(PrimLattice.G) to avoid recalculation
      tmphs = dot(transpose(PrimLattice.G),hess);
      hs = dot(tmphs,PrimLattice.G);
      PosType k = Orbitals[j]->kVec;
      TinyVector<complex<double>,3> ck;
      ck[0]=k[0];
      ck[1]=k[1];
      ck[2]=k[2];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      convert(e_mikr * u, vals(i,j));
      convert(e_mikr*(-eye*u*ck + gradu),grads(i,j));
      convert(e_mikr*(hs -u*outerProduct(ck,ck) - eye*outerProduct(ck,gradu) - eye*outerProduct(gradu,ck)),grad_grad_psi(i,j));
    }
  }
}

void EinsplineSetLocal::evaluate_notranspose(const ParticleSet& P, int first, int last,
    ValueMatrix_t& psi, GradMatrix_t& dpsi,
    HessMatrix_t& grad_grad_psi,
    GGGMatrix_t& grad_grad_grad_logdet)
{
  APP_ABORT(" EinsplineSetExtended<StorageType>::evaluate_notranspose not implemented for grad_grad_grad_logdet yet. \n");
}



SPOSetBase*
EinsplineSetLocal::makeClone() const
{
  return new EinsplineSetLocal(*this);
}

void
EinsplineSetLocal::resetParameters(const opt_variables_type& active)
{
}

}
