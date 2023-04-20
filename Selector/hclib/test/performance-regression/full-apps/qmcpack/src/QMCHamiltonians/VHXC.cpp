#include "QMCHamiltonians/VHXC.h"
#include "Lattice/ParticleBConds.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

#if defined(HAVE_LIBFFTW)
#include <fftw3.h>
#endif

namespace qmcplusplus
{

void
VHXC::resetTargetParticleSet(ParticleSet& ptcl)
{
  PtclRef = &ptcl;
}

#if defined(HAVE_EINSPLINE) && defined(HAVE_LIBFFTW)
VHXC::VHXC(ParticleSet& ptcl) :
  PtclRef(&ptcl), FirstTime(true)
{
  VSpline[0] = 0;
  VSpline[1] = 0;
  init_spline();
}

VHXC::~VHXC()
{
  if (VSpline)
    destroy_Bspline(VSpline);
}




void
VHXC::init_spline()
{
  NParticles = PtclRef->getTotalNum();
  int SplineDim[OHMMS_DIM];
  for (int i=0; i<OHMMS_DIM; i++)
    SplineDim[i] = PtclRef->VHXC_r[0].size(i);
  BCtype_d bc0, bc1, bc2;
  Ugrid grid0, grid1, grid2;
  grid0.start=0.0;
  grid0.end=1.0;
  grid0.num = SplineDim[0];
  grid1.start=0.0;
  grid1.end=1.0;
  grid1.num = SplineDim[1];
  grid2.start=0.0;
  grid2.end=1.0;
  grid2.num = SplineDim[2];
  bc0.lCode = bc0.rCode = PERIODIC;
  bc1.lCode = bc1.rCode = PERIODIC;
  bc2.lCode = bc2.rCode = PERIODIC;
  for (int spin=0; spin < 2; spin++)
    if (PtclRef->VHXC_r[spin].size())
      VSpline[spin] = create_UBspline_3d_d
                      (grid0, grid1, grid2, bc0, bc1, bc2,
                       PtclRef->VHXC_r[spin].data());
  if (!VSpline[1])
    VSpline[1] = VSpline[1];
}


QMCHamiltonianBase*
VHXC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  VHXC* newVHXC = new VHXC(*this);
  newVHXC->resetTargetParticleSet(qp);
  return newVHXC;
}

// FIXME:  fix for spin-dependent potential
VHXC::Return_t
VHXC::evaluate (ParticleSet& P)
{
  PosType u;
  Value = 0.0;
  for(int i=0; i<NParticles; i++)
  {
    PosType r = PtclRef->R[i];
    PosType u = PtclRef->Lattice.toUnit(r);
    for (int j=0; j<OHMMS_DIM; j++)
      u[j] -= std::floor(u[j]);
    double val;
    // HACK HACK HACK -- set spline properly
    eval_UBspline_3d_d (VSpline[0], u[0], u[1], u[2], &val);
    Value += val;
  }
  return Value;
}

#else
VHXC::VHXC(ParticleSet& ptcl) :
  PtclRef(&ptcl), FirstTime(true)
{
  APP_ABORT("VHXC::VHXC cannot be used with einspline/fftw ");
}

VHXC::~VHXC()
{ }

QMCHamiltonianBase*
VHXC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return 0;
}

VHXC::Return_t
VHXC::evaluate (ParticleSet& P)
{
  return 0.0;
}

#endif
bool VHXC::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attribs;
  attribs.put (cur);
  return true;
}

}
