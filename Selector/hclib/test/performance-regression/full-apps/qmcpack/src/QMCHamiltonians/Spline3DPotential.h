#ifndef QMCPLUSPLUS_MYSPLINEDDAS_H
#define QMCPLUSPLUS_MYSPLINEDDAS_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/Spline3D/Grid3D.h"
#include "Numerics/Spline3D/Spline3D.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"


namespace qmcplusplus
{


struct Spline3DPotential: public QMCHamiltonianBase
{

  /// pointer to the main grid (initialised in wavefunction)
  Grid3D* full_Grid;

  /// the spline to calculate the potential
  Spline3D* pot_m;

  /// Constructor
  Spline3DPotential(Grid3D* agrid, const string& fname)
  {
    full_Grid = agrid;
    pot_m = new Spline3D(agrid,agrid->npts_m);
    /// Create the spline from the given grid and initialise from the file
    cout << "Reading Potential File and initialising ... ";
    pot_m->read_data(fname.c_str());
    for(int i = 0; i < pot_m->f.size(); i++)
      pot_m->f[i] *= 0.036749033500418936;
    pot_m->d2fdr2();
    cout << "done! " << endl;
  }

  /// Destructor
  ~Spline3DPotential() { }

  /// evaluate the potential
  inline Return_t evaluate(ParticleSet& P)
  {
    ValueType e = 0.0;
    for(int i=0; i<P.getTotalNum(); i++)
    {
      e+=pot_m->evaluate(P.R[i]);
    }
    return e;
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t evaluate(ParticleSet& P, RealType& x)
  {
    return x=evaluate(P);
  }

  void evaluate(WalkerSetRef& P, ValueVectorType& LE)
  {
  }

};

}
#endif
