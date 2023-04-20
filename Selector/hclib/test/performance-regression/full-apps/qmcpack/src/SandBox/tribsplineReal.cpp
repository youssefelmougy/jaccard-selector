#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "Optimize/VarList.h"
#include "Numerics/TricubicBsplineSet.h"
#include "Numerics/OneDimGridBase.h"
#include "SandBox/TestFunc.h"
#include "Utilities/Clock.h"
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  double ri = 0.0;
  double rf = 1.0;
  std::vector<int> npts(3);
  npts[0]=101;
  npts[1]=101;
  npts[2]=101;
  double xcut=0.23;
  double ycut=0.67;
  const int nk0=1;
  const int nk1=1;
  const int nk2=1;
  typedef LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,npts[0]);
  gridY.set(ri,rf,npts[1]);
  gridZ.set(ri,rf,npts[2]);
  typedef double value_type;
  //Create an analytic function for assignment
  ComboFunc infunc;
  infunc.push_back(0.5 ,new TestFunc(1,1,1));
  infunc.push_back(0.3 ,new TestFunc(1,1,2));
  infunc.push_back(0.01,new TestFunc(5,3,2));
  infunc.push_back(0.01,new TestFunc(5,7,1));
  //infunc.push_back(0.1,new TestFunc(1,2,1));
  //infunc.push_back(0.01,new TestFunc(2,1,1));
  //infunc.push_back(0.01,new TestFunc(2,2,1));
  //infunc.push_back(0.001,new TestFunc(2,1,2));
  //infunc.push_back(0.001,new TestFunc(2,2,2));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //infunc.push_back(-0.3,new TestFunc(7,2,3));
  //infunc.push_back(0.01,new TestFunc(7,7,7));
  //infunc.push_back(0.001,new TestFunc(5,5,5));
  //Write to an array
  Array<value_type,3> inData(npts[0],npts[1],npts[2]);
  value_type *it=inData.data();
  for(int ix=0; ix<npts[0]; ix++)
  {
    double x(gridX(ix));
    for(int iy=0; iy<npts[1]; iy++)
    {
      double y(gridY(iy));
      for(int iz=0; iz<npts[2]; iz++)
      {
        TinyVector<double,3> pos(x,y,gridZ(iz));
        (*it)=infunc.f(pos);
        ++it;
      }
    }
  }
  TricubicBspline<value_type> aorb;
  aorb.setGrid(ri,rf,ri,rf,ri,rf,npts[0]-1,npts[1]-1,npts[2]-1);
  aorb.Init(inData);
  value_type lap,val,lap0,val0,g;
  TinyVector<value_type,3> grad,grad0;
  double x0=xcut;
  double y0=ycut;
  double z=-0.5, dz=0.01;
  std::ofstream fout("bs.value.dat");
  std::ofstream dfout("bs.grad.dat");
  std::ofstream d2fout("bs.lap.dat");
  fout.setf(std::ios::scientific, std::ios::floatfield);
  dfout.setf(std::ios::scientific, std::ios::floatfield);
  d2fout.setf(std::ios::scientific, std::ios::floatfield);
  while(z<1.5)
  {
    TinyVector<double,3> pos(x0,y0,z);
    val=aorb.evaluate(pos,grad,lap);
    val0=infunc.f(pos);
    grad0=infunc.df(pos);
    lap0=infunc.d2f(pos);
    fout << z
         << std::setw(20) << val0 << std::setw(20) << val  << endl;
    dfout << z
          << std::setw(20) << grad0[0]<< std::setw(20) << grad[0]
          << std::setw(20) << grad0[1]<< std::setw(20) << grad[1]
          << std::setw(20) << grad0[2]<< std::setw(20) << grad[2]
          << std::endl;
    d2fout << z
           << std::setw(20) << lap0 << std::setw(20) <<  lap << std::endl;
    z+=dz;
  }
  return 0;
}
