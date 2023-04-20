// -*- C++ -*-
/*! \author Jordan Vincent
 *  \author Curry Taylor
 *  \note  The original Prim was written in F90 by Tim Wilkens.
 */
#ifndef CLEBSCH_GORDAN_H2
#define CLEBSCH_GORDAN_H2

#include <blitz/array.h>

/**class Clebsch_Gordan
 *\brief Calculates the Clebsch-Gordan coefficients
 */
class Clebsch_Gordan
{
public:

  Clebsch_Gordan(const int lmax);
  ///destructor
  ~Clebsch_Gordan();

  ///maximum angular momentum
  int Lmax;
  ///array to store the Clebsch-Gordan coefficients
  blitz::Array<double,5> cg;
  ///returns \f$c_g(l1,l2,l3,m1,m2) \f$
  inline double operator()(int l1, int l2, int l3, int m1, int m2) const
  {
    return cg(l1,l2,l3,m1+Lmax,m2+Lmax);
  }

private:
  /// default constructor not implemented
  Clebsch_Gordan() { }

  void build_coefficients();
};

#endif


/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: Clebsch_Gordan.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

