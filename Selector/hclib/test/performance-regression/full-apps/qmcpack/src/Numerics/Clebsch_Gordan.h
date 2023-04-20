// -*- C++ -*-
/*! \author Jordan Vincent
 *  \author Curry Taylor
 *  \note  The original Prim was written in F90 by Tim Wilkens.
 */
#ifndef CLEBSCH_GORDAN_H
#define CLEBSCH_GORDAN_H

#include <vector>

/** Calculate and store the Clebsch-Gordan coefficients */
class Clebsch_Gordan
{
public:

  Clebsch_Gordan(const int lmax);
  ///destructor
  ~Clebsch_Gordan() { }

  ///vector to store the Clebsch-Gordan coefficients
  std::vector<double> CG_coeff;

  /**
   *@brief return \f$ index = l_1 + l_2 L1max + l_3 L1max^2
   + m_1 L1max^2 L2max + m_2 L1max^2 L2max^2 \f$
   *@note Lmax is internally added to m1 and m2 to offset the index
   */
  inline int index(int l1, int l2, int l3, int m1, int m2) const
  {
    return l1+l2*L1max+l3*L1max*L1max
           +(m1+Lmax)*L1max*L1max*L2max+(m2+Lmax)*L1max*L1max*L2max*L2max;
  }

  ///return \f$ \langle l_1 m_1 l_2 m_2 | l_3 (m_1+m_2) \rangle \f$
  inline double cg(int l1, int l2, int l3, int m1, int m2) const
  {
    return CG_coeff[index(l1,l2,l3,m1,m2)];
  }

private:

  ///maximum angular momentum
  int Lmax;
  //maximum for \f$ l_1 \f$ and \f$ l_2 \f$
  int L1max;
  //maximum for \f$ l_3, m_1 \f$ and \f$ m_2 \f$
  int L2max;
  /// default constructor not implemented
  Clebsch_Gordan() { }

};

#endif


/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: Clebsch_Gordan.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/

