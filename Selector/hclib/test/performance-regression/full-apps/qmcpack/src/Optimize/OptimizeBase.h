//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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

#ifndef QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H
#define QMCPLUSPLUS_OPTIMIZATIONFUNCION_BASE_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Optimize/LeastSquaredFit.h"

/** Base class for any cost function
 *
 * Based on K. Esler's MinimizeFunction.
 */
template<class T=double>
class CostFunctionBase
{

public:

  typedef T Return_t;

  /** boolean to indicate if the cost function is valid.
   *
   * Can be used by optimizers to stop optimization.
   */
  bool IsValid;

  CostFunctionBase():IsValid(true) { }

  virtual ~CostFunctionBase() {}

  virtual int NumParams() = 0;

  virtual Return_t& Params(int i) = 0;

  virtual Return_t Params(int i) const = 0;

  virtual Return_t Cost(bool needGrad=true) = 0;

  virtual void GradCost(vector<Return_t>& PGradient, const vector<Return_t>& PM, Return_t FiniteDiff=0) = 0;

  virtual void Report() = 0;

  /** find dl that minimize an object function in cg direction
   * @param x0 current parameters
   * @param cg direction
   * @param dl displacement
   * @param val optimal value
   * @param lambda_max maximum displacement
   * @return true, if lineoptimization is done by the derived class
   *
   * Default implementation returns false to perform a simple method
   */
  virtual bool lineoptimization(const std::vector<T>& x0, const std::vector<T>& gr,
                                Return_t val0,
                                Return_t& dl, Return_t& vopt, Return_t& lambda_max)
  {
    return false;
  }

  /** evaluate gradient, \f$ gr(\alpha_i)= -\frac{\partial F}{\partial \alhpa_i}\f$
   * @param gradients
   * @return true, if a cost function evaluates the gradients
   *
   * Default implementation returns false to perform a finite-difference method
   * for gradients.
   */
  virtual bool evaluateGradients(std::vector<T>& gr)
  {
    return false;
  }
};

/** base class for optimization(minimization) classes
 *
 * Template parameter T is the numerical type
 */
template<class T=double>
struct MinimizerBase
{

  /** stream to write intermediate message
   */
  ostream* msg_stream;

  /** typedef of the object function to be optimized
   */
  typedef CostFunctionBase<T> ObjectFuncType;

  /** default constructor */
  MinimizerBase():msg_stream(0) {}

  /** virtual destructor */
  virtual ~MinimizerBase() {}

  /** set msg_stream
   * @param os_ptr pointer to ostream
   */
  void setOstream(ostream* os_ptr)
  {
    msg_stream = os_ptr;
  }

  /** optimize an object function
   * @param atarget ObjectFuncType to be optimized
   * @return true if converged
   */
  virtual bool optimize(ObjectFuncType* atarget) = 0;

  virtual bool put(xmlNodePtr cur)=0;
};

#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: OptimizeBase.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
