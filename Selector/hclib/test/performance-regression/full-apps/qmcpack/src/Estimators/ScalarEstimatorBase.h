//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#define QMCPLUSPLUS_SCALAR_ESTIMATORBASE_H
#include <Particle/MCWalkerConfiguration.h>
#include <OhmmsData/RecordProperty.h>
#include <OhmmsData/HDFAttribIO.h>
#include <Estimators/accumulators.h>
#include <Estimators/TraceManager.h>

namespace qmcplusplus
{

class observable_helper;

/** Abstract class for an estimator of a scalar operator.
 *
 * ScalarEstimators derived from ScalarEstimatorBase  implement three main functions
 * - reset : reset the internal values so that observables can be accumulated
 * - accumulate : measure and accumulate its value and the square of the value
 * - report : evaluate the block average and variance
 * ScalarEstimatorBase and its derived classes do not perform any I/O function.
 */
struct ScalarEstimatorBase: public QMCTraits
{

  typedef accumulator_set<RealType> accumulator_type;
  typedef MCWalkerConfiguration::Walker_t       Walker_t;
  typedef MCWalkerConfiguration::const_iterator WalkerIterator;
  typedef RecordNamedProperty<RealType>         RecordListType;

  ///first index within an record of the first element handled by an object
  int FirstIndex;
  ///last index within an record of the first element handled by an object
  int LastIndex;
  ///scalars to be measured
  vector<accumulator_type> scalars;
  ///scalars saved
  vector<accumulator_type> scalars_saved;
//     RealType NSTEPS;
  ///whether traces are being collected
  bool tracing;

  inline ScalarEstimatorBase(): FirstIndex(0), LastIndex(0) {}

  virtual ~ScalarEstimatorBase() {}

  ///checkout trace arrays
  virtual void initialize_traces(TraceManager& tm) {};

  ///delete trace arrays
  virtual void finalize_traces() {};

  ///return average of the
  inline RealType average(int i=0) const
  {
    return scalars_saved[i].mean();
  }
  ///return a variance
  inline RealType variance(int i=0) const
  {
    return scalars_saved[i].variance();
  }
  ///retrun mean and variance
  inline pair<RealType,RealType> operator[](int i) const
  {
    return scalars[i].mean_and_variance();
  }

  ///return the size of scalars it manages
  inline int size() const
  {
    return scalars.size();
  }

  ///clear the scalars to collect
  inline void clear()
  {
    for(int i=0; i<scalars.size(); i++)
      scalars[i].clear();
  }

  /** take block average and write to a common container */
  template<typename IT>
  inline void takeBlockAverage(IT first)
  {
    first += FirstIndex;
    for(int i=0; i<scalars.size(); i++)
    {
      *first++ = scalars[i].mean();
      scalars_saved[i]=scalars[i]; //save current block
      scalars[i].clear();
    }
  }

  /** take block average and write to common containers for values and squared values
   * @param first starting iterator of values
   * @param first_sq starting iterator of squared values
   */
  template<typename IT>
  inline void takeBlockAverage(IT first, IT first_sq)
  {
    first += FirstIndex;
    first_sq += FirstIndex;
    for(int i=0; i<scalars.size(); i++)
    {
      *first++ = scalars[i].mean();
      *first_sq++ = scalars[i].mean2();
      scalars_saved[i]=scalars[i]; //save current block
      scalars[i].clear();
    }
  }

  /** a virtual function to accumulate observables or collectables
   * @param W const MCWalkerConfiguration
   * @param first const_iterator for the first walker
   * @param last const_iterator for the last walker
   * @param wgt weight
   *
   * Pass W along with the iterators so that the properties of W can be utilized.
   */
  virtual void accumulate(const MCWalkerConfiguration& W
                          , WalkerIterator first, WalkerIterator last , RealType wgt) = 0;

  /** add the content of the scalar estimator to the record
   * @param record scalar data list
   *
   * Each ScalarEstimatorBase object adds 1 to many accumulator_type
   */
  virtual void add2Record(RecordNamedProperty<RealType>& record) = 0;

  /** add descriptors of observables to utilize hdf5
   * @param h5desc descriptor of a data stored in a h5 group
   * @param gid h5 group to which each statistical data will be stored
   */
  virtual void registerObservables(vector<observable_helper*>& h5dec, hid_t gid)=0;

  ///clone the object
  virtual ScalarEstimatorBase* clone()=0;

  inline void setNumberOfBlocks(int nsamples)
  {
//       NSTEPS=1.0/ static_cast<RealType>(nsamples);
  }
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5976 $   $Date: 2013-09-13 14:39:44 -0400 (Fri, 13 Sep 2013) $
 * $Id: ScalarEstimatorBase.h 5976 2013-09-13 18:39:44Z jtkrogel $
 ***************************************************************************/
