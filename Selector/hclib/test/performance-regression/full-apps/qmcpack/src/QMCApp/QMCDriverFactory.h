//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#define QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "QMCApp/ParticleSetPool.h"
#include <bitset>

class Communicate;

namespace qmcplusplus
{

//forward declaration
class MCWalkerConfiguration;
class QMCDriver;
class WaveFunctionPool;
class HamiltonianPool;

struct QMCDriverFactory: public MPIObjectBase
{
  /*! enum for QMC Run Type */
  enum QMCRunType
  {
    DUMMY_RUN, /*!< dummy */
    VMC_RUN, /**< VMC type: vmc, vmc-ptcl, vmc-multiple, vmc-ptcl-multiple */
    CSVMC_RUN,
    DMC_RUN, /**< DMC type: dmc, dmc-ptcl*/
    RMC_RUN, /**< RMC type: rmc, rmc-ptcl */
    OPTIMIZE_RUN,/*!< Optimization */
    VMC_OPT_RUN, /*!< Optimization with vmc blocks */
    LINEAR_OPTIMIZE_RUN,
    CS_LINEAR_OPTIMIZE_RUN,
  };

  /*! enum to set the bit to determine the QMC mode */
  enum QMCModeEnum
  {
    UPDATE_MODE,  /**< bit for move: walker or pbyp */
    MULTIPLE_MODE, /**< bit for multple configuration */
    SPACEWARP_MODE, /**< bit for space-warping */
    ALTERNATE_MODE, /**< bit for performing various analysis and weird qmc methods */
    GPU_MODE,     /**< bit to use GPU driver */
    QMC_MODE_MAX=8
  };

  ///current QMC mode determined by curQmcModeBits
  unsigned long curQmcMode;

  ///8-bit (ALTERNATE_MODE,SPACEWARP_MODE, MULTIPLE_MODE, UPDATE_MODE)
  std::bitset<QMC_MODE_MAX> curQmcModeBits;

  ///type of qmcdriver
  QMCRunType curRunType;

  ///name of the current QMCriver
  std::string curMethod;

  /** current MCWalkerConfiguration
   */
  MCWalkerConfiguration *qmcSystem;

  /** current QMCDriver
   */
  QMCDriver *qmcDriver;

  /** ParticleSet Pool
   */
  ParticleSetPool* ptclPool;

  /** TrialWaveFunction Pool
   */
  WaveFunctionPool* psiPool;

  /** QMCHamiltonian Pool
   */
  HamiltonianPool* hamPool;

  /** default constructor **/
  QMCDriverFactory(Communicate* c);

  /** set the active qmcDriver */
  void putCommunicator(xmlNodePtr cur);

  /** set the active qmcDriver */
  bool setQMCDriver(int curSeries, xmlNodePtr cur);

  /** create a new QMCDriver
   */
  void createQMCDriver(xmlNodePtr cur);

  /** virtual destructor **/
  virtual ~QMCDriverFactory();
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: j.k.rofling@gmail.com $
 * $Revision: 5983 $   $Date: 2013-09-17 19:56:36 -0400 (Tue, 17 Sep 2013) $
 * $Id: QMCDriverFactory.h 5983 2013-09-17 23:56:36Z j.k.rofling@gmail.com $
 ***************************************************************************/

