//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
// -*- C++ -*-
#ifndef QMCPLUSPLUS_MULTIPLELOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_MULTIPLELOCALENERGYESTIMATOR_H

#include "Estimators/ScalarEstimatorBase.h"
#include "QMCDrivers/SpaceWarp.h"

namespace qmcplusplus
{

class QMCHamiltonian;
class TrialWaveFunction;

struct MultipleEnergyEstimator: public ScalarEstimatorBase
{

  enum {ENERGY_INDEX, ENERGY_SQ_INDEX, WEIGHT_INDEX, PE_INDEX, KE_INDEX, LE_INDEX};

  ///index to keep track how many times accumulate is called
  int CurrentWalker;

  ///total number of walkers
  int NumWalkers;

  ///starting column index of this estimator
  int FirstColumnIndex;

  ///index of the first Hamiltonian element
  int FirstHamiltonian;

  /** the number of duplicated wavefunctions/hamiltonians
   *
   * NumCopies is the row size of Walker_t::Properties
   */
  int NumCopies;

  /** the number of components of a Hamiltonian
   *
   * It also works as the column index for the local energy. For instacne,
   * UmbrellaSamples[walker_index]->operator()(copy_index,NumOperators)
   * is the local energy of the copy_index wavefunction for the walker_index walker
   */
  int NumOperators;

  ///the index of weight
  //int WeightIndex;

  /** the column size of EnergyContaner_t for UmbrellaSamples
   *
   * NumCols = NumOperators+2
   * 2 is for the local energy and weight
  int NumCols;
   */

  /** local energy data for each walker.
   */
  //Matrix<RealType> UmbrellaEnergy, UmbrellaWeight, RatioIJ;
  Matrix<RealType> RatioIJ;

  /** accumulated total local energies, weights
   *
   * The dimension of esum is NumCopies by LE_INDEX
   */
  Matrix<RealType> esum;

  /** accumulated local operators, e.g., Kinetic energy
   *
   * The dimension of elocal is NumCopies by NumOperators.
   */
  Matrix<RealType> elocal;

  /** the names of esum
   *
   * The dimension of esum_name is LE_INDEX by NumCopies. The columns
   * of the same field are bunched together.
   */
  Matrix<string> esum_name;

  /** the names of elocal
   *
   * The dimension of elocal_name is NumOperators by NumCopies. The columns
   * of the same field are bunched together.
   */
  Matrix<string> elocal_name;

  /** the names of energy differences
   */
  vector<string> ediff_name;

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  MultipleEnergyEstimator(QMCHamiltonian& h, int hcopy=1);
  MultipleEnergyEstimator(const MultipleEnergyEstimator& mest);

  inline RealType getUmbrellaWeight(int ipsi)
  {
    return esum(ipsi,WEIGHT_INDEX);
  }

  void accumulate(const Walker_t& awalker, RealType wgt);

  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(; first != last; ++first)
      accumulate(**first,wgt);
  }
  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void add2Record(RecordNamedProperty<RealType>& record);
  void registerObservables(vector<observable_helper*>& h5dec, hid_t gid);
  ScalarEstimatorBase* clone();

  /////reset all the cumulative sums to zero
  //void reset();

  ///** calculate the averages and reset to zero
  // *\param record a container class for storing scalar records (name,value)
  // *\param wgtinv the inverse weight
  // */
  //void report(RecordNamedProperty<RealType>& record, RealType wgtinv);
  //void report(RecordNamedProperty<RealType>& record, RealType wgtinv, BufferType& msg);

  /** initialize the multi-configuration data
   *
   * @param W MCWalkerConfiguration
   * @param h Collection of QMCHamiltonians*
   * @param psi Collection of TrialWaveFunctions*
   * @param tau time step for the initial weighted drift
   * @param require_register if true, use buffer for particle-by-particle
   */
  void initialize(MCWalkerConfiguration& W,
                  vector<QMCHamiltonian*>& h, vector<TrialWaveFunction*>& psi,
                  RealType tau, vector<RealType>& Norm, bool require_register=false);

  /** update the energy and weight for umbrella sampling
   * @param iw walker index
   * @param ipsi H/Psi index
   * @param e local energy of the iw-th walker for H[ipsi]/Psi[ipsi]
   * @param r umbrella weight
  inline void updateSample(int iw, int ipsi, RealType e, RealType invr) {
    //UmbrellaEnergy(iw,ipsi)=e;
    //UmbrellaWeight(iw,ipsi)=invr;
  }
   */

  void initialize(MCWalkerConfiguration& W, vector<ParticleSet*>& WW, SpaceWarp& Warp,
                  vector<QMCHamiltonian*>& h, vector<TrialWaveFunction*>& psi,
                  RealType tau, vector<RealType>& Norm, bool require_register=false);
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: MultipleEnergyEstimator.h 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
