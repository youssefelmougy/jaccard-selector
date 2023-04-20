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
#include "Estimators/LocalEnergyEstimator.h"
#include <Utilities/IteratorUtility.h>

namespace qmcplusplus
{

LocalEnergyEstimator::LocalEnergyEstimator(QMCHamiltonian& h, bool use_hdf5)
  :UseHDF5(use_hdf5),refH(h)
{
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
  scalars.resize(SizeOfHamiltonians+LE_MAX);
  scalars_saved.resize(SizeOfHamiltonians+LE_MAX);
}

ScalarEstimatorBase* LocalEnergyEstimator::clone()
{
  return new LocalEnergyEstimator(*this);
}

void  LocalEnergyEstimator::registerObservables(vector<observable_helper*>& h5desc, hid_t gid)
{
  if(!UseHDF5)
    return;
  int loc=h5desc.size();
  //add LocalEnergy and LocalPotential
  h5desc.push_back(new observable_helper("LocalEnergy"));
  h5desc.push_back(new observable_helper("LocalEnergy_sq"));
  h5desc.push_back(new observable_helper("LocalPotential"));
  std::vector<int> onedim(1,1);
  h5desc[loc]->set_dimensions(onedim,FirstIndex);
  h5desc[loc++]->open(gid);
  h5desc[loc]->set_dimensions(onedim,FirstIndex+1);
  h5desc[loc++]->open(gid);
  h5desc[loc]->set_dimensions(onedim,FirstIndex+2);
  h5desc[loc++]->open(gid);
  //hamiltonian adds more
  refH.registerObservables(h5desc,gid);
  int correction=FirstIndex+3;
  for(int i=loc; i<h5desc.size(); ++i)
    h5desc[i]->lower_bound += correction;
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 * @param record storage of scalar records (name,value)
 */
void LocalEnergyEstimator::add2Record(RecordListType& record)
{
  FirstIndex = record.size();
  int dumy=record.add("LocalEnergy");
  dumy=record.add("LocalEnergy_sq");
  dumy=record.add("LocalPotential");
  for(int i=0; i<SizeOfHamiltonians; ++i)
    record.add(refH.getObservableName(i));
  LastIndex=record.size();
  clear();
}

}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: LocalEnergyEstimator.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
