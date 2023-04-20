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
#include "QMCDrivers/DMC/WalkerReconfiguration.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/UtilityFunctions.h"
#include "Utilities/RandomGenerator.h"
using namespace qmcplusplus;

/** default constructor
 *
 * set SwapMode
 */
WalkerReconfiguration::WalkerReconfiguration(Communicate* c) :WalkerControlBase(c)
{
  SwapMode=1;
  UnitZeta=Random();
  //ofstream fout("check.dat");
}

int WalkerReconfiguration::getIndexPermutation(MCWalkerConfiguration& W)
{
  int nw(W.getActiveWalkers());
  if(Zeta.size()!=nw)
  {
    Zeta.resize(nw+1);
    IndexCopy.resize(nw);
    wConf.resize(nw);
  }
  //accumulate the energies
  RealType esum=0.0,e2sum=0.0,wtot=0.0,ecum=0.0;
  MCWalkerConfiguration::iterator it(W.begin());
  RealType r2_accepted=0.0,r2_proposed=0.0;
  for(int iw=0; iw<nw; iw++)
  {
    r2_accepted+=(*it)->Properties(R2ACCEPTED);
    r2_proposed+=(*it)->Properties(R2PROPOSED);
    RealType wgt((*it)->Weight);
    RealType e((*it)->Properties(LOCALENERGY));
    esum += wgt*e;
    e2sum += wgt*e*e;
    ecum += e;
    wtot += wConf[iw]=wgt;
    ++it;
  }
  curData[ENERGY_INDEX]=esum;
  curData[ENERGY_SQ_INDEX]=e2sum;
  curData[WALKERSIZE_INDEX]=nw;
  curData[WEIGHT_INDEX]=wtot;
  curData[EREF_INDEX]=ecum;
  curData[R2ACCEPTED_INDEX]=r2_accepted;
  curData[R2PROPOSED_INDEX]=r2_proposed;
  RealType nwInv=1.0/static_cast<RealType>(nw);
  RealType dstep=UnitZeta*nwInv;
  for(int iw=0; iw<nw; iw++)
  {
    Zeta[iw]=wtot*(dstep+static_cast<RealType>(iw)*nwInv);
  }
  Zeta[nw]=wtot+1.0;
  //for(int iw=0; iw<nw; iw++) {
  //  fout << iw << " " << Zeta[iw+1]-Zeta[iw] << " " << wConf[iw] << endl;
  //}
  //assign negative
  //std::fill(IndexCopy.begin(),IndexCopy.end(),-1);
  int ind=0;
  RealType wCur=0.0;
  //surviving walkers
  int icdiff=0;
  it=W.begin();
  vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++)
  {
    RealType tryp=wCur+abs(wConf[iw]);
    int ni=0;
    while(Zeta[ind]<tryp && Zeta[ind] >= wCur)
    {
      //IndexCopy[ind]=iw;
      ind++;
      ni++;
    }
    wCur+=abs(wConf[iw]);
    if(ni)
    {
      icdiff++;
    }
    ipip[iw]=ni;
  }
  //ofstream fout("check.dat", ios::app);
  //fout << wtot << " " << icdiff << endl;
  vector<int> plus,minus;
  for(int iw=0; iw<nw; iw++)
  {
    int m=ipip[iw];
    if(m>1)
      plus.insert(plus.end(),m-1,iw);
    else
      if(m==0)
        minus.push_back(iw);
  }
  curData[FNSIZE_INDEX]=nw-minus.size();
  curData[RNONESIZE_INDEX]=minus.size();
  for(int i=0; i<plus.size(); i++)
  {
    int im=minus[i],ip=plus[i];
    W[im]->makeCopy(*(W[ip]));
    W[im]->ParentID=W[ip]->ID;
    W[im]->ID=(++NumWalkersCreated)*NumContexts+MyContext;
  }
  //int killed = shuffleIndex(nw);
  //fout << "# Total weight " << wtot << " " << killed <<  endl;
  //cout << "<<<< CopyIndex " << endl;
  //std::copy(IndexCopy.begin(), IndexCopy.end(), ostream_iterator<int>(cout, " "));
  //cout << endl << "<<<<<<" << endl;
  //for(int iw=0; iw<nw; iw++) {
  //  if(IndexCopy[iw] != iw) {
  //    W[iw]->assign(*(W[IndexCopy[iw]]));
  //  }
  //}
  return icdiff;
}

int WalkerReconfiguration::shuffleIndex(int nw)
{
  vector<int> ipip(nw,0);
  for(int iw=0; iw<nw; iw++)
    ipip[IndexCopy[iw]]+=1;
  vector<int> indz;
  for(int iw=0; iw<nw; iw++)
  {
    if(ipip[iw]==0)
    {
      indz.push_back(iw);
    }
  }
  int ikilled=0;
  for(int iw=0; iw<nw; iw++)
  {
    if(ipip[iw] != 0)
    {
      IndexCopy[iw]=iw;
      for(int i=1; i<ipip[iw]; i++)
      {
        IndexCopy[indz[ikilled++]]=iw;
      }
    }
  }
  return indz.size();
}

int
WalkerReconfiguration::branch(int iter, MCWalkerConfiguration& W, RealType trigger)
{
  int nwkept = getIndexPermutation(W);
  //update EnsembleProperty
  measureProperties(iter);
  W.EnsembleProperty=EnsembleProperty;
  //W.EnsembleProperty.NumSamples=curData[WALKERSIZE_INDEX];
  //W.EnsembleProperty.Weight=curData[WEIGHT_INDEX];
  //RealType wgtInv(1.0/curData[WEIGHT_INDEX]);
  //accumData[ENERGY_INDEX]     += curData[ENERGY_INDEX]*wgtInv;
  //accumData[ENERGY_SQ_INDEX]  += curData[ENERGY_SQ_INDEX]*wgtInv;
  //accumData[WALKERSIZE_INDEX] += nwkept;
  ////accumData[WALKERSIZE_INDEX] += curData[WALKERSIZE_INDEX];
  //accumData[WEIGHT_INDEX]     += curData[WEIGHT_INDEX];
  //set Weight and Multiplicity to default values
  MCWalkerConfiguration::iterator it(W.begin()),it_end(W.end());
  while(it != it_end)
  {
    (*it)->Weight= 1.0;
    (*it)->Multiplicity=1.0;
    ++it;
  }
  //curData[WALKERSIZE_INDEX]=nwkept;
  return nwkept;
}

/***************************************************************************
 * $RCSfile: WalkerReconfiguration.cpp,v $   $Author: jnkim $
 * $Revision: 6034 $   $Date: 2013-10-28 10:23:38 -0400 (Mon, 28 Oct 2013) $
 * $Id: WalkerReconfiguration.cpp 6034 2013-10-28 14:23:38Z jnkim $
 ***************************************************************************/

