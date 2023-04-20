//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "QMCWaveFunctions/Jastrow/WMConstraints.h"
#include "QMCWaveFunctions/Jastrow/TruncatedPadeFunctor.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

WMConstraints::~WMConstraints()
{
  delete_iter(FuncList.begin(), FuncList.end());
}

bool WMConstraints::put(xmlNodePtr cur)
{
  //get generic parameters and grid information
  bool success=getVariables(cur);
  return success;
}

void WMConstraints::apply()
{
  for(int i=0; i<FuncList.size(); i++)
  {
    FuncList[i]->reset();
  }
}

void WMConstraints::addOptimizables(VarRegistry<RealType>& outVars)
{
  map<string,BasisSetType*>::iterator it(myBasisSet.begin()),it_end(myBasisSet.end());
  while(it != it_end)
  {
    BasisSetType& basis(*((*it).second));
    for(int ib=0; ib<basis.size(); ib++)
    {
      basis[ib]->addOptimizables(outVars);
    }
    ++it;
  }
  for(int i=0; i<InFuncList.size(); i++)
  {
    InFuncList[i]->addOptimizables(outVars);
  }
}

void WMConstraints::addBasisGroup(xmlNodePtr cur)
{
  string sourceOpt("e");
  string elementType("e");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sourceOpt,"source");
  aAttrib.add(elementType,"elementType");
  aAttrib.put(cur);
  RealType rcut(myGrid->rmax());
  map<string,BasisSetType*>::iterator it(myBasisSet.find(elementType));
  if(it == myBasisSet.end())
  {
    BasisSetType* newBasis=new BasisSetType;
    cur=cur->children;
    while(cur != NULL)
    {
      string cname((const char*)(cur->name));
      if(cname == "parameter")
      {
        //BasisType* a=new BasisType(1.0,rcut);
        WMFunctor<RealType>* a = new WMFunctor<RealType>(1.0,rcut);
        a->put(cur);
        newBasis->push_back(a);
      }
      cur=cur->next;
    }
    //add a new BasisSet
    myBasisSet[elementType]=newBasis;
  }
}


OrbitalBase* WMConstraints::createTwoBody(ParticleSet& target)
{
  InFuncType* infunc=0;
  setRadialGrid(target);
  xmlNodePtr cur=myNode->children;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "basisGroup")
    {
      addBasisGroup(cur);
    }
    else
      if(cname =="correlation")
      {
        string speciesA(target.getName());
        string speciesB(target.getName());
        OhmmsAttributeSet aAttrib;
        aAttrib.add(speciesA,"speciesA");
        aAttrib.add(speciesB,"speciesB");
        aAttrib.put(cur);
        if(speciesA==speciesB)
        {
          map<string,BasisSetType*>::iterator it(myBasisSet.find(speciesA));
          if(it == myBasisSet.end())
          {
            app_error() <<  "  WMBasisSet for " << speciesA << " does not exist." << endl;
            continue;
          }
          app_log() << "    Creating a correlation function = " << speciesA << "-" << speciesB << endl;
          infunc = createCorrelation(cur,(*it).second);
        }
      }
    cur=cur->next;
  }
  if(infunc==0)
    //we can create a default function
  {
    infunc=createDefaultTwoBody(myNode,target.getName());
  }
  //add analytic function
  InFuncList.push_back(infunc);
  //try to correct the cusp condition
  TruncatedPadeFunctor<RealType>* wrapFunc
  = new TruncatedPadeFunctor<RealType>(-0.5,infunc,myGrid->rmax());
  wrapFunc->reset();
  typedef TwoBodyJastrowOrbital<FuncType> JeeType;
  //create a Jastrow function
  JeeType *J2 = new JeeType(target);
  //numerical functor
  FuncType* nfunc= 0;
  if(wrapFunc->Rcut > 0.0)//pade truncation is good to go
  {
    //InFuncList.push_back(wrapFunc);
    nfunc= new FuncType(wrapFunc,myGrid);
  }
  else //use the original function
  {
    nfunc= new FuncType(infunc,myGrid);
    delete wrapFunc;
  }
#if !defined(HAVE_MPI)
  if(PrintTables)
  {
    ofstream fout("Jee.dat");
    nfunc->print(fout);
  }
#endif
  for(int i=0; i<4; i++)
    J2->addFunc(nfunc);
  FuncList.push_back(nfunc);
  return J2;
}

OrbitalBase* WMConstraints::createOneBody(ParticleSet& target, ParticleSet& source)
{
  vector<InFuncType*> jnSet;
  jnSet.resize(source.getSpeciesSet().getTotalNum(),0);
  xmlNodePtr cur=myNode->children;
  bool noOneBody=true;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "basisGroup")
    {
      addBasisGroup(cur);
    }
    else
      if(cname =="correlation")
      {
        string speciesA("e");
        string speciesB("e");
        OhmmsAttributeSet aAttrib;
        aAttrib.add(speciesA,"speciesA");
        aAttrib.add(speciesB,"speciesB");
        aAttrib.put(cur);
        if(speciesA != speciesB)
        {
          map<string,BasisSetType*>::iterator it(myBasisSet.find(speciesA));
          if(it == myBasisSet.end())
          {
            app_error() <<  "  WMBasisSet for " << speciesA << " does not exist." << endl;
            continue;
          }
          app_log() << "    Creating a correlation function = " << speciesA << "-" << speciesB << endl;
          int gid=source.getSpeciesSet().addSpecies(speciesA);
          jnSet[gid] = createCorrelation(cur,(*it).second);
          noOneBody=false;
        }
      }
    cur=cur->next;
  }
  if(noOneBody)
    return 0;
  typedef OneBodyJastrow<FuncType> JneType;
  JneType* jne=new JneType(source,target);
  for(int ig=0; ig<jnSet.size(); ig++)
  {
    if(jnSet[ig])
    {
      FuncType* nfunc= new FuncType(jnSet[ig],myGrid);
      jne->addFunc(ig,nfunc);
      FuncList.push_back(nfunc);
      InFuncList.push_back(jnSet[ig]);
    }
  }
  return jne;
}

WMConstraints::InFuncType*
WMConstraints::createCorrelation(xmlNodePtr cur,BasisSetType* basis)
{
  int nc=0;
  InFuncType* acombo=new InFuncType;
  cur=cur->children;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == "parameter")
    {
      string id("0");
      string ref("0");
      RealType c=1.0;
      OhmmsAttributeSet aAttrib;
      aAttrib.add(id,"id");
      aAttrib.add(ref,"ref");
      aAttrib.put(cur);
      putContent(c,cur);
      if(nc <basis->size())
        acombo->add((*basis)[nc++],c,id);
    }
    cur=cur->next;
  }
  if(nc)
    return acombo;
  else
  {
    delete acombo;
    return 0;
  }
}

WMConstraints::InFuncType*
WMConstraints::createDefaultTwoBody(xmlNodePtr cur, const string& tname)
{
  xmlNodePtr basisNode = xmlNewNode(NULL,(const xmlChar*)"basisGroup");
  xmlNewProp(basisNode,(const xmlChar*)"source",(const xmlChar*)tname.c_str());
  xmlNodePtr b0 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"4.55682");
  xmlNewProp(b0,(const xmlChar*)"id",(const xmlChar*)"eeB0");
  xmlNewProp(b0,(const xmlChar*)"name",(const xmlChar*)"B");
  xmlNodePtr b1 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"1.83679");
  xmlNewProp(b1,(const xmlChar*)"id",(const xmlChar*)"eeB1");
  xmlNewProp(b1,(const xmlChar*)"name",(const xmlChar*)"B");
  xmlNodePtr b2 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"4.59201");
  xmlNewProp(b2,(const xmlChar*)"id",(const xmlChar*)"eeB2");
  xmlNewProp(b2,(const xmlChar*)"name",(const xmlChar*)"B");
  xmlNodePtr b3 = xmlNewTextChild(basisNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0.226152");
  xmlNewProp(b3,(const xmlChar*)"id",(const xmlChar*)"eeB3");
  xmlNewProp(b3,(const xmlChar*)"name",(const xmlChar*)"B");
  xmlAddChild(cur,basisNode);
  addBasisGroup(basisNode);
  xmlNodePtr corrNode = xmlNewNode(NULL,(const xmlChar*)"correlation");
  xmlNewProp(corrNode,(const xmlChar*)"speciesA",(const xmlChar*)tname.c_str());
  xmlNewProp(corrNode,(const xmlChar*)"speciesB",(const xmlChar*)tname.c_str());
  xmlNodePtr c0 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"-0.137958");
  xmlNewProp(c0,(const xmlChar*)"id",(const xmlChar*)"eeC0");
  xmlNewProp(c0,(const xmlChar*)"name",(const xmlChar*)"C");
  xmlNodePtr c1 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"-1.09405");
  xmlNewProp(c1,(const xmlChar*)"id",(const xmlChar*)"eeC1");
  xmlNewProp(c1,(const xmlChar*)"name",(const xmlChar*)"C");
  xmlNodePtr c2 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"1.06578");
  xmlNewProp(c2,(const xmlChar*)"id",(const xmlChar*)"eeC2");
  xmlNewProp(c2,(const xmlChar*)"name",(const xmlChar*)"C");
  xmlNodePtr c3 = xmlNewTextChild(corrNode,NULL,(const xmlChar*)"parameter",(const xmlChar*)"0.96602");
  xmlNewProp(c3,(const xmlChar*)"id",(const xmlChar*)"eeC3");
  xmlNewProp(c3,(const xmlChar*)"name",(const xmlChar*)"C");
  xmlAddChild(cur,corrNode);
  map<string,BasisSetType*>::iterator it(myBasisSet.find("e"));
  return createCorrelation(corrNode,(*it).second);
}
}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5794 $   $Date: 2013-04-25 20:14:53 -0400 (Thu, 25 Apr 2013) $
 * $Id: WMConstraints.cpp 5794 2013-04-26 00:14:53Z jmcminis $
 ***************************************************************************/
