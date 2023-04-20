//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003- by Jeongnim Kim
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
#include <Configuration.h>
#include <Message/OpenMP.h>
#include <OhmmsData/AttributeSet.h>
#include <OhmmsApp/RandomNumberControl.h>
#include <Utilities/RandomGeneratorIO.h>
#include <Utilities/Timer.h>
#include <HDFVersion.h>
#include <io/hdf_archive.h>
#include <mpi/collectives.h>
#if defined(HAVE_LIBBOOST)
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>
#endif
#include <Utilities/SimpleParser.h>

namespace qmcplusplus
{
///initialize the static data members
PrimeNumberSet<RandomGenerator_t::uint_type> RandomNumberControl::PrimeNumbers;
std::vector<RandomGenerator_t*>  RandomNumberControl::Children;
RandomGenerator_t::uint_type RandomNumberControl::Offset=11u;

/// constructors and destructors
RandomNumberControl::RandomNumberControl(const char* aname)
  :OhmmsElementBase(aname), NeverBeenInitialized(true), myCur(NULL)//, Offset(5)
{ }

/// generic output
bool RandomNumberControl::get(std::ostream& os) const
{
  if(omp_get_max_threads()>1)
  {
    for(int ip=0; ip<omp_get_max_threads(); ip++)
    {
      Children[ip]->write(os);
      os << endl;
    }
  }
  else
  {
    Random.write(os);
  }
  return true;
}

/// generic input
bool RandomNumberControl::put(std::istream& is)
{
  return true;
}

/// reset the generator
void RandomNumberControl::reset()
{
  make_seeds();
}

/// reset the generator
void RandomNumberControl::make_seeds()
{
  int pid = OHMMS::Controller->rank();
  int nprocs = OHMMS::Controller->size();
  uint_type iseed=static_cast<uint_type>(std::time(0))%1024;
  mpi::bcast(*OHMMS::Controller,iseed);
  //OHMMS::Controller->bcast(iseed);//broadcast the seed
  Offset=iseed;
  vector<uint_type> mySeeds;
  RandomNumberControl::PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);
  Random.init(pid,nprocs,mySeeds[pid],Offset+pid);
  //change children as well
  make_children();
}

void RandomNumberControl::make_children()
{
  int nthreads=omp_get_max_threads();
  int n=nthreads-Children.size();
  while(n)
  {
    Children.push_back(new RandomGenerator_t);
    n--;
  }
  int rank=OHMMS::Controller->rank();
  int nprocs=OHMMS::Controller->size();
  int baseoffset=Offset+nprocs+nthreads*rank;
  vector<uint_type> myprimes;
  PrimeNumbers.get(baseoffset,nthreads,myprimes);
  for(int ip=0; ip<nthreads; ip++)
  {
    int offset=baseoffset+ip;
    Children[ip]->init(rank,nprocs,myprimes[ip],offset);
  }
  if(nprocs<4)
  {
    ostringstream o;
    o << "  Random seeds Node = " << rank << ":";
    for(int ip=0; ip<nthreads; ip++)
      o << setw(12) << myprimes[ip];
    cout << o.str() << endl;
  }
}

xmlNodePtr
RandomNumberControl::initialize(xmlXPathContextPtr acontext)
{
  xmlXPathObjectPtr rg_request
    = xmlXPathEvalExpression((const xmlChar*)"//random",acontext);
  if(xmlXPathNodeSetIsEmpty(rg_request->nodesetval))
    put(NULL);
  else
    put(rg_request->nodesetval->nodeTab[0]);
  xmlXPathFreeObject(rg_request);
  return myCur;
}

void RandomNumberControl::test()
{
  /* Add random number generator tester
  */
  int nthreads=omp_get_max_threads();
  vector<double> avg(nthreads),avg2(nthreads);
  #pragma omp parallel for
  for(int ip=0; ip<nthreads; ++ip)
  {
    const int n=1000000;
    double sum=0.0, sum2=0.0;
    RandomGenerator_t& myrand(*Children[ip]);
    for(int i=0; i<n; ++i)
    {
      double r=myrand.rand();
      sum +=r;
      sum2+= r*r;
    }
    avg[ip]=sum/static_cast<double>(n);
    avg2[ip]=sum2/static_cast<double>(n);
  }
  vector<double> avg_tot(nthreads*OHMMS::Controller->size()),avg2_tot(nthreads*OHMMS::Controller->size());
  mpi::gather(*OHMMS::Controller,avg,avg_tot);
  mpi::gather(*OHMMS::Controller,avg2,avg2_tot);
  double avg_g=0.0;
  double avg2_g=0.0;
  for(int i=0,ii=0; i<OHMMS::Controller->size(); ++i)
  {
    for(int ip=0; ip<nthreads; ++ip,++ii)
    {
      app_log() << "RNGTest " << setw(4) << i << setw(4) << ip
                << setw(20) << avg_tot[ii] << setw(20) << avg2_tot[ii]-avg_tot[ii]*avg_tot[ii] << endl;
      avg_g+=avg_tot[ii];
      avg2_g+=avg2_tot[ii];
    }
  }
  avg_g/=static_cast<double>(nthreads*OHMMS::Controller->size());
  avg2_g/=static_cast<double>(nthreads*OHMMS::Controller->size());
  app_log() << "RNGTest " << setw(4) << OHMMS::Controller->size() << setw(4) << nthreads
            << setw(20) << avg_g << setw(20) << avg2_g-avg_g*avg_g<< endl;
  app_log().flush();
}

bool RandomNumberControl::put(xmlNodePtr cur)
{
  if(NeverBeenInitialized)
  {
    bool init_mpi = true;
    int offset_in = -1; // default is to generate by Wall-clock
    if(cur != NULL)
    {
      std::string pname("yes");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(pname,"parallel");
      oAttrib.add(offset_in,"seed");
      oAttrib.put(cur);
      if(pname == "0" || pname == "false" || pname == "no")
        init_mpi=false;
    }
    int nprocs = 1;
    int pid = 0;
    if(init_mpi)
    {
      pid = OHMMS::Controller->rank();
      nprocs = OHMMS::Controller->size();
    }
    if(offset_in<0)
    {
      offset_in=static_cast<int>(static_cast<uint_type>(std::time(0))%1024);
      app_log() << "  Offset for the random number seeds based on time " << offset_in << endl;
      mpi::bcast(*OHMMS::Controller,offset_in);
    }
    else
      offset_in%=1024;
    Offset=offset_in;
    vector<uint_type> mySeeds;
    //allocate twice of what is required
    PrimeNumbers.get(Offset,nprocs*(omp_get_max_threads()+2), mySeeds);
    Random.init(pid,nprocs,mySeeds[pid],Offset+pid);
    app_log() << "  Random number offset = " << Offset
              <<  "  seeds = " << mySeeds[0] <<"-" << mySeeds[nprocs*omp_get_max_threads()] <<  endl;
    if(nprocs<4)
    {
      int imax=8*(mySeeds.size()/8);
      int jmax=std::min(std::size_t(8),mySeeds.size());
      for(int i=0; i<imax;)
      {
        for(int j=0; j<jmax; j++, i++)
          app_log() <<  std::setw(12) << mySeeds[i];
        app_log() << endl;
      }
      for(int i=imax; i<mySeeds.size(); i++)
        app_log() <<  std::setw(12) << mySeeds[i];
      app_log() << endl;
    }
    make_children();
    NeverBeenInitialized = false;
  }
  else
    reset();
  return true;
}

void RandomNumberControl::read(const string& fname, Communicate* comm)
{
  int nthreads=omp_get_max_threads();
  vector<uint_type> vt_tot, vt;
  vector<int> shape(2,0),shape_now(2,0);
  shape_now[0]=comm->size()*nthreads;
  shape_now[1]=Random.state_size();

  if(comm->rank()==0)
  {
#if defined(HAVE_LIBBOOST)
    using boost::property_tree::ptree;
    ptree pt;
    string xname=fname+".random.xml";
    read_xml(xname, pt);
    if(!pt.empty())
    {
      string engname=pt.get<string>("random.engine");
      if(engname==Random.EngineName)
      {
        istringstream dims(pt.get<string>("random.dims"));
        dims >> shape[0] >> shape[1];
        if(shape[0]==shape_now[0] && shape[1]==shape_now[1])
        {
          vt_tot.resize(shape[0]*shape[1]);
          istringstream v(pt.get<string>("random.states"));
          for(int i=0; i<vt_tot.size(); ++i) v>>vt_tot[i];
        }
        else
          shape[0]=shape[1]=0;
      }
    }
#else
    TinyVector<hsize_t,2> shape_t(0);
    shape_t[1]=Random.state_size();
    hyperslab_proxy<vector<uint_type>,2> slab(vt_tot,shape_t);
    string h5name=fname+".random.h5";
    hdf_archive hout(comm);
    hout.open(h5name,H5F_ACC_RDONLY);
    hout.push(hdf::main_state);
    hout.push("random");
    string engname;
    hout.read(slab,Random.EngineName);
    shape[0]=static_cast<int>(slab.size(0));
    shape[1]=static_cast<int>(slab.size(1));
#endif
  }

  mpi::bcast(*comm,shape);

  if(shape[0]!=shape_now[0] || shape[1] != shape_now[1])
  {
    app_log() << "Mismatched random number generators."
      << "\n  Number of streams     : old=" << shape[0] << " new= " << comm->size()*nthreads
      << "\n  State size per stream : old=" << shape[1] << " new= " << Random.state_size()
      << "\n  Using the random streams generated at the initialization." << endl;
    return;
  }

  app_log() << "  Restart from the random number streams from the previous configuration." << endl;
  vt.resize(nthreads*Random.state_size());

  if(comm->size()>1)
    mpi::scatter(*comm,vt_tot,vt);
  else
    std::copy(vt_tot.begin(),vt_tot.end(),vt.begin());

  {
    if(nthreads>1)
    {
      vector<uint_type>::iterator vt_it(vt.begin());
      for(int ip=0; ip<nthreads; ip++, vt_it += shape[1])
      {
        vector<uint_type> c(vt_it,vt_it+shape[1]);
        Children[ip]->load(c);
      }
    }
    else
      Random.load(vt);
  }
}

void RandomNumberControl::write(const string& fname, Communicate* comm)
{
  int nthreads=omp_get_max_threads();
  vector<uint_type> vt, vt_tot;
  vt.reserve(nthreads*1024);
  if(nthreads>1)
    for(int ip=0; ip<nthreads; ++ip)
    {
      vector<uint_type> c;
      Children[ip]->save(c);
      vt.insert(vt.end(),c.begin(),c.end());
    }
  else
    Random.save(vt);
  if(comm->size()>1)
  {
    vt_tot.resize(vt.size()*comm->size());
    mpi::gather(*comm,vt,vt_tot);
  }
  else
    vt_tot=vt;

  if(comm->rank()==0)
  {
#if defined(HAVE_LIBBOOST)
    using boost::property_tree::ptree;
    ptree pt;
    ostringstream dims,vt_o;
    dims<<comm->size()*nthreads << " " << Random.state_size();
    vector<uint_type>::iterator v=vt_tot.begin();
    for(int i=0; i<comm->size()*nthreads; ++i)
    {
      std::copy(v,v+Random.state_size(),std::ostream_iterator<uint_type>(vt_o," "));
      vt_o<<endl;
      v+=Random.state_size();
    }
    pt.put("random.engine", Random.EngineName);
    pt.put("random.dims",dims.str());
    pt.put("random.states",vt_o.str());
    string xname=fname+".random.xml";
    write_xml(xname, pt);
#else
    string h5name=fname+".random.h5";
    hdf_archive hout(comm);
    hout.create(h5name);
    hout.push(hdf::main_state);
    hout.push("random");
    TinyVector<hsize_t,2> shape(comm->size()*nthreads,Random.state_size());
    hyperslab_proxy<vector<uint_type>,2> slab(vt_tot,shape);
    hout.write(slab,Random.EngineName);
    hout.close();
#endif
  }
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 6100 $   $Date: 2013-12-04 11:42:14 -0500 (Wed, 04 Dec 2013) $
 * $Id: RandomNumberControl.cpp 6100 2013-12-04 16:42:14Z jnkim $
 ***************************************************************************/
