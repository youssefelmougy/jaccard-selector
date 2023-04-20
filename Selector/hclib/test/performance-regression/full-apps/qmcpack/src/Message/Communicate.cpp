/////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002, 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "qmc_common.h"
#include "Message/Communicate.h"
#include "Message/TagMaker.h"
#include <iostream>
#include <cstdio>
#include <Platforms/sysutil.h>
#include <tau/profiler.h>
#include <Utilities/UtilityFunctions.h>
#include <fstream>
#include <qmc_common.h>

#ifdef HAVE_ADIOS
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include "ADIOS/ADIOS_config.h"
#endif


//static data of TagMaker::CurrentTag is initialized.
int TagMaker::CurrentTag = 1000;

//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;

//default constructor: ready for a serial execution
Communicate::Communicate():
  myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0), d_ngroups(1)
{
}


Communicate::Communicate(int argc, char **argv)
{
  initialize(argc,argv);
}

//exclusive:  OOMPI, MPI or Serial
#ifdef HAVE_OOMPI

void
Communicate::set_world()
{
  myComm = OOMPI_COMM_WORLD;
  myMPI = myComm.Get_mpi();
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();
  d_groupid=0;
  d_ngroups=1;
}


Communicate::Communicate(const Communicate& comm, int nparts)
{
  qmcplusplus::qmc_common.mpi_groups=nparts;

  //this is a workaround due to the OOMPI bug with split
  if(nparts>1)
  {
    std::vector<int> nplist(nparts+1);
    int p=FairDivideLow(comm.rank(), comm.size(), nparts, nplist); //group
    int q=comm.rank()-nplist[p];//rank within a group
    //int n=comm.size()/nparts;
    //int p=comm.rank()/n;
    //int q=comm.rank()%n;
    MPI_Comm row;
    MPI_Comm_split(comm.getMPI(),p,q,&row);
    myComm=OOMPI_Intra_comm(row);
    d_groupid=p;
  }
  else
  {
    myComm=OOMPI_Intra_comm(comm.getComm());
    d_groupid=0;
  }
  myMPI = myComm.Get_mpi();
  d_mycontext=myComm.Rank();
  d_ncontexts=myComm.Size();
  d_ngroups=nparts;
}

Communicate::Communicate(const Communicate& comm, const std::vector<int>& jobs)
{
  //this is a workaround due to the OOMPI bug with split
  if(jobs.size()>1)
  {
    MPI_Comm row;
    std::vector<int> nplist(jobs.size()+1);
    int p=-1;
    nplist[0]=0;
    for(int i=0; i<jobs.size(); ++i)
    {
      nplist[i+1]=nplist[i]+jobs[i];
      if(comm.rank()>=nplist[i] && comm.rank()<nplist[i+1])
        p=i;
    }
    if(nplist[comm.size()] != comm.size())
    {
      APP_ABORT("Communicate::Communicate(comm,jobs) Cannot group MPI tasks. Mismatch of the tasks");
    }
    int q=comm.rank()-nplist[p];//rank within a group
    MPI_Comm_split(comm.getMPI(),p,q,&row);
    myComm=OOMPI_Intra_comm(row);
    d_groupid=p;
  }
  else
  {
    myComm=OOMPI_Intra_comm(comm.getComm());
    d_groupid=0;
  }
  myMPI = myComm.Get_mpi();
  d_mycontext=myComm.Rank();
  d_ncontexts=myComm.Size();
  d_ngroups=jobs.size();
}



//================================================================
// Implements Communicate with OOMPI library
//================================================================
Communicate::~Communicate()
{
}

void Communicate::initialize(int argc, char **argv)
{
  qmcplusplus::qmc_common.initialize(argc,argv);
  OOMPI_COMM_WORLD.Init(argc, argv);
  myComm = OOMPI_COMM_WORLD;
  myMPI = myComm.Get_mpi();
  d_mycontext = OOMPI_COMM_WORLD.Rank();
  d_ncontexts = OOMPI_COMM_WORLD.Size();
  d_groupid=0;
  d_ngroups=1;
#ifdef __linux__
  for (int proc=0; proc<OHMMS::Controller->size(); proc++)
  {
    if (OHMMS::Controller->rank() == proc)
    {
      fprintf (stderr, "Rank = %4d  Free Memory = %5ld MB\n", proc, freemem());
    }
    barrier();
  }
  barrier();
#endif
  std::string when="qmc."+getDateAndTime("%Y%m%d_%H%M");
}

void Communicate::finalize()
{
  static bool has_finalized=false;

  if(!has_finalized)
  {
#ifdef HAVE_ADIOS
    if(ADIOS::get_adios_init()){
      adios_read_finalize_method(ADIOS_READ_METHOD_BP);
      adios_finalize(OHMMS::Controller->rank());
    }
#endif
    OOMPI_COMM_WORLD.Finalize();
    has_finalized=true;
  }
}

void Communicate::cleanupMessage(void*)
{
}

void Communicate::abort()
{
  OOMPI_COMM_WORLD.Abort();
}

void Communicate::barrier()
{
  myComm.Barrier();
}

void Communicate::abort(const char* msg)
{
  std::cerr << msg << std::endl;
  OOMPI_COMM_WORLD.Abort();
}


#else

Communicate::~Communicate() {}

void Communicate::initialize(int argc, char **argv)
{
  std::string when="qmc."+getDateAndTime("%Y%m%d_%H%M");
}

void Communicate::set_world()
{
}

void Communicate::finalize()
{
}

void Communicate::abort()
{
  abort();
}

void Communicate::abort(const char* msg)
{
  std::cerr << msg << std::endl;
  abort();
}

void Communicate::barrier()
{
}

void Communicate::cleanupMessage(void*)
{
}

Communicate::Communicate(const Communicate& comm, int nparts)
  : myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0)
{
}

Communicate::Communicate(const Communicate& comm, const std::vector<int>& jobs)
  : myMPI(0), d_mycontext(0), d_ncontexts(1), d_groupid(0)
{
}

#endif // !HAVE_OOMPI
/***************************************************************************
 * $RCSfile$   $Author: cynthiagu $
 * $Revision: 6315 $   $Date: 2014-06-11 12:53:45 -0400 (Wed, 11 Jun 2014) $
 * $Id: Communicate.cpp 6315 2014-06-11 16:53:45Z cynthiagu $
 ***************************************************************************/
