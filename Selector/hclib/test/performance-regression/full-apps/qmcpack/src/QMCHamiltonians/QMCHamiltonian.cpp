/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
//   Department of Physics, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/NewTimer.h"
#ifdef QMC_CUDA
#include "Particle/MCWalkerConfiguration.h"
#endif

namespace qmcplusplus
{

/** constructor
*/
QMCHamiltonian::QMCHamiltonian()
  :myIndex(0),numCollectables(0),EnableVirtualMoves(false),
   id_sample(0),weight_sample(0),position_sample(0)
{
  streaming_position = false;
}

///// copy constructor is distable by declaring it as private
//QMCHamiltonian::QMCHamiltonian(const QMCHamiltonian& qh) {}

/** destructor
 */
QMCHamiltonian::~QMCHamiltonian()
{
  //@todo clean up H and auxH
}

bool QMCHamiltonian::get(std::ostream& os) const
{
  for(int i=0; i<H.size(); i++)
  {
    os.setf(ios::left);
    os << "  " << setw(16) << H[i]->myName;
    H[i]->get(os);
    os << "\n";
  }
  return true;
}

/** add a new Hamiltonian the the list of Hamiltonians.
 * @param h an operator
 * @param aname name of h
 * @param physical if true, a physical operator
 */
void
QMCHamiltonian::addOperator(QMCHamiltonianBase* h, const string& aname, bool physical)
{
  //change UpdateMode[PHYSICAL] of h so that cloning can be done correctly
  h->UpdateMode[QMCHamiltonianBase::PHYSICAL]=physical;
  if(physical)
  {
    for(int i=0; i<H.size(); ++i)
    {
      if(H[i]->myName == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to H, physical Hamiltonian " << endl;
    h->myName=aname;
    H.push_back(h);
    string tname="Hamiltonian::"+aname;
    NewTimer *atimer=new NewTimer(tname);
    myTimers.push_back(atimer);
    TimerManager.addTimer(atimer);
  }
  else
  {
    //ignore timers for now
    for(int i=0; i<auxH.size(); ++i)
    {
      if(auxH[i]->myName == aname)
      {
        app_warning() << "QMCHamiltonian::addOperator cannot " << aname << ". The name is already used" << endl;
        return;
      }
    }
    app_log() << "  QMCHamiltonian::addOperator " << aname << " to auxH " << endl;
    h->myName=aname;
    auxH.push_back(h);
  }
  
  EnableVirtualMoves|= h->getMode(QMCHamiltonianBase::VIRTUALMOVES);
}


void QMCHamiltonian::addOperatorType(const string& name, const string& type)
{
  app_log()<<"QMCHamiltonian::addOperatorType added type "<<type<<" named "<<name<<endl;
  operator_types[name] = type;
}


const string& QMCHamiltonian::getOperatorType(const string& name)
{
  map<string,string>::iterator type = operator_types.find(name);
  if(type==operator_types.end())
    APP_ABORT("QMCHamiltonian::getOperatorType\n  operator type not found for name "+name);
  return type->second;
}

///** remove a named Hamiltonian from the list
// *@param aname the name of the Hamiltonian
// *@return true, if the request hamiltonian exists and is removed.
// */
//bool
//QMCHamiltonian::remove(const string& aname)
//{
//  return false;
//}

  void QMCHamiltonian::update_source(ParticleSet& s)
  {
    for(int i=0; i<H.size(); ++i)
      H[i]->update_source(s);
    for(int i=0; i<auxH.size(); ++i)
      auxH[i]->update_source(s);
  }
  
/** add a number of properties to the ParticleSet
 * @param P ParticleSet to which multiple columns to be added
 *
 * QMCHamiltonian can add any number of properties to a ParticleSet.
 * Hindex contains the index map to the ParticleSet::PropertyList.
 * This enables assigning the properties evaluated by each QMCHamiltonianBase
 * object to the correct property column.
 */
//void
//QMCHamiltonian::addObservables(PropertySetType& plist)
//{
//  //first add properties to Observables
//  Observables.clear();
//  for(int i=0; i<H.size(); ++i) H[i]->addObservables(Observables);
//  for(int i=0; i<auxH.size(); ++i) auxH[i]->addObservables(Observables);
//
//  myIndex=plist.add(Observables.Names[0]);
//  for(int i=1; i<Observables.size(); ++i) plist.add(Observables.Names[i]);
//
//  app_log() << "  QMCHamiltonian::add2WalkerProperty added " << Observables.size() << " data to PropertyList" << endl;
//  app_log() << "    starting Index = " << myIndex << endl;
//}
//
int QMCHamiltonian::addObservables(ParticleSet& P)
{
  //first add properties to Observables
  Observables.clear();
  //ParticleSet::mcObservables (large data, e.g. density) are accumulated while evaluations
  P.Collectables.clear();
  P.Collectables.rewind();
  for(int i=0; i<H.size(); ++i)
    H[i]->addObservables(Observables,P.Collectables);
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->addObservables(Observables,P.Collectables);
  int last_obs;
  myIndex=P.PropertyList.add(Observables.Names[0]);
  for(int i=1; i<Observables.size(); ++i)
    last_obs=P.PropertyList.add(Observables.Names[i]);
  numCollectables=P.Collectables.size();
  app_log() << "\n  QMCHamiltonian::add2WalkerProperty added"
            << "\n    " << Observables.size()  << " to P::PropertyList "
            << "\n    " <<  P.Collectables.size() << " to P::Collectables "
            << "\n    starting Index of the observables in P::PropertyList = " << myIndex << endl;
  return Observables.size();
}

void QMCHamiltonian::resetObservables(int start, int ncollects)
{
  Observables.clear();
  BufferType collectables;
  collectables.rewind();
  for(int i=0; i<H.size(); ++i)
    H[i]->addObservables(Observables,collectables);
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->addObservables(Observables,collectables);
  if(collectables.size() != ncollects)
  {
    APP_ABORT("  QMCHamiltonian::resetObservables numCollectables != ncollects");
  }
  myIndex=start;
  numCollectables=ncollects;
}

void
QMCHamiltonian::registerObservables(vector<observable_helper*>& h5desc
                                    , hid_t gid)  const
{
  for(int i=0; i<H.size(); ++i)
    H[i]->registerObservables(h5desc,gid);
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->registerObservables(h5desc,gid);
}

void
QMCHamiltonian::registerCollectables(vector<observable_helper*>& h5desc
                                     , hid_t gid)  const
{
  //The physical operators cannot add to collectables
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->registerCollectables(h5desc,gid);
}


void QMCHamiltonian::initialize_traces(TraceManager& tm,ParticleSet& P)
{
  static bool first_init = true;
  if(first_init)
    app_log()<<"\n  Hamiltonian is initializing traces"<<endl;

  //fill string vectors for combined trace quantities
  vector<string> Eloc;
  vector<string> Vloc;
  vector<string> Vq,Vc,Vqq,Vqc,Vcc;
  for(int i=0; i<H.size(); ++i)
    Eloc.push_back(H[i]->myName);
  for(int i=1; i<H.size(); ++i)
    Vloc.push_back(H[i]->myName);
  for(int i=1; i<H.size(); ++i)
  {
    QMCHamiltonianBase& h = *H[i];
    if(h.is_quantum())
      Vq.push_back(h.myName);
    else if(h.is_classical())
      Vc.push_back(h.myName);
    else if(h.is_quantum_quantum())
      Vqq.push_back(h.myName);
    else if(h.is_quantum_classical())
      Vqc.push_back(h.myName);
    else if(h.is_classical_classical())
      Vcc.push_back(h.myName);
    else
      app_log()<<"  warning: potential named "<<h.myName<<" has not been classified according to its quantum domain (q,c,qq,qc,cc)\n    estimators depending on this classification may not function properly"<<endl;    
  }


  //make trace quantities available
  request.contribute_scalar("id",true);      //default trace quantity
  request.contribute_scalar("step",true);    //default trace quantity
  request.contribute_scalar("weight",true);  //default trace quantity
  request.contribute_array("position");
  for(int i=0; i<H.size(); ++i)
    H[i]->contribute_trace_quantities();
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->contribute_trace_quantities();
  

  //note availability of combined quantities
  request.contribute_combined("LocalEnergy",Eloc,true);
  request.contribute_combined("LocalPotential",Vloc,true,true);
  if(Vq.size()>0)
    request.contribute_combined("Vq",Vq,true,true);
  if(Vc.size()>0)
    request.contribute_combined("Vc",Vc,true,true);
  if(Vqq.size()>0)
    request.contribute_combined("Vqq",Vqq,true,true);
  if(Vqc.size()>0)
    request.contribute_combined("Vqc",Vqc,true,true);
  if(Vcc.size()>0)
    request.contribute_combined("Vcc",Vcc,true,true);


  //collect trace requests
  vector<TraceRequest*> requests;
  //  Hamiltonian request (id, step, weight, positions)
  requests.push_back(&request);
  //  requests from Hamiltonian components
  for(int i=0; i<H.size(); ++i)
    requests.push_back(&H[i]->request);
  //  requests from other observables
  for(int i=0; i<auxH.size(); ++i)
    requests.push_back(&auxH[i]->request);

  //collect trace quantity availability/requests from contributors/requestors
  for(int i=0;i<requests.size();++i)
    tm.request.incorporate(*requests[i]);

  //balance requests with availability, mark quantities as streaming/writing
  tm.request.determine_stream_write();

  //relay updated streaming information to all contributors/requestors
  for(int i=0;i<requests.size();++i)
    tm.request.relay_stream_info(*requests[i]);

  //set streaming/writing traces in general
  tm.update_status();

  // setup traces, if any quantities should be streaming
  bool tracing = request.streaming();
  if(tracing!=tm.streaming_traces)
    APP_ABORT("QMCHamiltonian::initialize_traces  trace request failed to initialize properly");
  if(!tracing)
  {
    if(first_init)
      app_log()<<"    no traces streaming"<<endl;
  }
  else
  {
    if(first_init)
      app_log()<<"    traces streaming"<<endl;
    //checkout trace quantities 
    //(requested sources checkout arrays to place samples in for streaming)
    //  checkout walker trace quantities
    streaming_position = request.streaming_array("position");
    if(request.streaming_default_scalars)
    {
      id_sample     = tm.checkout_int<1>("id");
      step_sample   = tm.checkout_int<1>("step");
      weight_sample = tm.checkout_real<1>("weight");
    }
    if(streaming_position)
      position_sample = tm.checkout_real<2>("position",P,DIM);
    //  checkout observable trace quantities
    for(int i=0; i<H.size(); ++i)
    {
      if(first_init)
        app_log()<<"    QMCHamiltonianBase::checkout_trace_quantities  "<<H[i]->myName<<endl;
      H[i]->checkout_trace_quantities(tm);
    }
    for(int i=0; i<auxH.size(); ++i)
    {
      if(first_init)
        app_log()<<"    QMCHamiltonianBase::checkout_trace_quantities  "<<auxH[i]->myName<<endl;
      auxH[i]->checkout_trace_quantities(tm);
    }
    //setup combined traces that depend on H information
    //  LocalEnergy, LocalPotential, Vq, Vc, Vqq, Vqc, Vcc
    if(Vloc.size()>0 && request.streaming("LocalPotential"))
      tm.make_combined_trace("LocalPotential",Vloc);
    if(Eloc.size()>0 && request.streaming("LocalEnergy"))
      tm.make_combined_trace("LocalEnergy",Eloc);
    if(Vq.size()>0 && request.streaming("Vq"))
      tm.make_combined_trace("Vq",Eloc);
    if(Vc.size()>0 && request.streaming("Vc"))
      tm.make_combined_trace("Vc",Eloc);
    if(Vqq.size()>0 && request.streaming("Vqq"))
      tm.make_combined_trace("Vqq",Eloc);
    if(Vqc.size()>0 && request.streaming("Vqc"))
      tm.make_combined_trace("Vqc",Eloc);
    if(Vcc.size()>0 && request.streaming("Vcc"))
      tm.make_combined_trace("Vcc",Eloc);

    //all trace samples have been created (streaming instances)
    //  mark the ones that will be writing also
    tm.screen_writes();

    //observables that depend on traces check them out
    if(first_init)
      app_log()<<"\n  Hamiltonian is fulfilling trace requests from observables"<<endl;
    for(int i=0; i<auxH.size(); ++i)
    {
      if(first_init)
        app_log()<<"    QMCHamiltonianBase::get_required_traces  "<<auxH[i]->myName<<endl;
      auxH[i]->get_required_traces(tm);
    }
    //report
  
    //write traces status to the log
    if(first_init)
      tm.user_report();

    first_init=false;
  }
}


void QMCHamiltonian::collect_walker_traces(Walker_t& walker,int step)
{
  if(request.streaming_default_scalars)
  {
    (*id_sample)(0)     = walker.ID;
    (*step_sample)(0)   = step;
    (*weight_sample)(0) = walker.Weight;
  }
  if(streaming_position)
    for(int i=0; i<walker.R.size(); ++i)
      for(int d=0; d<DIM; ++d)
        (*position_sample)(i,d) = walker.R[i][d];
}


void QMCHamiltonian::finalize_traces()
{
  if(request.streaming_default_scalars)
  {
    delete id_sample;
    delete step_sample;
    delete weight_sample;
  }
  if(streaming_position)
    delete position_sample;
  if(request.streaming())
  {
    for(int i=0; i<H.size(); ++i)
      H[i]->delete_trace_quantities();
    for(int i=0; i<auxH.size(); ++i)
      auxH[i]->delete_trace_quantities();
  }
  streaming_position = false;
}


/** Evaluate all the Hamiltonians for the N-particle  configuration
 *@param P input configuration containing N particles
 *@return the local energy
 */
QMCHamiltonian::Return_t
QMCHamiltonian::evaluate(ParticleSet& P)
{
  //if(EnableVirtualMoves) P.initVirtualMoves();
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    LocalEnergy += H[i]->evaluate(P);
    H[i]->setObservables(Observables);
    H[i]->collect_scalar_traces();
    myTimers[i]->stop();
    H[i]->setParticlePropertyList(P.PropertyList,myIndex);
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  // auxHevaluate(P);
  return LocalEnergy;
}
    
QMCHamiltonian::RealType 
QMCHamiltonian::evaluateValueAndDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    vector<RealType>& dlogpsi,
    vector<RealType>& dhpsioverpsi,
    bool compute_deriv)
{
  //if(EnableVirtualMoves) P.initVirtualMoves();

  LocalEnergy=KineticEnergy=H[0]->evaluate(P);
  if(compute_deriv)
    for(int i=1; i<H.size(); ++i)
      LocalEnergy += H[i]->evaluateValueAndDerivatives(P,optvars,dlogpsi,dhpsioverpsi);
  else
    for(int i=1; i<H.size(); ++i)
      LocalEnergy += H[i]->evaluate(P);
  return LocalEnergy;
}

QMCHamiltonian::RealType 
QMCHamiltonian::evaluateVariableEnergy(ParticleSet& P, bool free_nlpp)
{
  //if(EnableVirtualMoves) P.initVirtualMoves();
  RealType nlpp=0.0;
  RealType ke=H[0]->evaluate(P);
  if(free_nlpp)
    for(int i=1; i<H.size(); ++i)
    {
      if(H[i]->isNonLocal()) nlpp+=H[i]->evaluate(P);
    }
  return ke+nlpp;
}

void QMCHamiltonian::auxHevaluate(ParticleSet& P )
{
  for(int i=0; i<auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
    auxH[i]->collect_scalar_traces();
    auxH[i]->setParticlePropertyList(P.PropertyList,myIndex);
    //H[i]->setParticlePropertyList(P.PropertyList,myIndex);
  }
}

///This is more efficient. Only calculate auxH elements if moves are accepted.
void QMCHamiltonian::auxHevaluate(ParticleSet& P, Walker_t& ThisWalker)
{
  collect_walker_traces(ThisWalker,P.current_step);
  for(int i=0; i<auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
    auxH[i]->collect_scalar_traces();
    auxH[i]->setParticlePropertyList(P.PropertyList,myIndex);
  }
}

void QMCHamiltonian::rejectedMove(ParticleSet& P, Walker_t& ThisWalker )
{
  // weight should be 0 from DMC
  //   since other traced properties will be invalid
  //   (they will be from the walker moved before this one)
  collect_walker_traces(ThisWalker,P.current_step);
//   ThisWalker.rejectedMove();
  for(int i=0; i<auxH.size(); ++i)
  {
    auxH[i]->setHistories(ThisWalker);
    RealType sink = auxH[i]->rejectedMove(P);
  }
}

QMCHamiltonian::Return_t
QMCHamiltonian::evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
{
  //if(EnableVirtualMoves) P.initVirtualMoves();
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    LocalEnergy += H[i]->evaluate(P,Txy);
    H[i]->setObservables(Observables);
    H[i]->collect_scalar_traces();
    myTimers[i]->stop();
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
//   auxHevaluate(P);
  return LocalEnergy;
}


QMCHamiltonian::Return_t QMCHamiltonian::getEnsembleAverage()
{
  Return_t sum=0.0;
  for(int i=0; i<H.size(); i++)
    sum += H[i]->getEnsembleAverage();
  return sum;
}

/** return pointer to the QMCHamtiltonian with the name
 *@param aname the name of Hamiltonian
 *@return the pointer to the named term.
 *
 * If not found, return 0
 */
QMCHamiltonianBase* QMCHamiltonian::getHamiltonian(const string& aname)
{
  for(int i=0; i<H.size(); ++i)
    if(H[i]->myName == aname)
      return H[i];
  for(int i=0; i<auxH.size(); ++i)
    if(auxH[i]->myName == aname)
      return auxH[i];
  return 0;
}

void QMCHamiltonian::resetTargetParticleSet(ParticleSet& P)
{
  for(int i=0; i<H.size(); i++)
    H[i]->resetTargetParticleSet(P);
  for(int i=0; i<auxH.size(); i++)
    auxH[i]->resetTargetParticleSet(P);
}

void QMCHamiltonian::setRandomGenerator(RandomGenerator_t* rng)
{
  for(int i=0; i<H.size(); i++)
    H[i]->setRandomGenerator(rng);
  for(int i=0; i<auxH.size(); i++)
    auxH[i]->setRandomGenerator(rng);
}

QMCHamiltonian* QMCHamiltonian::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  QMCHamiltonian* myclone=new QMCHamiltonian;
  for(int i=0; i<H.size(); ++i)
    H[i]->add2Hamiltonian(qp,psi,*myclone);
  for(int i=0; i<auxH.size(); ++i)
    auxH[i]->add2Hamiltonian(qp,psi,*myclone);
  //sync indices
  myclone->resetObservables(myIndex,numCollectables);
  //Hamiltonian needs to make sure qp.Collectables are the same as defined by the original Hamiltonian
  if(numCollectables)
  {
    qp.Collectables.clear();
    qp.Collectables.resize(numCollectables);
  }
  //Assume tau is correct for the Kinetic energy operator and assign to the rest of the clones
  //Return_t tau = H[0]->Tau;
  //myclone->setTau(tau);
  return myclone;
}

QMCHamiltonian::Return_t QMCHamiltonian::registerData(ParticleSet& P, BufferType& buffer)
{
  LocalEnergy=0.0;
  for(int i=0; i<H.size(); ++i)
  {
    LocalEnergy+=H[i]->registerData(P,buffer);
    H[i]->setObservables(Observables);
  }
  buffer.add(LocalEnergy);
  return LocalEnergy;
}

QMCHamiltonian::Return_t QMCHamiltonian::updateBuffer(ParticleSet& P, BufferType& buffer)
{
  LocalEnergy=0.0;
  for(int i=0; i<H.size(); ++i)
  {
    LocalEnergy+=H[i]->updateBuffer(P,buffer);
    H[i]->setObservables(Observables);
  }
  buffer.add(LocalEnergy);
  return LocalEnergy;
}

void QMCHamiltonian::copyFromBuffer(ParticleSet& P, BufferType& buffer)
{
  for(int i=0; i<H.size(); ++i)
    H[i]->copyFromBuffer(P,buffer);
  buffer.get(LocalEnergy);
}

QMCHamiltonian::Return_t QMCHamiltonian::evaluate(ParticleSet& P, BufferType& buffer)
{
  LocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    LocalEnergy += H[i]->Value;
    H[i]->setObservables(Observables);
    H[i]->collect_scalar_traces();
  }
  KineticEnergy=H[0]->Value;
  P.PropertyList[LOCALENERGY]=LocalEnergy;
  P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for(int i=0; i<auxH.size(); ++i)
  {
    RealType sink = auxH[i]->evaluate(P);
    auxH[i]->setObservables(Observables);
    auxH[i]->collect_scalar_traces();
  }
  for(int i=0; i<H.size(); ++i)
    H[i]->copyToBuffer(P,buffer);
  buffer.put(LocalEnergy);
  return LocalEnergy;
}

QMCHamiltonian::Return_t QMCHamiltonian::evaluatePbyP(ParticleSet& P, int active)
{
  NewLocalEnergy = 0.0;
  for(int i=0; i<H.size(); ++i)
    NewLocalEnergy +=H[i]->evaluatePbyP(P,active);
  return NewLocalEnergy;
}
void QMCHamiltonian::acceptMove(int active)
{
  for(int i=0; i<H.size(); ++i)
    H[i]->acceptMove(active);
  LocalEnergy=NewLocalEnergy;
  for(int i=0; i<H.size(); ++i)
    H[i]->setObservables(Observables);
}

void QMCHamiltonian::rejectMove(int active)
{
  for(int i=0; i<H.size(); ++i)
    H[i]->rejectMove(active);
}


#ifdef QMC_CUDA
void
QMCHamiltonian::evaluate(MCWalkerConfiguration &W,
                         vector<RealType> &energyVector)
{
  vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  if (LocalEnergyVector.size() != nw)
  {
    LocalEnergyVector.resize(nw);
    KineticEnergyVector.resize(nw);
    AuxEnergyVector.resize(nw);
  }
  if (energyVector.size() != nw)
    energyVector.resize(nw);
  for (int i=0; i<LocalEnergyVector.size(); i++)
    LocalEnergyVector[i] = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    H[i]->addEnergy(W, LocalEnergyVector);
    //H[i]->setObservables(Observables);
    myTimers[i]->stop();
  }
  //KineticEnergyVector=H[0]->ValueVector;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[LOCALENERGY] = LocalEnergyVector[iw];
    walkers[iw]->getPropertyBase()[LOCALPOTENTIAL] =
      LocalEnergyVector[iw] - walkers[iw]->getPropertyBase()[NUMPROPERTIES];
  }
  energyVector = LocalEnergyVector;
  // P.PropertyList[LOCALENERGY]=LocalEnergy;
  // P.PropertyList[LOCALPOTENTIAL]=LocalEnergy-KineticEnergy;
  for(int i=0; i<auxH.size(); ++i)
  {
    auxH[i]->addEnergy(W, AuxEnergyVector);
    //auxH[i]->setObservables(Observables);
  }
}



void
QMCHamiltonian::evaluate(MCWalkerConfiguration &W,
                         vector<RealType> &energyVector,
                         vector<vector<NonLocalData> > &Txy)
{
  vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  if (LocalEnergyVector.size() != nw)
  {
    LocalEnergyVector.resize(nw);
    KineticEnergyVector.resize(nw);
    AuxEnergyVector.resize(nw);
  }
  if (energyVector.size() != nw)
    energyVector.resize(nw);
  std::fill(LocalEnergyVector.begin(),LocalEnergyVector.end(),0.0);
  //for (int i=0; i<LocalEnergyVector.size(); i++)
  //  LocalEnergyVector[i] = 0.0;
  for(int i=0; i<H.size(); ++i)
  {
    myTimers[i]->start();
    H[i]->addEnergy(W, LocalEnergyVector, Txy);
    myTimers[i]->stop();
  }
  KineticEnergyVector=H[0]->ValueVector;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    walkers[iw]->getPropertyBase()[LOCALENERGY] = LocalEnergyVector[iw];
    walkers[iw]->getPropertyBase()[LOCALPOTENTIAL] =
      LocalEnergyVector[iw] - walkers[iw]->getPropertyBase()[NUMPROPERTIES];
  }
  energyVector = LocalEnergyVector;

  if(auxH.size())
  {
    std::fill(AuxEnergyVector.begin(),AuxEnergyVector.end(),0.0);
    for(int i=0; i<auxH.size(); ++i)
      auxH[i]->addEnergy(W, AuxEnergyVector);
  }
}
#else
void
QMCHamiltonian::evaluate(MCWalkerConfiguration &W,
                         vector<RealType> &energyVector)
{
}

void
QMCHamiltonian::evaluate(MCWalkerConfiguration &W,
                         vector<RealType> &energyVector,
                         vector<vector<NonLocalData> > &Txy)
{
}
#endif
}



/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 6366 $   $Date: 2014-10-03 15:43:46 -0400 (Fri, 03 Oct 2014) $
 * $Id: QMCHamiltonian.cpp 6366 2014-10-03 19:43:46Z jtkrogel $
 ***************************************************************************/

