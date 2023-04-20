//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <OhmmsData/AttributeSet.h>
#include <QMCApp/ParticleSetPool.h>
#include <QMCApp/WaveFunctionPool.h>
#include <QMCApp/HamiltonianPool.h>
#include <Estimators/PostProcessor.h>
#include <Estimators/SpinDensityPostProcessor.h>

namespace qmcplusplus
{

  PostProcessor::PostProcessor(const string& id,int ss,int se)
  {
    app_log()<<"  PostProcessor Engine created"<<endl;
    set(id,ss,se);
  }


  void PostProcessor::put(xmlNodePtr cur,ParticleSetPool& ppool,
                          WaveFunctionPool& wpool,HamiltonianPool& hpool)
  {
    app_log()<<"  Initializing PostProcessor"<<endl;
    
    string snone = "";

    string Pqname  = snone;
    string Pcname  = snone;
    string Psiname = snone;
    string Hname   = snone;

    OhmmsAttributeSet attrib;
    attrib.add(Pqname ,"Pq" );
    attrib.add(Pcname ,"Pc" );
    attrib.add(Psiname,"Psi");
    attrib.add(Hname  ,"H"  );
    attrib.put(cur);

    bool notPq  =  Pqname==snone;
    bool notPc  =  Pcname==snone;
    bool notPsi = Psiname==snone;
    bool notH   =   Hname==snone;
    if(notPq||notPc||notPsi||notH)
    {
      if(notPq)
        app_log()<<"PostProcessor::put  attribute Pq  (quantum ParticleSet) is missing"<<endl;
      if(notPc)
        app_log()<<"PostProcessor::put  attribute Pc  (classical ParticleSet) is missing"<<endl;
      if(notPsi)
        app_log()<<"PostProcessor::put  attribute Psi (WaveFunction) is missing"<<endl;
      if(notH)
        app_log()<<"PostProcessor::put  attribute H   (Hamiltonian) is missing"<<endl;
      APP_ABORT("PostProcessor::put  xml input is incomplete, see messages above.");
    }

    ParticleSet&       Pq  = *ppool.getParticleSet(Pqname);
    ParticleSet&       Pc  = *ppool.getParticleSet(Pcname);
    TrialWaveFunction& Psi = *wpool.getWaveFunction(Psiname);
    QMCHamiltonian&    H   = *hpool.getHamiltonian(Hname);

    notPq  =  &Pq==0 ||  Pq.getName()!=Pqname;
    notPc  =  &Pc==0 ||  Pc.getName()!=Pcname;
    notPsi = &Psi==0 || Psi.getName()!=Psiname;
    notH   =   &H==0 ||   H.getName()!=Hname;
    if(notPq||notPc||notPsi||notH)
    {
      if(notPq)
        app_log()<<"PostProcessor::put ParticleSet Pq="<<Pqname<<" cannot be found"<<endl;
      if(notPc)
        app_log()<<"PostProcessor::put ParticleSet Pc="<<Pcname<<" cannot be found"<<endl;
      if(notPsi)
        app_log()<<"PostProcessor::put WaveFunction Psi="<<Psiname<<" cannot be found"<<endl;
      if(notH)
        app_log()<<"PostProcessor::put Hamiltonian H="<<Hname<<" cannot be found"<<endl;
      APP_ABORT("PostProcessor::put  cannot find qmcsystem objects");
    }

    cur = cur->children;
    while(cur != NULL)
    { 
      string cname((const char*)cur->name);
      bool not_text = cname!="text";
      PostProcessorBase* pp = 0;
      if(cname=="spindensity")
        pp = new SpinDensityPostProcessor(Pq,Pc,H);
      else if(not_text)
        APP_ABORT("PostProcessor::put  "+cname+" is not a valid element of postprocess");
      if(not_text)
      {
        if(pp)
        {
          app_log()<<"    adding "+cname+" postprocessor"<<endl;
          pp->set(cname,project_id,series_start,series_end);
          pp->put(cur);
          add(pp);
        }
        else
          APP_ABORT("PostProcessor::put  creation of "+cname+" postprocessor failed");
      }
      cur = cur->next;
    }
  }


  void PostProcessor::postprocess()
  {
    app_log()<<"  Postprocessing requested data"<<endl;
    for(int i=0;i<postprocessors.size();++i)
    {
      PostProcessorBase& pp = *postprocessors[i];
      app_log()<<"    postprocessing "+pp.type<<endl;
      pp.postprocess();
    }
    app_log()<<"  Postprocessing completed\n"<<endl;
  }

}
