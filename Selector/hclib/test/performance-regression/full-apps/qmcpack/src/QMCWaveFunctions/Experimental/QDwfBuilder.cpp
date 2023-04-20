#include "OhmmsData/libxmldefs.h"
#include "QMCWaveFunctions/QDwfBuilder.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{

bool QDwfBuilder::put(xmlNodePtr cur)
{
  typedef DiracDeterminant<SPOSet_t> DDet_t;
  typedef SlaterDeterminant<SPOSet_t> SDet_t;
  map<string,QDwf*> qdwfs;
  int nbasis=0;
  SDet_t* sdet = new SDet_t;
  cur = cur->xmlChildrenNode;         /// cur->name = DeterminantSet
  while(cur!=NULL)
  {
    string cname((const char*)(cur->name));
    if(cname =="SlaterDeterminant")
    {
      int first = 0;
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        string cname2((const char*)(tcur->name));
        if(cname2 == "Determinant")
        {
          SPOSet_t* swfs = new SPOSet_t;
          int norb
          =atoi((const char*)(xmlGetProp(tcur, (const xmlChar *)"orbitals"))) ;
          XMLReport("The number of orbitals for a Dirac Determinant " << norb)
          xmlNodePtr t = tcur->xmlChildrenNode;
          while(t != NULL)
          {
            string oname((const char*)(t->name));
            if(oname == "Orbital")
            {
              string
              orbname((const char*)(xmlGetProp(t, (const xmlChar *)"name")));
              map<string,QDwf*>::iterator it = qdwfs.find(orbname);
              if(it  == qdwfs.end())
              {
                XMLReport("Adding a new orbital " << orbname)
                QDwf *neworb = new QDwf();
                neworb->put(t,wfs_ref.VarList);
                swfs->add(neworb);
                qdwfs[oname] = neworb;
              }
              else
              {
                XMLReport("Using an existing " << orbname)
                swfs->add((*it).second);
              }
            }
            t = t->next;
          }
          XMLReport("The number of basis functions " << swfs->size())
          DDet_t* ddet = new DDet_t(*swfs,first);
          ddet->set(first,swfs->size());
          XMLReport("Adding a determinant to the SlaterDeterminant of size "
                    << swfs->size() << " begining with orbital " << first)
          first += swfs->size();
          sdet->add(ddet);
        } //dirac-determinant
        tcur = tcur->next;
      }
    }//slater-determinant
    cur = cur->next;
  }
  wfs_ref.add(sdet);
  return true;
}
}
