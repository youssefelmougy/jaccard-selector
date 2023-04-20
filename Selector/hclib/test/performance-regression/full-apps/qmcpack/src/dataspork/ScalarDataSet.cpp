#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include "ScalarDataSet.h"
#include "ScalarObservable.h"
#include <boost/tokenizer.hpp>
using namespace std;

ScalarDataSet::ScalarDataSet(const string& fname, int first, int last):Grouped(false)
{
  read(fname,first,last);
}

ScalarDataSet::~ScalarDataSet()
{
  RawDataSetType::iterator oit(Data.begin()), oit_end(Data.end());
  while(oit != oit_end)
  {
    delete (*oit).second;
    ++oit;
  }
}

ScalarDataSet*
ScalarDataSet::clone(int nc)
{
  if(nc<0)
    nc=NumRows;
  ScalarDataSet* c= new ScalarDataSet;
  RawDataSetType::const_iterator it(Data.begin()), it_end(Data.end());
  while(it != it_end)
  {
    ScalarObservable* s=new ScalarObservable;
    s->resize(nc,0.0);
    c->Data[(*it).first]=s;
    ++it;
  }
  return c;
}

void ScalarDataSet::read(const string& fname, int first, int last)
{
  ifstream fin(fname.c_str());
  if(!fin.is_open())
  {
    cerr << " Cannot open " << fname << endl;
    return;
  }
  string aline;
  vector<std::string> columnList;
  bool find_data=true;
  int nskip=0;
  while(find_data && !fin.eof())
  {
    streampos cur=fin.tellg();
    if(std::getline(fin,aline))
    {
      typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
      boost::char_separator<char> sep(" ");
      tokenizer tok(aline, sep);
      tokenizer::iterator it=tok.begin();
      if(isdigit((*it)[0]))
      {
        find_data=false;
        fin.seekg(cur);
      }
      else
      {
        nskip++;
        columnList.clear();
        if((*it).find("#") == 0)
        {
          ++it;
        }
        for(; it!=tok.end(); it++)
        {
          columnList.push_back(*it);
        }
      }
    }
  }
  //temporary vector for the columns so that we can use index
  vector<ScalarObservable*> dtemp;
  vector<std::string>::iterator cit(columnList.begin()), cit_end(columnList.end());
  while(cit != cit_end)
  {
    map<string,ScalarObservable*>::iterator dit(Data.find(*cit));
    if(dit == Data.end())
      //add to newSet
    {
      ScalarObservable* newO=new ScalarObservable;
      dtemp.push_back(newO);
      Data[*cit]=newO;
    }
    else
    {
      //this should be optional for append
      (*dit).second->Self.clear();
      dtemp.push_back((*dit).second);
    }
    ++cit;
  }
  int nc=columnList.size();
  vector<double> arow(nc);
  double firstColumn=-1;
  if(last<0)
    last=numeric_limits<int>::max();
  int nr=0;
  while(nr<first)
  {
    std::getline(fin,aline);
    nr++;
  }
  //treat the first column as string
  string dummy,prev("invalid");
  while(nr < last && !fin.eof())
  {
    fin >> dummy;
    for(int ic=1; ic<nc; ic++)
      fin >> arow[ic];
    if(dummy != prev)
    {
      dtemp[0]->push_back(nr);
      for(int ic=1; ic<nc; ic++)
        dtemp[ic]->push_back(arow[ic]);
      nr++;
    }
    prev=dummy;
  }
  FirstRow=first;
  LastRow=nr;
  NumRows=dtemp[0]->size();
}

void ScalarDataSet::accumulate(const RawDataSetType& raw)
{
  RawDataSetType::const_iterator sit(raw.begin()), sit_end(raw.end());
  RawDataSetType::iterator tit(Data.begin());
  while(sit != sit_end)
    //loop over the ScalarObservables
  {
    ScalarObservable::iterator in((*sit).second->begin());
    ScalarObservable::iterator out((*tit).second->begin());
    ScalarObservable::iterator out_end((*tit).second->end());
    while(out!=out_end)
    {
      *out++ += *in++;
    }
    ++sit;
    ++tit;
  }
}

void
ScalarDataSet::group(const map<string,ScalarDataSet*>& dlist)
{
  Grouped=true;
  FirstRow=0;
  LastRow=NumRows=(*Data.begin()).second->size();
  map<string,ScalarDataSet*>::const_iterator first(dlist.begin()), last(dlist.end());
  double dn=0;
  while(first != last)
  {
    accumulate((*first).second->Data);
    dn+=1.0;
    ++first;
  }
  //divide by the sample numbers
  dn=1.0/dn;
  RawDataSetType::iterator tit(Data.begin()), tit_end(Data.end());
  while(tit != tit_end)
    //loop over the ScalarObservables
  {
    ScalarObservable::iterator out((*tit).second->begin());
    ScalarObservable::iterator out_end((*tit).second->end());
    while(out!=out_end)
    {
      *out++ *= dn;
    }
    ++tit;
  }
}

void
ScalarDataSet::writeCollect(ostream& os, const string& cname)
{
  map<string,ScalarObservable*>::iterator dit(Data.find(cname));
  if(dit != Data.end())
  {
    os << "    <collect name=\""<<(*dit).first << "\""
       << " value=\"" <<(*dit).second->get_average()
       << "\"/>\n";
    cout << (*dit).first << " " << (*dit).second->get_average() << endl;
  }
}

void
ScalarDataSet::writeObservable(ostream& os, const string& oname)
{
  map<string,ScalarObservable*>::iterator dit(Data.find(oname));
  if(dit != Data.end())
  {
    cout.setf(ios::left,ios::adjustfield);
    cout << setw(16) << (*dit).first;
    cout.setf(ios::right,ios::adjustfield);
    os << "    <observable name=\""<<(*dit).first << "\">\n";
    (*dit).second->get_stat();
    (*dit).second->write(os);
    os << "    </observable>\n";
  }
}
