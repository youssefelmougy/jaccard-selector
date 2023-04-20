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
#ifndef OHMMS_RANDOMNUMBERCONTROL_H__
#define OHMMS_RANDOMNUMBERCONTROL_H__
#include "OhmmsData/OhmmsElementBase.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/PrimeNumberSet.h"

class Communicate;

namespace qmcplusplus
{

/**class RandomNumberControl
 *\brief Encapsulate data to initialize and save the status of the random number generator
 *
 * Default:  myName = "random"
 * 2007-12-01
 *   Use PrimeNumbers to generate random seeds.
 */
class RandomNumberControl : public OhmmsElementBase
{

public:

  typedef RandomGenerator_t::uint_type uint_type;
  static PrimeNumberSet<uint_type> PrimeNumbers;
  //children random number generator
  static std::vector<RandomGenerator_t*>  Children;

  /// constructors and destructors
  RandomNumberControl(const char* aname="random");

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();
  static void test();

  static void make_seeds();
  static void make_children();

  xmlNodePtr initialize(xmlXPathContextPtr);

  /** read random state from a hdf file
   * @param fname file name
   * @param comm communicator so that everyone reads its own data
   */
  static void read(const string& fname, Communicate* comm);
  /** write random state to a hdf file
   * @param fname file name
   * @param comm communicator so that everyone writes its own data
   */
  static void write(const string& fname, Communicate* comm);

private:

  bool NeverBeenInitialized;
  xmlNodePtr myCur;
  static uint_type Offset;
};
}

#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 6100 $   $Date: 2013-12-04 11:42:14 -0500 (Wed, 04 Dec 2013) $
 * $Id: RandomNumberControl.h 6100 2013-12-04 16:42:14Z jnkim $
 ***************************************************************************/
