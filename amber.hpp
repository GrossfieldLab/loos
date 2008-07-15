/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2008, Tod D. Romo, Alan Grossfield
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/




#if !defined(AMBER_HPP)
#define AMBER_HPP




#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tokenizer.hpp>

using namespace std;


#include <loos_defs.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>
#include <utils.hpp>

//! Class for reading in AMBER parmtop/coord files...
/** Notes:
 *o No check is made to make sure pointers have been read prior to other flags...
 */

class Amber : public AtomicGroup {
private:

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

public:

  Amber() : natoms(0), nres(0), nbonh(0), mbona(0) { }
  virtual ~Amber() { }

  explicit Amber(const string fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + fname));
    read(ifs);
  }

  explicit Amber(const char* fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(fname);
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + string(fname)));
    read(ifs);
  }

  explicit Amber(ifstream& ifs) : natoms(0), nres(0), nbonh(0), mbona(0) {
    read(ifs);
  }

  Amber copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    Amber p(grp);

    return(p);
  }

  void read(istream& is);



private:

  Amber(const AtomicGroup& grp) : AtomicGroup(grp), natoms(0), nres(0), nbonh(0), mbona(0) { }

  void verifyFormat(istream&, const string);
  void parseCharges(istream&);
  void parseMasses(istream&);
  void parseResidueLabels(istream&);
  void parseResiduePointers(istream&);
  void assignResidues(void);
  void parseBonds(istream&, const int);
  void parsePointers(istream&);
  void parseTitle(istream&);
  void parseAtomNames(istream&);


private:


  string _title;

  // These are internal and are used for parsing the parmtop info...
  uint natoms, nres, nbonh, mbona;

  vector<string> residue_labels;
  vector<uint> residue_pointers;


};





#endif


