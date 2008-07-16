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
/*!
 * This class is largely geared towards reading parmtop files.  It
 * only parses a subset of the spec and follows more the format as
 * defined from example files and VMD than from the Amber website.
 */

class Amber : public AtomicGroup {
private:

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

public:

  Amber() : natoms(0), nres(0), nbonh(0), mbona(0) { }
  virtual ~Amber() { }

  //! Read in a parmtop file
  explicit Amber(const string fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + fname));
    read(ifs);
  }

  //! Read in a parmtop file and assign coordinates from the \a crds file.
  explicit Amber(const string parm, const string crds) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(parm.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + parm));
    read(ifs);
    readCoords(crds);
  }

  //! Read in a parmtop file
  explicit Amber(const char* fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(fname);
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + string(fname)));
    read(ifs);
  }

  //! Read in a parmtop file and assign coordinates from the \a crds file.
  explicit Amber(const char* parm, const char* crds) : natoms(0), nres(0), nbonh(0), mbona(0) {
    ifstream ifs(parm);
    if (!ifs)
      throw(runtime_error("Cannot open Amber parmtop file " + string(parm)));
    read(ifs);
    readCoords(crds);
  }


  explicit Amber(ifstream& ifs) : natoms(0), nres(0), nbonh(0), mbona(0) {
    read(ifs);
  }

  //! Deep copy
  Amber copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    Amber p(grp);

    return(p);
  }

  //! Parse the parmtop file
  void read(istream& is);

  //! Read in a coord (or restart file)
  void readCoords(const char* fname) {
    ifstream ifs(fname);
    if (!ifs)
      throw(runtime_error("Cannot open Amber crdinp file " + string(fname)));
    readCoords(ifs);
  }

  //! Read in a coord (or restart) file.
  void readCoords(const string fname) {
    ifstream ifs(fname.c_str());
    if (!ifs)
      throw(runtime_error("Cannot open Amber crdinp file " + fname));
    readCoords(ifs);
  }


  //! Parses a coord/restart file from the given stream...
  /*!
   * This member function will auto-detect whether or not the file is
   * a coordinates file or a restart file and will handle either
   * appropriately.  Will also auto-detect box parameters and update
   * them, if present.
   */
  void readCoords(istream& is);



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


