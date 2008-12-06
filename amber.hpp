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
  explicit Amber(const std::string fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    std::ifstream ifs(fname.c_str());
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber parmtop file " + fname));
    read(ifs);
  }

  //! Read in a parmtop file and assign coordinates from the \a crds file.
  explicit Amber(const std::string parm, const std::string crds) : natoms(0), nres(0), nbonh(0), mbona(0) {
    std::ifstream ifs(parm.c_str());
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber parmtop file " + parm));
    read(ifs);
    readCoords(crds);
  }

  //! Read in a parmtop file
  explicit Amber(const char* fname) : natoms(0), nres(0), nbonh(0), mbona(0) {
    std::ifstream ifs(fname);
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber parmtop file " + std::string(fname)));
    read(ifs);
  }

  //! Read in a parmtop file and assign coordinates from the \a crds file.
  explicit Amber(const char* parm, const char* crds) : natoms(0), nres(0), nbonh(0), mbona(0) {
    std::ifstream ifs(parm);
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber parmtop file " + std::string(parm)));
    read(ifs);
    readCoords(crds);
  }


  explicit Amber(std::ifstream& ifs) : natoms(0), nres(0), nbonh(0), mbona(0) {
    read(ifs);
  }

  //! Clones an object for polymorphism...
  virtual Amber* clone(void) const {
    return(new Amber(*this));
  }

  //! Deep copy
  Amber copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    Amber p(grp);

    return(p);
  }

  //! Parse the parmtop file
  void read(std::istream& is);

  //! Read in a coord (or restart file)
  void readCoords(const char* fname) {
    std::ifstream ifs(fname);
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber crdinp file " + std::string(fname)));
    readCoords(ifs);
  }

  //! Read in a coord (or restart) file.
  void readCoords(const std::string& fname) {
    std::ifstream ifs(fname.c_str());
    if (!ifs)
      throw(std::runtime_error("Cannot open Amber crdinp file " + fname));
    readCoords(ifs);
  }


  //! Parses a coord/restart file from the given stream...
  /*!
   * This member function will auto-detect whether or not the file is
   * a coordinates file or a restart file and will handle either
   * appropriately.  Will also auto-detect box parameters and update
   * them, if present.
   */
  void readCoords(std::istream& is);



private:

  Amber(const AtomicGroup& grp) : AtomicGroup(grp), natoms(0), nres(0), nbonh(0), mbona(0) { }

  void verifyFormat(std::istream&, const std::string);
  void parseCharges(std::istream&);
  void parseMasses(std::istream&);
  void parseResidueLabels(std::istream&);
  void parseResiduePointers(std::istream&);
  void assignResidues(void);
  void parseBonds(std::istream&, const int);
  void parsePointers(std::istream&);
  void parseTitle(std::istream&);
  void parseAtomNames(std::istream&);


private:


  std::string _title;

  // These are internal and are used for parsing the parmtop info...
  uint natoms, nres, nbonh, mbona;

  std::vector<std::string> residue_labels;
  std::vector<uint> residue_pointers;


};





#endif


