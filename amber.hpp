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




#if !defined(LOOS_AMBER_HPP)
#define LOOS_AMBER_HPP



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  //! Class for reading in AMBER parmtop/coord files...
  /*!
   * This class is largely geared towards reading parmtop files.  It
   * only parses a subset of the spec and follows more the format as
   * defined from example files and VMD than from the Amber website.
   *
   * Atomic numbers will be deduced from the masses.  No error is
   * generated if an atomic mass is unknown to LOOS.  In order to
   * verify that all atoms have an assigned mass, use the following,
   *\code
   * bool ok = amber.allHaveProperty(Atom::anumbit);
   *\endcode
   *
   */

  class Amber : public AtomicGroup {
  private:

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  public:

    Amber() : natoms(0), nres(0), nbonh(0), mbona(0),
              _lineno(0), _unget(false) { }
    virtual ~Amber() { }

    //! Read in a parmtop file
    explicit Amber(const std::string fname) : natoms(0), nres(0), nbonh(0), mbona(0),
                                              _lineno(0), _unget(false) {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
        throw(std::runtime_error("Cannot open Amber parmtop file " + fname));
      read(ifs);
    }

    //! Read in a parmtop file
    explicit Amber(const char* fname) : natoms(0), nres(0), nbonh(0), mbona(0),
                                        _lineno(0), _unget(false) {
      std::ifstream ifs(fname);
      if (!ifs)
        throw(std::runtime_error("Cannot open Amber parmtop file " + std::string(fname)));
      read(ifs);
    }

    explicit Amber(std::ifstream& ifs) : natoms(0), nres(0), nbonh(0), mbona(0),
                                         _lineno(0), _unget(false) {
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

    //! Return the title
    std::string title() const { return(_title); }

  private:

    Amber(const AtomicGroup& grp) : AtomicGroup(grp), natoms(0), nres(0), nbonh(0), mbona(0) { }

    void getNextLine(std::istream& is);

    void verifyFormat(std::istream&, const std::string&, const std::string&);
    void parseCharges(std::istream&);
    void parseMasses(std::istream&);
    void parseResidueLabels(std::istream&);
    void parseResiduePointers(std::istream&);
    void assignResidues(void);
    void parseBonds(std::istream&, const uint);
    void parsePointers(std::istream&);
    void parseTitle(std::istream&);
    void parseAtomNames(std::istream&);
    void parseAmoebaRegularBondNumList(std::istream&);
    void parseAmoebaRegularBondList(std::istream&, const uint);


    template<typename T>
    std::vector<T> readBlock(std::istream& is, const int field_width) {
      std::vector<T> data;
      while (true) {
        getNextLine(is);
        if (is.eof())
          break;
        if (is.fail())
          throw(FileParseError("Error while reading block of data from amber file", _lineno));
        if (_current_line[0] == '%') {
          _unget = true;
          break;
        }
        std::istringstream iss(_current_line);
        T d;
        while (iss >> std::setw(field_width) >> d)
          data.push_back(d);
      }

      return(data);
    }

  private:


    std::string _title;

    // These are internal and are used for parsing the parmtop info...
    uint natoms, nres, nbonh, mbona, _amoeba_regular_bond_num_list;

    std::vector<std::string> residue_labels;
    std::vector<uint> residue_pointers;


    std::string _current_line;
    uint _lineno;
    bool _unget;
  };


}



#endif


