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



#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <LineReader.hpp>


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


    struct FormatSpec {
      FormatSpec() : repeat(1), type('?'), width(0), precision(0) { }
      int repeat;
      char type;
      int width;
      int precision;
    };

    struct AmberLineReader : public LineReader {
      AmberLineReader() : LineReader() { }
      AmberLineReader(std::istream& is) : LineReader(is) { }
      AmberLineReader(std::istream& is, std::string& s) : LineReader(is, s) { }

      virtual void stripComment(std::string& s) const {
        if (s.compare(0, 8, "%COMMENT") == 0)
          s = "";
      }
    };

  public:

    Amber() : natoms(0), nres(0), nbonh(0), mbona(0)  { }
    virtual ~Amber() { }

    //! Read in a parmtop file
    explicit Amber(const std::string& fname)
      : natoms(0), nres(0), nbonh(0), mbona(0) {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
	throw(FileOpenError(fname));
      reader.stream(ifs);
      reader.name(fname);
      read(ifs);
    }

    explicit Amber(std::istream& ifs)
      : natoms(0), nres(0), nbonh(0), mbona(0), reader(ifs) {
      read(ifs);
    }

    static pAtomicGroup create(const std::string& fname) {
      return(pAtomicGroup(new Amber(fname)));
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
    void read(std::istream& ifs);

    //! Return the title
    std::string title() const { return(_title); }

  private:

    Amber(const AtomicGroup& grp) : AtomicGroup(grp), natoms(0), nres(0), nbonh(0), mbona(0) { }

    FormatSpec parseFormat(const std::string& expected_types, const std::string& where);

    void parseCharges();
    void parseMasses();
    void parseResidueLabels();
    void parseResiduePointers();
    void assignResidues(void);
    void parseBonds(const uint);
    void parsePointers();
    void parseTitle();
    void parseAtomNames();
    void parseBoxDimensions();
    void parseAmoebaRegularBondNumList();
    void parseAmoebaRegularBondList(const uint);


    // Reads in a "block" of data.  Reading terminates on the first
    // line that begins with a '%'.

    template<typename T>
    std::vector<T> readBlock(const int field_width) {
      std::vector<T> data;
      while (reader.getNext()) {
        std::string line = reader.line();
        if (line[0] == '%') {
          reader.push_back(line);
          break;
        }
        std::istringstream iss(line);
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

    AmberLineReader reader;

  };


}



#endif
