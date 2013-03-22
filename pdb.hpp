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




#if !defined(LOOS_PDB_HPP)
#define LOOS_PDB_HPP




#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <map>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <pdb_remarks.hpp>
#include <cryst.hpp>
#include <utils.hpp>
#include <utils_structural.hpp>



namespace loos {

  //! PDB reading/writing class
  /**This class models a basic PDB file format.  Special handling is
   *included for dealing with periodic boundary conditions.  If there is
   *a special REMARK header, then the box size is extracted from this
   *and stored in the parent AtomicGroup.  If not, but a CRYST1 record
   *is present, then the a, b, c parameters are used to set the periodic
   *box.
   *
   * This class can handle some minor variations in the PDB format
   * (since there are so many different standards out there).  This is
   * determined by the strictness_policy setting.  Be default, a weak
   * strictness is used.  This will allow variations like frame-shifts
   * in the resid field to be accepted.  If you would rather have the
   * strict formatting honored, you'll need to set each PDB object to
   * strict:
   \code
   PDB pdb;
   pdb.strict(true);
   \endcode
  */
  class PDB : public AtomicGroup {
    enum BigNumberFormatType { HYBRID36, HEX, NONE };

  public:
    class BadConnectivity : public std::runtime_error {
    public:
      explicit BadConnectivity(const std::string& msg) : runtime_error(msg) { };
    };

  public:
    PDB() : _show_charge(false), _auto_ter(true), _output_format(NONE), _input_format(HYBRID36),
            _has_cryst(false), strictness_policy(false),
            _format_warned(false) { }
    virtual ~PDB() {}

    //! Read in PDB from a filename
    explicit PDB(const std::string fname) : _show_charge(false), _auto_ter(true),
                                            _output_format(NONE), _input_format(HYBRID36),
                                            _has_cryst(false), strictness_policy(false),
                                            _format_warned(false)
    {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
        throw(std::runtime_error("Cannot open PDB file " + fname));
      read(ifs);
    }

    //! Read in a PDB from a filename
    explicit PDB(const char* fname) : _show_charge(false), _auto_ter(true),
                                      _output_format(NONE), _input_format(HYBRID36),
                                      _has_cryst(false), strictness_policy(false),
                                      _format_warned(false)
    {
      std::ifstream ifs(fname);
      if (!ifs)
        throw(std::runtime_error("Cannot open PDB file " + std::string(fname)));
      read(ifs);
    }

    //! Read in a PDB from an ifstream
    explicit PDB(std::ifstream& ifs) : _show_charge(false), _auto_ter(true),
                                       _output_format(NONE), _input_format(HYBRID36),
                                       _has_cryst(false), strictness_policy(false),
                                       _format_warned(false)
    { read(ifs); }


    //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
    virtual PDB* clone(void) const;

    //! Creates a deep copy (see AtomicGroup::copy() for more info)
    PDB copy(void) const;

    //! Class method for creating a PDB from an AtomicGroup
    /** There should probably be some internal checks to make sure we
     *  have enough info to actually write out a PDB, but currently
     *  there are no such checks...
     */
    static PDB fromAtomicGroup(const AtomicGroup&);

    BigNumberFormatType inputFormat() const { return(_input_format); }
    void outputFormat(const BigNumberFormatType t) { _output_format = t; }
    BigNumberFormatType outputFormat() const { return(_output_format); }

    bool showCharge(void) const;
    //! Special handling for charges since the PDB form is daft
    void showCharge(bool b = true);

    bool strict(void) const;
    //! Determins how strict the input parser is to the '96 PDB standard.
    void strict(const bool b);

    bool autoTerminate(void) const;

    //! Automatically insert TER record at end of output.
    /*! Controls whether or not a "TER" record is automatically added
      when printing out the group/PDB
    */
    void autoTerminate(bool b = true);

    //! Accessor for the remarks object...
    Remarks& remarks(void);
    //! Accessor for the remarks object...
    void remarks(const Remarks&);

    const UnitCell& unitCell(void);
    void unitCell(const UnitCell&);

    //! Output as a PDB
    friend std::ostream& operator<<(std::ostream&, const PDB&);

    //! Set how to handle large numbers...
    void setOutputNumberFormatAsHybrid36() { _output_format = HYBRID36; }
    void setOutputNumberFormatAsHex() { _output_format = HEX; }


    //! Read in PDB from an ifstream
    void read(std::istream& is);

  private:
    class ComparePatoms {
      bool operator()(const pAtom& a, const pAtom& b) { return(a->id() < b->id()); }
    };


    //! Create a PDB from an AtomicGroup (i.e. upcast)
    PDB(const AtomicGroup& grp) : AtomicGroup(grp), _show_charge(false), _auto_ter(true), _output_format(HEX), _has_cryst(false) { }

    bool emptyString(const std::string&);

    // These will modify the PDB upon a successful parse...
    void parseRemark(const std::string&);
    void parseAtomRecord(const std::string&);
    void parseConectRecord(const std::string&);
    void parseCryst1Record(const std::string&);

    // Convert an Atom to a string representation in PDB format...
    std::string atomAsString(const pAtom p) const;

    std::string numberToString(const int i, const int width) const;

    friend std::ostream& FormatConectRecords(std::ostream&, const PDB&);

    pAtom findAtom(const int i);
    void uniqueBonds();
    int stringToNumber(const std::string& s, const uint start, const uint width);

  private:
    bool _show_charge;
    bool _auto_ter;
    BigNumberFormatType _output_format;
    BigNumberFormatType _input_format;
    bool _has_cryst;
    bool strictness_policy;
    bool _format_warned;
    Remarks _remarks;
    UnitCell cell;
    std::map<int, pAtom> _atomid_to_patom;
  };

}

#endif

