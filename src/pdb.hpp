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
#include <exceptions.hpp>
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
    public:
        PDB() : _max_index(0), _show_charge(false), _auto_ter(true), _has_cryst(false),
                strictness_policy(false), _missing_q(false), _missing_b(false),
                _missing_segid(false), _fname("<not set>") { }
        virtual ~PDB() {}
      
        //! Read in PDB from a filename
        explicit PDB(const std::string& fname) 
            : _max_index(0), _show_charge(false), _auto_ter(true),
              _has_cryst(false), strictness_policy(false),
              _missing_q(false), _missing_b(false), _missing_segid(false),
              _fname(fname)
        {
            std::ifstream ifs(fname.c_str());
            if (!ifs)
                throw(FileOpenError(fname));
            read(ifs);
        }
      
        //! Read in a PDB from an ifstream
        explicit PDB(std::istream& ifs) 
            : _max_index(0), _show_charge(false), _auto_ter(true),
              _has_cryst(false), strictness_policy(false),
              _missing_q(false), _missing_b(false), _missing_segid(false),              
              _fname("stream")
        {
            read(ifs);
        }
      
        static pAtomicGroup create(const std::string& fname) {
            return(pAtomicGroup(new PDB(fname)));
        }

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

        //! Special handling for charges since the PDB form is daft
        void showCharge(bool b);


        bool showCharge(void) const;
        bool strict(void) const;
        bool autoTerminate(void) const;

        //! Determins how strict the input parser is to the '96 PDB standard.
        void strict(const bool b);

        //! Automatically insert TER record at end of output.
        /*! Controls whether or not a "TER" record is automatically added
          when printing out the group/PDB
        */
        void autoTerminate(bool b);

        //! Accessor for the remarks object...
        Remarks& remarks(void);
        //! Accessor for the remarks object...
        void remarks(const Remarks&);

        const UnitCell& unitCell(void);
        void unitCell(const UnitCell&);

        //! Output as a PDB
        friend std::ostream& operator<<(std::ostream&, const PDB&);

        std::string asString() const;

        //! Read in PDB from an ifstream
        void read(std::istream& is);

    private:
        class ComparePatoms {
            bool operator()(const pAtom& a, const pAtom& b) { return(a->id() < b->id()); }
        };


        //! Create a PDB from an AtomicGroup (i.e. upcast)
        PDB(const AtomicGroup& grp) : AtomicGroup(grp), _show_charge(false), _auto_ter(true), _has_cryst(false) { }

        bool emptyString(const std::string&);

    bool isMissingOccupancies() const { return(_missing_q); }
    bool isMissingBFactors() const { return(_missing_b); }
    bool isMissingSegids() const { return(_missing_segid); }
    bool isMissingFields() const { return(_missing_q || _missing_b || _missing_segid); }


        // These will modify the PDB upon a successful parse...
        void parseRemark(const std::string&);
        void parseAtomRecord(const std::string&);
        void parseConectRecord(const std::string&);
        void parseCryst1Record(const std::string&);

        // Convert an Atom to a string representation in PDB format...
        std::string atomAsString(const pAtom p) const;

#if !defined(SWIG)
        friend std::ostream& FormatConectRecords(std::ostream&, const PDB&);
#endif

        pAtom findAtom(const int i);
        void uniqueBonds();

    private:
        uint _max_index;
        bool _show_charge;
        bool _auto_ter;
        bool _has_cryst;
        bool strictness_policy;
        bool _missing_q, _missing_b, _missing_segid;
        std::string _fname;
        Remarks _remarks;
        UnitCell cell;
        std::map<int, pAtom> _atomid_to_patom;
    };

}

#endif

