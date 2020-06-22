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




#if !defined(LOOS_CIF_HPP)
#define LOOS_CIF_HPP

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <map>

// openbabel includes
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/residue.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/generic.h>
#include <openbabel/obiter.h>
#include <openbabel/oberror.h>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <cryst.hpp>
#include <exceptions.hpp>


namespace loos {

//! PDBX/mmcif reading/writing class
/** This class reads and writes PDBX/mmcif format by wrapping openbabel. It will
 *  use the cell records in the file (if present) to set the periodic box
 */
class MMCIF : public AtomicGroup {
public:
    MMCIF() : _fname("<not set>"), _has_cryst(false)   { }
    virtual ~MMCIF() {}

    //! Create an mmcif given a filename
    explicit MMCIF(const std::string& fname)
        : _fname(fname), _has_cryst(false)
    {
        std::ifstream ifs(fname.c_str());
        if (!ifs)
        {
            throw(FileOpenError(fname));
        }
        read(ifs);
    }

    //! Create an mmcif given an ifstream
    explicit MMCIF(std::istream& ifs)
        : _fname("stream"), _has_cryst(false)
    {
        read(ifs);
    }


   static pAtomicGroup create(const std::string& fname) {
       return(pAtomicGroup(new MMCIF(fname)));
   }

   //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
   virtual MMCIF* clone(void) const;

   //! Creates a deep copy (see AtomicGroup::copy() for more info)
   MMCIF copy(void) const;

   //! Class method for creating an MMCIF from an AtomicGroup
   /** There should probably be some internal checks to make sure we
    *  have enough info to actually write out an MMCIF, but currently
    *  there are no such checks...
    */
   static MMCIF fromAtomicGroup(const AtomicGroup&);

   //! Create an MMCIF from an AtomicGroup (i.e. upcast)
   MMCIF(const AtomicGroup& grp) : AtomicGroup(grp) { }

    //! Read in a mmcif file from an istream
    void read(std::istream& ifs);

    const UnitCell& unitCell(void);
    void unitCell(const UnitCell& c);

    //! Create on OpenBabel object, primarily on the way to writing out
    /**
     * This should probably be private, but for testing purposes it's easier
     * to have it here.
     */
    OpenBabel::OBMol * toOpenBabel(void) const;


#if !defined(SWIG)
    //! Output as an MMCIF
        friend std::ostream& operator<<(std::ostream&, const MMCIF&);
#endif

private:
    std::string _fname;
    bool _has_cryst;
    std::map<uint, pAtom> _atomid_to_patom;
    UnitCell cell;

    pAtom findAtom(const uint id);
};



}

#endif
