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



#if !(defined PSF_HPP)
#define PSF_HPP

#include <ios>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdexcept>
#include <ctype.h>

#include <loos_defs.hpp>
#include <Atom.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  //! Class for reading a subset of the PSF format
  /**
   * Notes:
   *
   * PSF files don't have coordinates, so the coordinates property won't be
   * set.  As a result, you'll get an UnsetProperty exception if you try to 
   * access the coordinates.  The idea is you use a PSF with something else
   * which will supply the coordinates (usually some kind of Trajectory).
   *
   * This code will read both NAMD and CHARMM PSF files, since it doesn't use
   * the atom type information anyway.
   *
   * The code extracts the atom/residue/segment names and numbers, plus mass, 
   * partial charge, and connectivity.  Higher order connectivity (angles, 
   * dihedrals, etc) are ignored.
   *
  */
  class PSF : public AtomicGroup {
  public:
    PSF() { }
    virtual ~PSF() {}

    explicit PSF(const std::string fname) {
      std::ifstream ifs(fname.c_str());
      if (!ifs) {
        throw(std::runtime_error("Cannot open PSF file " + std::string(fname)));
      }
      read(ifs);
    }

    explicit PSF(std::ifstream &ifs) {
      read(ifs);
    }

    //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
    virtual PSF* clone(void) const {
      return(new PSF(*this));
    }

    //! Creates a deep copy (see AtomicGroup::copy() for more info)
    PSF copy(void) const {
      AtomicGroup grp = this->AtomicGroup::copy();
      PSF p(grp);

      // Add PSF specific member data copies here...
      return(p);
    }

    void read(std::istream& is);  


  private:

    PSF(const AtomicGroup& grp) : AtomicGroup(grp) { }


    void parseAtomRecord(const std::string s);  
  
    //! Use the mass to deduce the atomic number of the atom
    int deduceAtomicNumber(pAtom pa);

  };


}

#endif
