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



#if !(defined LOOS_CHARMM_HPP)
#define LOOS_CHARMM_HPP

#include <ios>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdexcept>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  //! Class for reading a CHARMM coordinate file
  /**
   *
   * The code extracts the atom/residue/segment names and numbers, plus the weight,
   * which is assigned into the occupancy record.
   * The code supports both the small and large CHARMM crd formats.
   *
  */
  class CHARMM : public AtomicGroup {
  public:
    CHARMM() { }
    virtual ~CHARMM() {}

    explicit CHARMM(const std::string fname) {
      std::ifstream ifs(fname.c_str());
      if (!ifs) {
        throw(std::runtime_error("Cannot open CHARMM coordinate file " + 
                                 std::string(fname)));
      }
      read(ifs);
    }

    explicit CHARMM(std::istream &ifs) {
      read(ifs);
    }

    //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
    virtual CHARMM* clone(void) const;

    //! Creates a deep copy (see AtomicGroup::copy() for more info)
    CHARMM copy(void) const;

    void read(std::istream& is);  


  private:

    CHARMM(const AtomicGroup& grp) : AtomicGroup(grp) { }

    //void parseAtomRecord(const std::string s);  
  

  };


}

#endif
