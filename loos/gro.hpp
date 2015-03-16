/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


#if !defined(LOOS_GRO_HPP)
#define LOOS_GRO_HPP


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <loos/loos_defs.hpp>
#include <loos/AtomicGroup.hpp>

namespace loos {

  //! Implements a GROMACS model file (.gro)
  class Gromacs : public AtomicGroup {
  public:

    Gromacs() { }
    
    explicit Gromacs(const std::string& fname) : _filename(fname), _max_index(0) {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
        throw(FileOpenError(fname));
      read(ifs);
    }

    explicit Gromacs(std::istream& ifs) : _filename("stream"), _max_index(0) { read(ifs); }

    static pAtomicGroup create(const std::string& fname) {
      return(pAtomicGroup(new Gromacs(fname)));
    }


    std::string title(void) const { return(title_); }

  private:
    std::string _filename;
    std::string title_;
    uint _max_index;

  private:
    
    void read(std::istream& ifs);

  };



}


#endif
