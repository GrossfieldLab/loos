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

// openbabel includes
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/residue.h>
#include <openbabel/atom.h>
#include <openbabel/obiter.h>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <exceptions.hpp>


namespace loos {

//! PDBX/mmcif reading/writing class
/** This class reads and writes PDBX/mmcif format by wrapping openbabel. It will
 *  use the cell records in the file (if present) to set the periodic box
 */
class MMCIF : public AtomicGroup {
public:
    MMCIF() : _fname("<not set>")   { }
    virtual ~MMCIF() {}

    //! Create an mmcif given a filename
    explicit MMCIF(const std::string& fname)
        : _fname(fname)
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
        : _fname("stream")
    {
        read(ifs);
    }

    //! Read in a mmcif file from an istream
    void read(std::istream& ifs);

private:
    std::string _fname;

};



}

#endif
