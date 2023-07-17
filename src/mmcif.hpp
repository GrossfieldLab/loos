/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2023, Tod D. Romo, Alan Grossfield
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

#if !defined(LOOS_MMCIF_HPP)
#define LOOS_MMCIF_HPP

#include <loos.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>    // cif::Document -> Structure
#include <gemmi/gz.hpp> 

namespace loos {

    //! Class to read pdbx/mmcif files

    class MMCIF : public AtomicGroup {
    public:
        MMCIF() {}
        virtual ~MMCIF() {}

        explicit MMCIF(const std::string& filename) : _filename(filename) {
            read(filename);
        }

        static pAtomicGroup create(const std::string& filename) {
            return pAtomicGroup(new MMCIF(filename));
        }
    
        //! Clones an object for polymorphism (see AtomicGroup::clone() for more info)
        virtual PDB* clone(void) const;

        void read(const std::string& filename);

    private:
        std::string _filename;
    };

}
#endif