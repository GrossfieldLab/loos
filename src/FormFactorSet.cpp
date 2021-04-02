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

#include <loos.hpp>
#include <FormFactorSet.hpp>

namespace loos {

    void FormFactorSet::setup() {
        // Add the nuclei we know to the map
        // H
        _map.insert(std::pair<uint,FormFactor>(1, FormFactor(1)));
        // C
        _map.insert(std::pair<uint,FormFactor>(6, FormFactor(6)));
        // N
        _map.insert(std::pair<uint,FormFactor>(7, FormFactor(7)));
        // O
        _map.insert(std::pair<uint,FormFactor>(8, FormFactor(8)));
        // P
        _map.insert(std::pair<uint,FormFactor>(15, FormFactor(15)));
        // S
        _map.insert(std::pair<uint,FormFactor>(16, FormFactor(16)));
    }

}
