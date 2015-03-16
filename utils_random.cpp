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


#include <boost/random.hpp>
#include <loos/utils_random.hpp>

namespace loos {
  base_generator_type& rng_singleton(void) {
    static base_generator_type rng;

    return(rng);
  }


  // Seeding based on the block is not the best method, but probably
  // sufficient for our purposes...
  uint randomSeedRNG(void) {
    base_generator_type& rng = rng_singleton();

    uint seedval = static_cast<uint>(time(0));

    rng.seed(seedval);
    return(seedval);
  }
};
