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


#if !defined(LOOS_UTILS_RANDOM_HPP)
#define LOOS_UTILS_RANDOM_HPP

#include <boost/random.hpp>
#include <loos_defs.hpp>

namespace loos {

  typedef boost::mt19937 base_generator_type;

  //! Suite-wide random number generator singleton
  /**
   * LOOS makes no assumptions about how the random number generator
   * gets seeded.  It is up to the tool-writer to seed it with a known
   * value,
\code
rng_singleton().seed(seed_value);
\endcode
   * or call randomSeedRNG() to randomly seed the random number
   * generator...
   */
  base_generator_type& rng_singleton(void);

  //! Randomly seeds the RNG
  /**Currently uses time(3) to seed the RNG obtained from the singleton...
   * Returns the seed used.
   */
  uint randomSeedRNG(void);

};


#endif
