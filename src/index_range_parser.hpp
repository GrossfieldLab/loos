/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2016, Tod D. Romo, Alan Grossfield
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



#if !defined(LOOS_INDEX_RANGE_PARSER)
#define LOOS_INDEX_RANGE_PARSER


#include <vector>
#include <string>

#include <loos_defs.hpp>

namespace loos {

  //! Use boost::Spirit to parse a new-style range list
  /**
   * Similar to parseRangeList(), but uses boost::spirit can requires a
   * maximum index size.  This allows you to specify ranges without knowing
   * the endpoint, e.g. 10:2:, means start at 10, skipping 2, until the value
   * of maxsize.  This is useful for trajectories, otherwise you would have to
   * know how big they are before using a range to specify a skip and stride.
   */
  std::vector<uint> parseIndexRange(const std::string& input, const uint maxsize);

}


#endif 
