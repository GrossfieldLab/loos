/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019, Alan Grossfield
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

#if !defined(LOOS_CHAINSTATE_HPP)
#define LOOS_CHAINSTATE_HPP

#include "loos.hpp"

namespace loos {

class ChainState {
public:

ChainState() {}

ChainState(const uint segs, const uint bins) : _num_segs(segs),
                                               _num_bins(bins)
                                               {
    _bin_width = 2.0 / _num_bins;
}

void computeChainState(const AtomicGroup &group,
                       const GCoord &normal,
                       std::vector<uint> &segs);
double getStateProb(const std::vector<uint> &segs);

private:
    uint _num_segs;
    uint _num_bins;
    double _bin_width;
    std::map<std::vector<uint>, uint> _state_counts;
    uint _counts;

};

}

#endif
