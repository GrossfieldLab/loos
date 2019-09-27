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

typedef std::vector<uint> StateVector;
typedef std::function<bool(std::pair<StateVector, uint>, std::pair<StateVector, uint>)> Comparator;

class ChainState {
public:

ChainState() { }

ChainState(const uint segs, const uint bins) : _num_segs(segs),
                                               _num_bins(bins),
                                               counts(0)
                                               {
    _bin_width = 2.0 / _num_bins;
    compareStateProbs =
                [](std::pair<StateVector, uint> elem1,
                   std::pair<StateVector, uint> elem2)
                {
                return elem1.second > elem2.second;
                };
}

//! Compute the state of the chain and store it
void computeChainState(const AtomicGroup &group,
                       const GCoord &normal,
                       StateVector &segs);
void computeChainState(const AtomicGroup &group,
                       const GCoord &normal);


//! Return the probability of the state specified by segs
double getStateProb(const StateVector &segs);

//! Return all state probabilities
std::set<std::pair<StateVector, uint>, Comparator > getAllProbs() const;

//! compute Shannon entropy from the state distribution
double entropy();

//! Number of entries in the histogram
uint num_counts() const {
    return counts;
};

//! Number of states observed
uint num_states() const {
    return state_counts.size();
};

//! Number of possible states
uint total_states() const {
    float t = pow(static_cast<float>(_num_bins), _num_segs);
    return static_cast<uint>(t);
}

//! Relative entropy of the state probability distribution relative to reference
double relative_entropy(const std::map<StateVector, double> &ref);

private:
    uint _num_segs;
    uint _num_bins;
    double _bin_width;
    //boost::unordered_map<StateVector, uint> state_counts;
    std::map<StateVector, uint> state_counts;
    uint counts;
    Comparator compareStateProbs;

};

class RefChainDist {
public:
    std::map<StateVector, double> state_dist;

    RefChainDist() { }

    RefChainDist(const std::string &filename);

    RefChainDist(const ChainState &chainState);

    //! Read in reference distribution
    void readInput(const std::string &filename);

    //! Compute relative entropy of other wrt to this distribution
    double relative_entropy(RefChainDist &other);

private:
};

}

#endif
