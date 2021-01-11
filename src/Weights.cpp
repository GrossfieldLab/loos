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

#include "Weights.hpp"
#include <sstream>

namespace loos {
//! Weights class to handle reweighting values computed from a trajectory

//! Normalize the weights so they sum to 1
inline void Weights::normalize() {
  double sum = 0.0;
  for (uint i = 0; i < _weights.size(); ++i) {
    sum += _weights[i];
  }
  // TODO : Really should check for underflow to prevent div by 0
  for (uint i = 0; i < _weights.size(); ++i) {
    _weights[i] /= sum;
  }
}

//! Keep track of total weight used
inline void Weights::accumulate() {
  _total += _weights.at(_traj->currentFrame());
  _totalTraj += _weights.at(_traj->currentFrame());
}

inline void Weights::accumulate(const uint index) {
  _total += _weights.at(index);
  _totalTraj += _weights.at(index);
}

//! Return the totalWeight, as tracked using accumulate
inline const double Weights::totalWeight() { return _total; }

//! Return the weight of the current trajectory
inline const double Weights::trajWeight() { return _totalTraj; }

//! Return the weight for the current frame of the trajectory
inline const double Weights::get() {
  current_frame = _traj->currentFrame();
  return _weights.at(current_frame);
}

//! Return the weight for frame index of the trajectory
inline const double Weights::get(const uint index) { return _weights.at(index); }

//! calling nomenclature wraps get
inline const double Weights::operator()() { return get(); }

inline const double Weights::operator()(const uint index) { return get(index); }

//! Bind a new weight to the current frame
inline void Weights::set(double newWeight) { _weights.at(current_frame) = newWeight; }

//! Bind a new weight to a particular frame
inline void Weights::set(double newWeight, const uint index) {
  _weights.at(index) = newWeight;
}

//! binding nomenclature wraps set
inline void Weights::operator()(double newWeight) { set(newWeight); }
inline void Weights::operator()(double newWeight, const uint index) {
  set(newWeight, index);
}

//! set all weights from passed vector
inline void Weights::operator()(std::vector<double> &newWeights) {
  if (_weights.size() != newWeights.size())
    throw(LOOSError(std::string(
      "Number of weights in class is " + std::to_string(_weights.size())
      + " number inserted is " + std::to_string(newWeights.size())
      + " these must match."
    )));
  // note that operator= for stl::vectors does in fact copy their contents.
  // made this a move operation to speed up transfer; obviously this leaves argument empty.
  _weights = std::move(newWeights);
}

//! bind the provided pTraj to this instance of Weights 
void Weights::addTraj(pTraj &traj){
  _traj = traj;
  _totalTraj = 0.0;
}

//! Return the number of weights
inline uint Weights::size() { return _num_weights; }

//! Return the vector of weights
inline std::vector<double> Weights::weights() { return _weights; }
} // namespace loos
