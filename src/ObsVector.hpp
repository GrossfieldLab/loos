/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2019, Louis G. Smith
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

// ObsVector.hpp
// A template for keeping the observables calculated at each frame in ram
#if !defined(LOOS_FRAME_VECTOR_HPP)
#define LOOS_FRAME_VECTOR_HPP
#include <loos_defs.hpp>
#include <vector>

namespace loos {
// template class that can contain the data from a trajectory as a vector
template <typename Measurement>
class ObsVector {
  pTraj _traj;
  std::vector<Measurement> _obsVector;

 public:
  ObsVector(pTraj& traj) : _traj(traj), { std::vector<Measurement> _obsVector; }
  ~ObsVector();
  // ith element corresponds to _traj->current_frame() of ith obs.
  Measurement get(const uint index) {return _obsVector[index];}
  Measurement obs(const uint index) {return get(index);}
  void push_back(const Measurement datum) {return _obsVector.push_back(const Measurement datum);}
  std::vector<Measurement> operator()(void) {return _obsVector;}
  Measurement operator(const uint index) {return get(index);}
  pTraj getTraj(void) {return _traj;}
  void frame(const uint index) {_traj->readFrame(index);}
};

}  // namespace loos
#endif