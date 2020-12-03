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
uint Weights::read_weights(const std::string &filename) {
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    std::cerr << "Cannot open weights file: " << filename << std::endl;
    throw(FileOpenError(filename));
  }

  std::string input;
  while (getline(ifs, input)) {
    // skip blank lines and comments beginning with "#"
    if ((input.length() == 0) || (input[0] == '#')) {
      // do nothing
    }
    // TODO: we should really let it be any column
    else {
      double value = parseStringAs<double>(input);
      _weights.push_back(value);
    }
  }
  return _weights.size();
}

//! Read in a list of files matching weights files to trajectory files
uint Weights::read_weights_list(const std::string &filename) {
  uint num_weights_files = 0;
  _has_list = true;

  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    throw(FileOpenError(filename));
  }

  std::string line, traj_file, weights_file;
  std::istringstream iss;
  while (getline(ifs, line)) {
    // skip blank lines and comments beginning with "#"
    if ((line.length() == 0) || (line[0] == '#')) {
      // do nothing
    } else {
      num_weights_files++;
      iss.str(line);
      iss >> traj_file >> weights_file;
      _weights_files[traj_file] = weights_file;
    }
  }
  return num_weights_files;
}

//! Normalize the weights so they sum to 1
void Weights::normalize() {
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
void Weights::accumulate() {
  _total += _weights.at(_traj->currentFrame());
  _totalTraj += _weights.at(_traj->currentFrame());
}

void Weights::accumulate(const uint index) {
  _total += _weights.at(index);
  _totalTraj += _weights.at(index);
}

//! Return the totalWeight, as tracked using accumulate
const double Weights::totalWeight() { return _total; }

//! Return the weight of the current trajectory
const double Weights::trajWeight() { return _totalTraj; }

//! Add trajectory to class and verify size match with existing Weights
void Weights::add_traj(pTraj &traj) {
  _traj = traj;
  // If we have a list of weights files, read the correct one
  // TODO: need to check to make sure the filename is in the map
  if (_has_list) {
    _filename = _weights_files[_traj->filename()];
  }
  _num_weights = read_weights(_filename);
  // # of weights must match number of frames in the associated traj
  if (_num_weights != _traj->nframes()) {
    throw(LOOSError(std::string(
        "Number of weights must match the length of the trajectory")));

    // Zero out the weight of the trajectory
    _totalTraj = 0.0;
  }
}

//! Return the weight for the current frame of the trajectory
const double Weights::get() {
  current_frame = _traj->currentFrame();
  return _weights.at(current_frame);
}

//! Return the weight for frame index of the trajectory
const double Weights::get(const uint index) { return _weights.at(index); }

//! calling nomenclature wraps get
const double Weights::operator()() { return get(); }

const double Weights::operator()(const uint index) { return get(index); }

//! Bind a new weight to the current frame
void Weights::set(double newWeight) { _weights.at(current_frame) = newWeight; }

//! Bind a new weight to a particular frame
void Weights::set(double newWeight, const uint index) {
  _weights.at(index) = newWeight;
}

//! binding nomenclature wraps set
void Weights::operator()(double newWeight) { set(newWeight); }
void Weights::operator()(double newWeight, const uint index) {
  set(newWeight, index);
}

//! set all weights from passed vector
void Weights::operator()(std::vector<double> &newWeights) {
  if (_weights.size() != newWeights.size())
    throw(LOOSError(std::string(
      "Number of weights in class is " + std::to_string(_weights.size())
      + " number inserted is " + std::to_string(newWeights.size())
      + " these must match."
    )));
  // note that operator= for stl::vectors does in fact copy their contents.
  _weights = newWeights;
}

//! Return the number of weights
uint Weights::size() { return _num_weights; }

//! Return the vector of weights
std::vector<double> Weights::weights() { return _weights; }
} // namespace loos
