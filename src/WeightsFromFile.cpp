
/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2020, Tod D. Romo, Alan Grossfield, Louis Smit
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

#include "WeightsFromFile.hpp"
#include <sstream>

namespace loos {
//! WeightsFromFile class to handle reweighting values computed from a
//! trajectory
uint WeightsFromFile::read_weights(const std::string &filename) {
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
uint WeightsFromFile::read_weights_list(const std::string &filename) {
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

//! Add trajectory to class and verify size match with existing WeightsFromFile
void WeightsFromFile::add_traj(pTraj &traj) {
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

} // namespace loos