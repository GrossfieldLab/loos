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

namespace loos {
    //! Weights class to handle reweighting values computed from a trajectory
    uint Weights::read_weights(const std::string& filename)
    {
        std::ifstream ifs(filename.c_str());
        if (!ifs) {
            throw(FileOpenError(filename));
        }

        std::string input;
        while (getline(ifs, input)) {
            // skip blank lines and comments beginning with "#"
            if ( (input.length() == 0) || (input[0] == '#' ) ) {
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


    //! Normalize the weights so they sum to 1
    void Weights::normalize() {
        double sum = 0.0;
        for (uint i=0; i<_weights.size(); ++i) {
            sum += _weights[i];
        }
        std::cerr << "Total weight = " << sum << std::endl;
        // TODO : Really should check for underflow to prevent div by 0
        for (uint i=0; i<_weights.size(); ++i) {
            _weights[i] /= sum;
        }
    }

    //! Keep track of total weight used
    void Weights::accumulate() {
        _total += _weights.at(_traj->currentFrame());
    }

    void Weights::accumulate(const uint index) {
        _total += _weights.at(index);
    }

    //! Return the totalWeight, as tracked using accumulate
    const double Weights::totalWeight() {
        return _total;
    }

    //! Add trajectory to class and verify size match
    void Weights::add_traj(pTraj const traj) {
        _traj = traj;
        _num_weights = read_weights(_filename);
        // # of weights must match number of frames in the associated traj
        if (_num_weights != _traj->nframes()) {
            throw(LOOSError(std::string("Number of weights must match the length of the trajectory")));
        }

    }

    //! Return the weight for the current frame of the trajectory
    const double Weights::get() {
        current_frame = _traj->currentFrame();
        return _weights.at(current_frame);
    }

    //! Return the weight for frame index of the trajectory
    const double Weights::get(const uint index) {
        return _weights.at(index);
    }

    //! calling nomenclature wraps get
    const double Weights::operator()() {
        return get();
    }

    const double Weights::operator()(const uint index) {
        return get(index);
    }

    //! Return the number of weights
    uint Weights::size() {
        return _num_weights;
    }

    //! Return the vector of weights
    std::vector<double> Weights::weights() {
        return _weights;
    }
}
