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

#include <Weights.hpp>

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
        // TODO : Really should check for underflow to prevent div by 0
        for (uint i=0; i<_weights.size(); ++i) {
            _weights[i] /= sum;
        }
    }

    //! Return the weight for the current frame of the trajectory
    double Weights::get() {
        current_frame = _traj->currentFrame();
        return _weights.at(current_frame);
    }

    //! Return the weight for frame index of the trajectory
    double Weights::get(uint index) {
        return _weights.at(index);
    }

    double Weights::operator()() {
        return get();
    }

    double Weights::operator()(uint index) {
        return get(index);
    }
}
