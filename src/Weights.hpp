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

#if !defined(LOOS_WEIGHTS_HPP)
#define LOOS_WEIGHTS_HPP

#include <loos_defs.hpp>
#include <Trajectory.hpp>
#include <exceptions.hpp>
#include <iostream>
#include <string>
#include <stdexcept>
#include <map>

namespace loos {

    class Weights {
    public:
        Weights(const std::string &filename, pTraj const traj ):
                                        current_frame(0),
                                        _filename(filename)
                                       {
            add_traj(traj);
        };

        Weights(const std::string &filename): current_frame(0),
                                             _filename(filename) {

        }

        Weights() {

        };
        ~Weights() { };

        const double get();
        const double get(const uint index);
        uint current_frame;
        uint size();

        void normalize();
        void accumulate();
        void accumulate(const uint index);
        const double totalWeight();
        const double trajWeight();
        void add_traj(pTraj const traj);
        void new_traj(pTraj const traj);
        const double operator()();
        const double operator()(const uint index);

        std::vector<double> weights();




    private:
        uint read_weights(const std::string &filename);
        uint read_weights_list(const std::string &filename);
        uint _num_weights;
        pTraj _traj;
        std::vector<double> _weights;
        std::map<std::string, std::string> _weights_files;
        double _total;
        double _totalTraj;
        std::string _filename;
        bool _has_list;
    };

}

#endif
