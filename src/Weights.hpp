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
        uint current_frame;
    private:
        double _total;

    public:
        virtual const double get();
        virtual const double get(const uint index);
        virtual void set(double newWeight);
        virtual void set(double newWeight, const uint index);
        virtual uint size();

        virtual void normalize();
        virtual void accumulate();
        virtual void accumulate(const uint index);
        virtual const double totalWeight();
        virtual const double trajWeight();
        virtual void add_traj(pTraj&  traj)=0;
        virtual const double operator()();
        virtual const double operator()(const uint index);
        virtual void operator()(double newWeight);
        virtual void operator()(double newWeight, const uint index);
        virtual void operator()(std::vector<double>& newWeights);
        
        virtual std::vector<double> weights();

    private:
        uint _num_weights;
        pTraj _traj;
        std::vector<double> _weights;
        double _totalTraj;

    public:
        Weights(const std::vector<double> &weightsvec, pTraj& traj ):
                                        _weights(weightsvec),
                                        _num_weights(weightsvec.size()),
                                        current_frame(0),
                                        _total(0.0),
                                        _traj{traj} { };

        Weights(const std::vector<double> &weightsvec): 
                                        current_frame(0),                                            
                                        _weights(weightsvec),
                                        _num_weights(weightsvec.size()),
                                        _total(0.0) { };

        Weights() : current_frame(0),
                    _total(0.0) { };

        // define virtual destructor inline to ensure vtable gets made correctly.
        virtual ~Weights() { }

    };

}

#endif
