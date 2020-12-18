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

#include <Trajectory.hpp>
#include <exceptions.hpp>
#include <iostream>
#include <loos_defs.hpp>
#include <map>
#include <stdexcept>
#include <string>

namespace loos {

class Weights {
public:
  uint current_frame;

protected:
  std::vector<double> _weights;
  pTraj _traj;
  uint _num_weights;
  double _total;
  double _totalTraj;

public:
  // all of these public methods have a definition
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
  virtual const double operator()();
  virtual const double operator()(const uint index);
  virtual void operator()(double newWeight);
  virtual void operator()(double newWeight, const uint index);
  virtual void operator()(std::vector<double> &newWeights);

  virtual std::vector<double> weights();
  virtual void addTraj(pTraj &traj);

public:
  Weights(const std::vector<double> &weightsvec, pTraj &traj)
      : current_frame(0), _weights(weightsvec), _traj{traj},
        _num_weights(weightsvec.size()), _total(0.0) {};

  Weights(const std::vector<double> &weightsvec)
      : current_frame(0), _weights(weightsvec), _num_weights(weightsvec.size()),
        _total(0.0){};
  //! mostly here for function-based weights instances, such as UniformWeight
  //! (constant function).
  Weights(pTraj &traj)
      : current_frame(0), _traj(traj), _num_weights{traj->nframes()},
        _total(0.0){};

  Weights() : current_frame(0), _total(0.0){};

  // define virtual destructor inline to ensure vtable gets made correctly.
  virtual ~Weights() {}
};

} // namespace loos

#endif
