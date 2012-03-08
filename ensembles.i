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



%header %{
#include <ensembles.hpp>
%}

namespace loos {

  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble);

  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble, const std::vector<XForm>& xforms);

  AtomicGroup averageStructure(const AtomicGroup&, const std::vector<XForm>&, pTraj& traj, const std::vector<uint>& indices);

  AtomicGroup averageStructure(const AtomicGroup&, const std::vector<XForm>&, pTraj& traj);

  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(std::vector<AtomicGroup>& ensemble, greal threshold = 1e-6, int maxiter=1000);

  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj& traj, const std::vector<uint>& frame_indices, greal threshold = 1e-6, int maxiter=1000);

  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj& traj, greal threshold = 1e-6, int maxiter=1000);


  void applyTransforms(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms);

  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory);
  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory, std::vector<uint>& frames);
};

