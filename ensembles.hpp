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



#if !defined(ENSEMBLES_HPP)
#define ENSEMBLES_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include <MatrixImpl.hpp>


namespace loos {
  class XForm;

  typedef Math::Matrix<float, Math::ColMajor> RealMatrix;
  typedef Math::Matrix<double, Math::ColMajor> DoubleMatrix;

  //! Compute the average structure of a set of AtomicGroup objects
  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble);

  //! Compute the average structure by reading through a trajectory
  AtomicGroup averageStructure(const AtomicGroup&, const std::vector<XForm>&, pTraj);

  //! Compute an iterative superposition (a la Alan)
  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(std::vector<AtomicGroup>& ensemble, greal threshold = 1e-6, int maxiter=1000);

  //! Compute an iterative superposition by reading in frames from the Trajectory.
  /*!
   * This function will internally cache an AtomicGroup copy for each
   * frame of the trajectory.  This could chew up a lot of memory, but
   * we make the assumption that you will usually be aligning against
   * a fairly small subset of each frame...
   */
  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj, greal threshold = 1e-6, int maxiter=1000);


  void applyTransforms(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms);

  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory);
  RealMatrix extractCoords(std::vector<AtomicGroup>& ensemble);
  RealMatrix extractCoords(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms);
  void subtractAverage(RealMatrix& M);
  
  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(RealMatrix& M);
  boost::tuple<DoubleMatrix, DoubleMatrix, DoubleMatrix> svd(DoubleMatrix& M);

  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(std::vector<AtomicGroup>& ensemble);

};



#endif
