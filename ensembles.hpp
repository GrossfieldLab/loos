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



#if !defined(LOOS_ENSEMBLES_HPP)
#define LOOS_ENSEMBLES_HPP

#include <vector>
#include <boost/tuple/tuple.hpp>

#include <loos_defs.hpp>
#include <MatrixImpl.hpp>
#include <MatrixOps.hpp>

#include <ProgressCounters.hpp>
#include <ProgressTriggers.hpp>

#include <AtomicGroup.hpp>
#include <Trajectory.hpp>

namespace loos {
  class XForm;

  //! Compute the average structure of a set of AtomicGroup objects
  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble);

  //! Compute the average structure of a set of AtomicGroup objects
  /**
   * Takes into consideration the passed set of transforms...
   */
  AtomicGroup averageStructure(const std::vector<AtomicGroup>& ensemble, const std::vector<XForm>& xforms);

  //! Compute the average structure from a trajectory reading only certain frames
  /**
   * Note that the trajectory is NOT stored in memory.  Frames will be read as needed.
   * The trajectory "iterator" will be left pointing to the next frame after the last
   * frame index passed.
  */
  AtomicGroup averageStructure(const AtomicGroup&, const std::vector<XForm>&, pTraj& traj, const std::vector<uint>& indices);

  //! Compute the average structure using all frames in a trajectory
    /**
     * This version only reads a frame at a time from the trajectory.  The trajectory
     * iterator will be left pointing to the end of the trajectory.
    */
  AtomicGroup averageStructure(const AtomicGroup&, const std::vector<XForm>&, pTraj& traj);

#if !defined(SWIG)
  //! Compute an iterative superposition (a la Alan)
  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(std::vector<AtomicGroup>& ensemble, greal threshold = 1e-6, int maxiter=1000);

  //! Compute an iterative superposition by reading in frames from the Trajectory.
  /**
   * The iterativeAlignment() functions that take a trajectory as an argument do
   * NOT cache frames of the trajectory internally.  This means that the trajectory
   * will be read as many times as is necessary for the alignment to
   * converge.  In practice, the OS-specific caching will likely result
   * in decent performance.  If speed is essential, then consider
   * using the iterativeAlignment() version that takes a
   * \p std::vector<AtomicGroup>& as argument instead.
   */
  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj& traj, const std::vector<uint>& frame_indices, greal threshold = 1e-6, int maxiter=1000);

    //! Compute an iterative superposition for an entire trajectory
  boost::tuple<std::vector<XForm>, greal, int> iterativeAlignment(const AtomicGroup& g, pTraj& traj, greal threshold = 1e-6, int maxiter=1000);
#endif // !defined(SWIG)

  void applyTransforms(std::vector<AtomicGroup>& ensemble, std::vector<XForm>& xforms);

  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory);
  void readTrajectory(std::vector<AtomicGroup>& ensemble, const AtomicGroup& model, pTraj trajectory, std::vector<uint>& frames);




  
#if !defined(SWIG)
  RealMatrix extractCoords(const std::vector<AtomicGroup>& ensemble);
  RealMatrix extractCoords(const std::vector<AtomicGroup>& ensemble, const std::vector<XForm>& xforms);

  void subtractAverage(RealMatrix& M);

  //! Compute the SVD of an ensemble with optional alignment (note RSVs returned are transposed)
  /**
   * Returns the U, S, and V' of the SVD of the passed ensemble.  If align is true, then the ensemble
   * is iteratively aligned prior to computing the SVD.
   */
  boost::tuple<RealMatrix, RealMatrix, RealMatrix> svd(std::vector<AtomicGroup>& ensemble, const bool align = true);



#endif   // !defined(SWIG)



  std::vector< std::vector<double> > readCoords(AtomicGroup& model, pTraj& traj,
                                                const std::vector<uint>& indices,
                                                const bool updates);
};



#endif
