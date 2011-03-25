/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011, Tod D. Romo, Alan Grossfield
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



#if !defined(TRAJECTORY_ITERATOR_HPP)
#define TRAJECTORY_ITERATOR_HPP


#include <Trajectory.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace loos {

  //! An iterable Trajectory class
  /**
   * This class wraps a Trajectory and an AtomicGroup.  It can operate
   * similar to any C++ const iterator.
   */
  class TrajectoryIterator : public boost::iterator_facade<
    TrajectoryIterator,
    AtomicGroup,
    boost::random_access_traversal_tag,
    const AtomicGroup&
    >

  {
  public:

    typedef AtomicGroup                       value_type;
    typedef TrajectoryIterator                const_iterator;


    TrajectoryIterator(AtomicGroup model, pTraj traj) : model_(model.copy()),
                                                        trajectory_(traj),
                                                        current_frame_number_(0)
    {
      trajectory_->rewind();
    }


    TrajectoryIterator(AtomicGroup model, pTraj traj, const long n) : model_(model.copy()),
                                                                      trajectory_(traj),
                                                                      current_frame_number_(n)
    {
      trajectory_->rewind();
    }


    TrajectoryIterator(const TrajectoryIterator& o) : model_(o.model_),
                                                      trajectory_(o.trajectory_),
                                                      current_frame_number_(o.current_frame_number_)
    { }



    const_iterator begin() const {
      const_iterator iter = TrajectoryIterator(*this);
      iter.current_frame_number_ = 0;
      iter.trajectory_->rewind();
      return(iter);
    }


    /**
     * This is likely going to be an expensive operation since it
     * could be called at each loop iteration.  Not sure how best to
     * optimize this, other than to suggest caching the returned
     * iterator if you <em>really</em> need speed...
     */
    const_iterator end() const {
      const_iterator iter = TrajectoryIterator(*this);
      iter.current_frame_number_ = iter.trajectory_->nframes();
      return(iter);
    }


    bool equal(const const_iterator& other) const {
      return(current_frame_number_ == other.current_frame_number_);
    }


    void increment() { 
      ++current_frame_number_;
    }

    void decrement() {
      --current_frame_number_;
    }

    void advance(const long i) {
      current_frame_number_ += i;
    }

    const AtomicGroup& dereference() const {
      if (current_frame_number_ < 0 || current_frame_number_ >= trajectory_->nframes())
        throw(std::range_error("TrajectoryIterator index out of bounds"));

      trajectory_->seekFrame(current_frame_number_);
      trajectory_->parseFrame();
      trajectory_->updateGroupCoords(model_);
      return(model_);
    }

    long distance_to(const TrajectoryIterator& other) const {
      return(other.current_frame_number_ - current_frame_number_);
    }


  private:

    TrajectoryIterator(const long n) : model_(),
                                       trajectory_(),
                                       current_frame_number_(n) { }


    mutable AtomicGroup model_;
    pTraj trajectory_;
    long current_frame_number_;
  };




};




#endif
