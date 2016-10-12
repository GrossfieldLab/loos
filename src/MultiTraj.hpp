/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2016, Tod D. Romo, Alan Grossfield
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

#if !defined(LOOS_MULTITRAJ_HPP)
#define LOOS_MULTITRAJ_HPP

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>
#include <sfactories.hpp>




namespace loos {

  class MultiTrajectory : public Trajectory {
  public:

    MultiTrajectory(const AtomicGroup& model)
      : _nframes(0), _skip(0), _stride(1), _curtraj(0), _curframe(0), _model(model)
    { }

    MultiTrajectory(const std::vector<std::string>& filenames,
		     const AtomicGroup& model)
      : _nframes(0), _skip(0), _stride(1), _curtraj(0), _curframe(0), _model(model)
    {
      initWithList(filenames, model);
    }


    MultiTrajectory(const std::vector<std::string>& filenames,
		      const AtomicGroup& model,
		      const uint skip,
		      const uint stride)
      : _nframes(0), _skip(skip), _stride(stride), _curtraj(0), _curframe(skip), _model(model)
    {
      initWithList(filenames, model);
    }


    void addTrajectory(const std::string& filename) {
      pTraj traj = createTrajectory(filename, _model);
      _trajectories.push_back(traj);
      _nframes += (traj->nframes() - _skip) / _stride;
    }

    virtual std::string description() const { return("virtual-trajectory"); }

    virtual uint natoms() const { return(_model.size()); }
    
    virtual uint nframes() const { return(_nframes); }

    uint nframes(const uint i) const { return( (_trajectories[i]->nframes() - _skip + _stride - 1) / _stride ); }

    uint size() const { return(_trajectories.size()); }
    pTraj operator[](const uint i) const { return(_trajectories[i]); }
    
    // Ignore timesteps (for now)
    virtual float timestep() const { return(0.0); }

    virtual bool hasPeriodicBox() const {
      return(_trajectories[_curtraj]->hasPeriodicBox());
    }

    virtual GCoord periodicBox() const {
      return(_trajectories[_curtraj]->periodicBox());
    }

    virtual std::vector<GCoord> coords() {
      return(_trajectories[_curtraj]->coords());
    }

    

    uint currentTrajectoryIndex() const { return(_curtraj); }
    uint currentFrameIndex() const { return(_curframe); }
    
    
  private:

    virtual void rewindImpl();
    virtual void seekNextFrameImpl();
    virtual void seekFrameImpl(const uint i);
    virtual bool parseFrame();
    virtual void updateGroupCoordsImpl(AtomicGroup& g);



    
    // Make these private so you can't accidently try to use them...
    MultiTrajectory(const std::string& s) { }
    MultiTrajectory(const std::istream& fs) { }
    
    void initWithList(const std::vector<std::string>& filenames, const AtomicGroup& model);
    
  private:
    uint _nframes;
    uint _skip, _stride;
    uint _curtraj, _curframe;
    AtomicGroup _model;
    std::vector<pTraj> _trajectories;

  };


}



#endif // !defined(LOOS_MULTITRAJ_HPP)
