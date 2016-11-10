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

	//! Combine multiple trajectories (pTraj's) into one large virtual trajectory
	/**
	 * This class can be used just about anywhere a regular Trajectory/pTraj can be used.
	 * Note that the skip and stride settings are applied to each sub-trajectory (as opposed
	 * to the composite trajectory).  They are also set ONLY at instantiation.
	 *
	 */
	class MultiTrajectory : public Trajectory {
	public:
		typedef std::pair<uint, uint>   Location;

		MultiTrajectory()
			: _nframes(0), _skip(0), _stride(1), _curtraj(0), _curframe(0)
		{ }

		//! instantiate a new empty MultiTrajectory
		MultiTrajectory(const AtomicGroup& model)
			: _nframes(0), _skip(0), _stride(1), _curtraj(0), _curframe(0), _model(model)
		{ }

		MultiTrajectory(const AtomicGroup& model, const uint skip, const uint stride)
			: _nframes(0), _skip(skip), _stride(stride), _curtraj(0), _curframe(0), _model(model)
		{ }


		//! Instantiate a new MultiTrajectory using the passed filenames
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


		//! Add a trajectory (by filename)
		void addTrajectory(const std::string& filename) {
			pTraj traj = createTrajectory(filename, _model);
			_trajectories.push_back(traj);
			if (traj->nframes() > _skip)
				_nframes += (traj->nframes() - _skip) / _stride;
		}


		virtual std::string description() const { return("virtual-trajectory"); }

		virtual uint natoms() const { return(_model.size()); }

		//! Total number of frames in composite trajectory
		virtual uint nframes() const { return(_nframes); }

		//! Number of frames in the ith trajectory
		uint nframes(const uint i) const {
			if (i >= _trajectories.size())
				throw(LOOSError("Requesting trajectory size for non-existent trajectory in MultiTraj"));

			if (_trajectories[i]->nframes() <= _skip)
				return 0;
			return( (_trajectories[i]->nframes() - _skip + _stride - 1) / _stride );
		}

		//! Number of trajectories contained
		uint size() const { return(_trajectories.size()); }

		//! Access the individual trajectories
		pTraj operator[](const uint i) const {
			if (i >= _trajectories.size())
				throw(LOOSError("MultiTraj trajectory index out of bounds"));
			return(_trajectories[i]);
		}

		//! Ignore timesteps (for now)
		virtual float timestep() const { return(0.0); }

		//! Whether or not the current sub-trajectory has a periodic box
		virtual bool hasPeriodicBox() const {
			uint i = eof() ? _trajectories.size()-1 : _curtraj;
			return(_trajectories[i]->hasPeriodicBox());
		}

		//! The periodic box of the current sub-trajectory
		virtual GCoord periodicBox() const {
			uint i = eof() ? _trajectories.size()-1 : _curtraj;
			return(_trajectories[i]->periodicBox());
		}

		//! Whether or not the current sub-trajectory has a periodic box
		virtual bool hasVelocities() const {
			uint i = eof() ? _trajectories.size()-1 : _curtraj;
			return(_trajectories[i]->hasVelocities());
		}


		//! Coordinates from the most recently read frame
		virtual std::vector<GCoord> coords() {
			uint i = eof() ? _trajectories.size()-1 : _curtraj;
			return(_trajectories[i]->coords());
		}



		//! Index into the trajectory list for the trajectory currently used
		uint currentTrajectoryIndex() const { return(_curtraj); }

		//! Raw index into the current trajectory for the current frame (i.e. with skip & stride applied)
		uint currentFrameIndex() const { return(_curframe); }

		Location frameIndexToLocation(const uint i);

		bool eof() const {
			return _curtraj >= _trajectories.size();
		}

	private:

		virtual void rewindImpl();
		virtual void seekNextFrameImpl();
		virtual void seekFrameImpl(const uint i);
		virtual bool parseFrame();
		virtual void updateGroupCoordsImpl(AtomicGroup& g);
		virtual void updateGroupVelocitiesImpl(AtomicGroup& g);

		void findNextUsableTraj();


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
