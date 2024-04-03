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

#if !defined(LOOS_TRAJECTORY_HPP)
#define LOOS_TRAJECTORY_HPP

#include <istream>
#include <string>
#include <stdexcept>
#include <vector>

#include <boost/utility.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/filesystem.hpp>

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>

#include <AtomicGroup.hpp>


namespace loos {


	//! Base-class for polymorphic trajectories.
	/** This is the interface for trajectories in LOOS.  It is expected
	 *  that at least one frame of coordinates will be buffered internally
	 *  at any given time.
	 *
	 *  Additionally, this class is not designed to provide output, only
	 *  input...
	 *
	 *  +IMPORTANT NOTE+
	 *  The derived classes MUST read in and cache the first frame as
	 *  part of their initialization.  This prevents problems where
	 *  updateGroupCoords() is called prior to the class reading any
	 *  trajectory data (which can occur with some formats, such as
	 *  DCD's, that only have to read a header to configure the internal
	 *  data...  However, just inserting a readFrame(0) in the
	 *  constructor will leave the trajectory iterator in an incorrect
	 *  state--the first call to readFrame() will return the 2nd frame,
	 *  not the first, which is probably not the desired behavior.  The
	 *  derived class must also then set the cached_first flag to true
	 *  after the readFrame(0).  See the DCD class for an example of
	 *  this.
	 */

	class Trajectory  {
	public:
		typedef boost::shared_ptr<std::istream>      pStream;


		Trajectory() : cached_first(false), _filename("unset"), _current_frame(0) { }

		//! Automatically open the file named \a s
		Trajectory(const std::string& s)
			: cached_first(false), _filename(s), _current_frame(0)
		{
			setInputStream(s);
		}

		//! Open using the given stream...
		Trajectory(std::istream& fs) : cached_first(false), _filename("istream"), _current_frame(0)
		{
			setInputStream(fs);
		}


		Trajectory(const Trajectory& t) : ifs(t.ifs), cached_first(t.cached_first), _filename(t._filename), _current_frame(t._current_frame)
		{
		}

		virtual ~Trajectory() { }


		//! Return a string describing trajectory format
		virtual std::string description() const { return("No Description Available"); }

		//! Return the stored filename
		virtual std::string filename() const { return(_filename); }

		//! # of atoms per frame
		virtual uint natoms(void) const =0;
		//! Timestep per frame
		virtual float timestep(void) const =0;
		//! Number of frames in the trajectory
		virtual uint nframes(void) const =0;

		//! Whether or not the trajectory format supports velocities
		virtual bool hasVelocities() const { return(false); }

		//! Conversion applied to velocities to get to \AA/ps
		virtual double velocityConversionFactor() const { return(1.0); }

		//! Rewinds the readFrame() iterator
		bool rewind(void) {
			cached_first = true;
			rewindImpl();
			_current_frame = 0;
			return(parseFrame());
		}

		//! Tests whether or not the given frame/trajectory has periodic
		//! boundary information.
		/** The presence of periodic box information does not necessarily
		 * indicate that said information has been read in yet.  For
		 * example, the presence of crystal data is in the header so this
		 * can be detected before any frame is read, but the crystal data
		 * itself is only read when a frame is read in.
		 */
		virtual bool hasPeriodicBox(void) const =0;
		//! Returns the periodic box for the current frame/trajectory
		virtual GCoord periodicBox(void) const =0;

		//! Returns the current frames coordinates as a vector of GCoords
		/** Some formats, notably DCDs, do not interleave their
		 * coordinates.  This means that this could be a potentially
		 * expensive operation.
		 */
		virtual std::vector<GCoord> coords(void) const =0;

		//! Update the coordinates in an AtomicGroup with the current frame.
		/** The Atom::index() property is used as an index into the
		 * current frame for retrieving coordinates.  The index property
		 * is typically set when a model is read in by the corresponding
		 * class (e.g. PDB, Amber, etc) and is determined by the ordering
		 * of the atoms in that file.  This is not the same as an atomid.
		 * For example, the first atom in a structure could have an atomid
		 * of 1000, but its index would be 0.
		 *
		 * updateGroupCoords() normally assumes that the passed
		 * AtomicGroup has valid indices.  As a safety check, the
		 * createTrajectory() function will check that the AtomicGroup
		 * passed to it has indices.  Only the first atom is tested,
		 * unless the DEBUG compile-flig is turned on, which will force updateGroupCoords() to
		 * validate the entire AtomicGroup every time, with
		 * correspondingly poorer performance.
		 *
		 * Also note that the declaration of this function has changed in
		 * release 2.1.0.  It is no longer virtual, using the NVI-idiom
		 * instead.  Derived classes should override the
		 * updateGroupCoordsImpl() function.
		 */
		void updateGroupCoords(AtomicGroup& g)
		{
#if defined(DEBUG)
			if (! g.allHaveProperty(Atom::indexbit))
				throw(LOOSError("Atoms in AtomicGroup have unset index properties and cannot be used to read a trajectory."));
#else
			if (! g.empty())
				if (! g[0]->checkProperty(Atom::indexbit))
					throw(LOOSError("Atoms in AtomicGroup have unset index properties and cannot be used to read a trajectory."));
#endif

			updateGroupCoordsImpl(g);
		}



		//! Returns the current frame's velocities as a vector of GCoords
		/**
		 * If the trajectory format supports velocities "natively", then those will
		 * be returned.  If not, the velocities will be assumed to be stored in the
		 * coordinates.  Those will be scaled by the velcotiyConversionFactor() and
		 * then returned.
		 */
		virtual std::vector<GCoord> velocities(void) const {
			if (hasVelocities())
				return(velocitiesImpl());
			else
			{
				std::vector<GCoord> vels = coords();
				for (uint i=0; i<vels.size(); ++i)
					vels[i] *= velocityConversionFactor();
				return(vels);
			}
		}



		void updateGroupVelocities(AtomicGroup& g) {
#if defined(DEBUG)
			if (! g.allHaveProperty(Atom::indexbit))
				throw(LOOSError("Atoms in AtomicGroup have unset index properties and cannot be used to read a trajectory."));
#else
			if (! g.empty())
				if (! g[0]->checkProperty(Atom::indexbit))
					throw(LOOSError("Atoms in AtomicGroup have unset index properties and cannot be used to read a trajectory."));
#endif

			if (hasVelocities())
				updateGroupVelocitiesImpl(g);
			else
				g.copyVelocitiesWithIndex( velocities() );
		}


		//! Seek to the next frame in the sequence (used by readFrame() when
		//! operating as an iterator).
		void seekNextFrame(void) {
			cached_first = false;
			++_current_frame;
			if (!atEnd())
				seekFrameImpl(_current_frame);
		}

		//! Seek to a specific frame, be it in the same contiguous file or
		//! in separate files.
		void seekFrame(const uint i) {
			cached_first = false;
			_current_frame = i;
			seekFrameImpl(i);
		}

		//! Parse an actual frame.
		/** parseFrame() is expected to read in a frame through the
		 * Trajectory's StreamWrapper.  It returns a bool indicating whether
		 * or not it was able to actually read a frame (i.e. false indicates
		 * EOF).
		 */
		virtual bool parseFrame(void) =0;

		//! Reads the next frame in a trajectory, returning false if at the end.
		bool readFrame(void) {
			bool b = true;

			if (atEnd())
				return false;

			if (!cached_first) {
				seekNextFrame();
				b = parseFrame();
			} else
				cached_first = false;

			return(b);
		}

		//! Reads a specific frame in a trajectory.
		/** Reading a specific frame also resets the readFrame() iterator
		 * version so it will continue where readFrame(i) left off...
		 */
		bool readFrame(const int i) {
			bool b = true;

			if (!(i == 0 && cached_first)) {
				seekFrame(i);
				b = parseFrame();
			}
			cached_first = false;
			return(b);
		}

		bool atEnd() const {
			return(_current_frame >= nframes());
		}

		uint currentFrame() const {
			return(_current_frame);
		}

	protected:
		void setInputStream(const std::string& fname)
		{
			_filename = fname;
			// Check if the file is 0-sized
			uintmax_t fileSize;
			try {
				fileSize = boost::filesystem::file_size(fname);
			} catch (boost::filesystem::filesystem_error& e) {
				throw(FileOpenError(fname, std::string("Error computing file size: ") + e.what()));
			}
			if (fileSize == 0) {
				throw(FileReadError(fname, std::string("File is empty")));
			}
			
			ifs = pStream(new std::fstream(fname.c_str(), std::ios_base::in | std::ios_base::binary));
			if (!ifs->good())
				throw(FileOpenError(fname));
		}


		void setInputStream(std::istream& fs)
		{
			_filename = "istream";
			ifs = pStream(&fs, boost::lambda::_1);    // lambda function makes a NOOP deallocator
		}




		pStream ifs;
		bool cached_first;    // Indicates that the first frame is cached by
		// the subclass...

		std::string _filename;   // Remember filename (if passed)
		uint _current_frame;

	private:

		//! NVI implementation for seeking next frame
		virtual void seekNextFrameImpl() =0;

		//! NVI implementation for seeking a specific frame
		virtual void seekFrameImpl(const uint) =0;

		//! NVI implementation of rewind
		virtual void rewindImpl(void) =0;

		//! NVI implementation of updateGroupCoords() for derived classes to override
		virtual void updateGroupCoordsImpl(AtomicGroup& g) =0;

		virtual void updateGroupVelocitiesImpl(AtomicGroup& g) {
			throw(LOOSError("No velocity update implementation defined but trajectory supports it"));
		}

		virtual std::vector<GCoord> velocitiesImpl() const { return(std::vector<GCoord>()); }

	};

}

#endif
