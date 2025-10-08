/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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


/*
  NOTE:  This code is based on the xdrfile library authored by:
  David van der Spoel <spoel@gromacs.org>
  Erik Lindahl <lindahl@gromacs.org>
  and covered by the GLPL-v3 license
*/


#if !defined(LOOS_TRR_HPP)
#define LOOS_TRR_HPP
#include <ios>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <loos_defs.hpp>

#include <Coord.hpp>
#include <xdr.hpp>
#include <AtomicGroup.hpp>
#include <Trajectory.hpp>

#include <boost/format.hpp>


namespace loos {

	//! Class representing the GROMACS TRR trajectory files
	/**
	 * The TRR format can technically have differing numbers of atoms or
	 * even components (ie coords, velocities, forces) per frame.  This
	 * is supported by LOOS, to an extent.  If the components vary, then
	 * so long as there are coordinates, there should be no problem.  If
	 * the number of atoms vary and you're using the updateGroupCoords()
	 * function, then care must be taken that the group matches the
	 * frame.  Recall that the atomids are actually indices into the
	 * trajectory frame (ie atomid 1 takes its coords from the first
	 * trajectory coordinate, etc).
	 *
	 * The TRR format can be written in either a single or double
	 * precision format.  Internally, LOOS will read in either and
	 * translate it into whatever the default size is for LOOS (ie
	 * GCoord, which by default is a double).  Metadata, such as lambda,
	 * or box sizes, are always upconverted to double by LOOS.
	 *
	 * Since the TRR frame size is not fixed, the entire
	 * trajectory will be quickly scanned to build up an index of where
	 * the frames begin (see the loos::XTC class for more information).
	 *
	 * Finally, note that GROMACS stores data in nm whereas LOOS uses
	 * angstroms, so coordinate/box data will be automatically scaled by
	 * LOOS.
	 */

	class TRR : public Trajectory {
		static const int magic;       // various internal constants
		// GROMACS uses...
		static const int DIM;

		// This is the internal header for a TRR frame
		struct Header {
			bool bDouble;
			int ir_size;
			int e_size;
			int box_size;
			int vir_size;
			int pres_size;
			int top_size;
			int sym_size;
			int x_size;
			int v_size;
			int f_size;

			int natoms;
			int step;
			int nre;
			float tf;
			float lambdaf;
			double td;
			double lambdad;
		};

	public:
		explicit TRR(const std::string& s) : Trajectory(s), xdr_file(ifs.get()) {
			init();
		}

		explicit TRR(const char* p) : Trajectory(p), xdr_file(ifs.get()) {
			init();
		}
		explicit TRR(std::istream& is) : Trajectory(is), xdr_file(ifs.get()) {
			init();
		}

		std::string description() const { return("Gromacs TRR"); }
		static pTraj create(const std::string& fname, const AtomicGroup& model) {
			return(pTraj(new TRR(fname)));
		}


		uint natoms(void) const { return(hdr_.natoms); }

		/**
		 * GROMACS trajectories do not define a timestep value as DCD does
		 */
		float timestep(void) const { return(0); }

		uint nframes(void) const { return(frame_indices.size()); }
		bool hasPeriodicBox(void) const { return(hdr_.box_size != 0); }
		GCoord periodicBox(void) const { return(box); }


		std::vector<GCoord> coords(void) const { return(coords_); }

		// TRR specific attributes...
		std::vector<double> virial(void) const { return(vir_); }
		std::vector<double> pressure(void) const { return(pres_); }
		std::vector<GCoord> forces(void) const { return(forc_); }

		bool isDouble(void) const { return(hdr_.bDouble); }
		bool hasVirial(void) const { return(hdr_.vir_size != 0); }
		bool hasPressure(void) const { return(hdr_.pres_size != 0); }
		bool hasCoords(void) const { return(hdr_.x_size != 0); }
		bool hasVelocities(void) const { return(hdr_.v_size != 0); }
		bool hasForces(void) const { return(hdr_.f_size != 0); }

		double time(void) const { return( hdr_.bDouble ? hdr_.td : hdr_.tf ); }
		double lambda(void) const { return( hdr_.bDouble ? hdr_.lambdad : hdr_.lambdaf); }

		int step(void) const { return(hdr_.step); }

		bool parseFrame(void) {
			if (!readHeader(hdr_))  // Catch EOF
				return(false);

			if (hdr_.bDouble)
				return(readRawFrame<double>());

			return(readRawFrame<float>());
		}


	private:
		void init(void);
		int floatSize(Header& h);
		bool readHeader(Header& h);

		// These simplify reading of blocks of data from the TRR file.
		// Templatized since GROMACS can store data in either single or
		// double-precision...

		template<typename T>
		void readBlock(std::vector<double>& v, const uint n, const std::string& msg) {
			std::vector<T> buf(n);
			uint i = xdr_file.read(buf.data(), n);
			if (i != n)
				throw(FileReadError(_filename, "Unable to read " + msg));

			for (i=0; i<n; ++i)
				v.push_back(buf[i]);
		}


		// This assumes the block of data are triplets and converts them
		// into GCoords, scaling from nm to Angstroms along the way...
		template<typename T>
		void readBlock(std::vector<GCoord>& v, const uint n, const std::string& msg) {
			T* buf = new T[n];
			if (buf == 0)
				throw(std::runtime_error("Out of memory"));

			if (xdr_file.read(buf, n) != n) {
				delete[] buf;
				throw(FileReadError(_filename, "Unable to read " + msg));
			}
			for (uint i=0; i<n; i += DIM)
				v.push_back(GCoord(buf[i], buf[i+1], buf[i+2]) * 10.0);

			delete[] buf;
		}


		// Note: Assumes that the object Header has already been read...
		template<typename T>
		bool readRawFrame() {

			// Clear data first...
			box_.clear();
			vir_.clear();
			pres_.clear();
			velo_.clear();
			forc_.clear();
			coords_.clear();

			if (hdr_.box_size) {
				readBlock<T>(box_, DIM*DIM, "box");
				box = GCoord(box_[0], box_[4], box_[8]) * 10.0;   // Convert
				// to angstroms
			}

			if (hdr_.vir_size)
				readBlock<T>(vir_, DIM*DIM, "virial");
			if (hdr_.pres_size)
				readBlock<T>(pres_, DIM*DIM, "pressure");

			if (hdr_.x_size)
				readBlock<T>(coords_, hdr_.natoms * DIM, "Coordinates");

			if (hdr_.v_size)
				readBlock<T>(velo_, hdr_.natoms * DIM, "Velocities");

			if (hdr_.f_size)
				readBlock<T>(forc_, hdr_.natoms * DIM, "Forces");


			return(! ((xdr_file.get())->fail() || (xdr_file.get())->eof()) );
		}

		void rewindImpl(void) { ifs->clear(); ifs->seekg(0, std::ios_base::beg); }
		void seekNextFrameImpl(void) { }
		void seekFrameImpl(uint);
		void updateGroupCoordsImpl(AtomicGroup& g);
		void updateGroupVelocitiesImpl(AtomicGroup& g);
		std::vector<GCoord> velocitiesImpl() const { return(velo_); }


	private:
		internal::XDRReader xdr_file;
		std::vector<GCoord> coords_;
		GCoord box;
		std::vector<size_t> frame_indices;   // Index into file for start
		// of frame header

		std::vector<double> box_;
		std::vector<double> vir_;
		std::vector<double> pres_;
		std::vector<GCoord> velo_;
		std::vector<GCoord> forc_;

		Header hdr_;
	};



};


#endif

