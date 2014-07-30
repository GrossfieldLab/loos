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


#if !defined(LOOS_XTC_HPP)
#define LOOS_XTC_HPP
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


  //! Class representing GROMACS reduced precision, compressed trajectories
  /**
   * The XTC format does not have even-sized frames nor does it have
   * any kind of frame index or number of frames metadata.  The entire
   * trajectory has to be scanned in order to determine the number of
   * frames and to build an index that allows seeking to specific
   * frames.  This is done by reading only enough of each frame header
   * to permit building the index, so it should be a pretty fast
   * operation.
   */
  class XTC : public Trajectory {

    // Systems with this size or smaller will not be compressed.
    static const uint min_compressed_system_size;
      

    // The frame header
    struct Header {
      uint magic;
      uint natoms, step;
      float time, box[9];
    };

    // Globals required by the decoding routines...
    static const int magicints[];
    static const int firstidx, lastidx;

    static const int magic;


    //! Size of the data stored in the XTC file
    typedef float    xtc_t;

  public:
    explicit XTC(const std::string& s) : Trajectory(s), xdr_file(ifs()),natoms_(0) {
      init();
    }

    explicit XTC(const char* p) : Trajectory(p), xdr_file(ifs()),natoms_(0) {
      init();
    }

    explicit XTC(std::istream& is) : Trajectory(is), xdr_file(ifs()), natoms_(0) {
      init();
    }



    uint natoms(void) const { return(natoms_); }
    float timestep(void) const { return(timestep_); }
    uint nframes(void) const { return(frame_indices.size()); }
    bool hasPeriodicBox(void) const { return(true); }
    GCoord periodicBox(void) const { return(box); }


    std::vector<GCoord> coords(void) { return(coords_); }

    //! Return the stored file's precision
    double precision(void) const { return(precision_); }

  private:

    void init(void) {
      scanFrames();
      coords_.reserve(natoms_);
      if (!parseFrame())
        throw(std::logic_error("Unable to read in the first frame"));
      cached_first = true;
    }      

    internal::XDRReader xdr_file;
    std::vector<size_t> frame_indices;
    uint natoms_;
    GCoord box;
    double precision_;
    std::vector<GCoord> coords_;
    double timestep_;
    Header current_header_;
    
    bool parseFrame(void);

  private:

    int sizeofint(int);
    int sizeofints(uint*, const uint);
    int decodebits(int*, uint);
    void decodeints(int*, const int, int, uint*, int*);
    bool readFrameHeader(Header&);
    void scanFrames(void);
    
    void seekNextFrameImpl(void) { }
    void seekFrameImpl(uint);
    void rewindImpl(void) { ifs()->clear(); ifs()->seekg(0); }
    void updateGroupCoordsImpl(AtomicGroup& g);
    bool readCompressedCoords(void);
    bool readUncompressedCoords(void);
  };

}



#endif
