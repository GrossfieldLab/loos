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



#if !defined(LOOS_DCD_HPP)
#define LOOS_DCD_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <exception>

#include <loos_defs.hpp>

#include <StreamWrapper.hpp>
#include <Trajectory.hpp>


namespace loos {


  //! Class for reading DCD files
  /**
   *
   *Instantiating a DCD object with either a filename or an fstream
   *only reads the header from the file, not any frames.  When a frame
   *is read, the x,y,z coordinates are stored internally in a vector.
   *This must be copied out to the caller or used to update the
   *coordinates for an AtomicGroup.
   *
   *Notes:
   *
   *  - Does NOT support fixed atoms
   *
   *  - Does NOT support velocity format
   *
   *  - Reorders the crystal parameters (if present) so they are in
   *    a more sensible order (i.e. a, b, c, alpha, beta, gamma)
   *    
   *  - [Almost] everything returned is a copy
   *
   *  - Endian detection is based on the expected size of the header
   */
  class DCD : public Trajectory {
    static bool suppress_warnings;

    
    // Use a union to convert data to appropriate type...
    typedef union { unsigned int ui; int i; char c[4]; float f; } DataOverlay;

  public:
    class EndOfFile : public LOOSError {
    public:
      EndOfFile() : LOOSError("unexpected end of file while reading DCD") { }
    };

    //! Begin reading from the file named s
    explicit DCD(const std::string s) :  Trajectory(s), _natoms(0), _nframes(0),
                                         qcrys(std::vector<double>(6)),
                                         frame_size(0), first_frame_pos(0),
                                         swabbing(false) { initTrajectory(); }

    //! Begin reading from the file named s
    explicit DCD(const char* s) :  Trajectory(s), _natoms(0), _nframes(0), 
                                   qcrys(std::vector<double>(6)), frame_size(0),
                                   first_frame_pos(0), swabbing(false) { initTrajectory(); }

    //! Begin reading from the stream ifs
    explicit DCD(std::istream& fs) : Trajectory(fs), _natoms(0), _nframes(0),
                                     qcrys(std::vector<double>(6)), frame_size(0), first_frame_pos(0),
                                     swabbing(false) { initTrajectory(); };

      std::string description() const { return("CHARMM/NAMD DCD"); }

      static pTraj create(const std::string& fname, const AtomicGroup& model) {
	return(pTraj(new DCD(fname)));
      }



    // Accessor methods...

    virtual uint natoms(void) const;
    virtual bool hasPeriodicBox(void) const;
    virtual GCoord periodicBox(void) const;

    std::vector<std::string> titles(void) const;

    int icntrl(const int) const;
    void icntrl(const int, const int);

    // * legacy *
    std::vector<double> crystalParams(void) const;
    bool hasCrystalParams(void) const;

    virtual float timestep(void) const;
    virtual uint nframes(void) const;

    //! Return the raw coords...
    std::vector<dcd_real> xcoords(void) const;
    //! Return the raw coords...
    std::vector<dcd_real> ycoords(void) const;
    //! Return the raw coords...
    std::vector<dcd_real> zcoords(void) const;

    // The following track CHARMm names (more or less...)
    unsigned int nsteps(void) const;
    float delta(void) const;
    int nsavc(void) const;
    int nfile(void) const;
    int nfixed(void) const;

    //! Returns true if the DCD file being read is in the native endian format
    bool nativeFormat(void) const;

    //! Auto-interleave the coords into a vector of GCoord()'s.
    /*!  This can be a pretty slow operation, so be careful. */
    virtual std::vector<GCoord> coords(void);

    //! Interleave coords, selecting entries indexed by map
    // This is slated to go away...
    std::vector<GCoord> mappedCoords(const std::vector<int>& map);



    static void setSuppression(const bool b) { suppress_warnings = b; }

    //! Parse a frame of the DCD
      virtual bool parseFrame(void);
    
  private:

    //! Read in the header from the stored stream
    void readHeader(void);

    void initTrajectory();

    uint calculateNumberOfFrames();
    
    // Trajectory member functions we must provide...
    virtual void seekNextFrameImpl(void) { }    // DCD frames are always contiguous, so do nothing...
    //! Calculate offset into DCD file for frame and seek to it.
    virtual void seekFrameImpl(const uint);
    
    //! Rewind the file to the first DCD frame.
    virtual void rewindImpl(void);
    
    //! Update an AtomicGroup coordinates with the currently-read frame.
    virtual void updateGroupCoordsImpl(AtomicGroup& g);



    void allocateSpace(const int n);
    bool readCrystalParams(void);
    bool readCoordLine(std::vector<float>& v);

    void endianMatch(StreamWrapper& fsw);

    // For reading F77 I/O
    unsigned int readRecordLen(void);
    DataOverlay* readF77Line(unsigned int *len);



  private:
    int _icntrl[20];          // DCD header data
    uint _natoms;              // # of atoms
    uint _nframes;
    std::vector<std::string> _titles;   // Vector of title lines from DCD
    std::vector<double> qcrys;     // Crystal params
    float _delta;             // Timestep (extracted from _icntrl)

    long frame_size;          // *Internal* size (in bytes) of each frame
    long first_frame_pos;     // *Internal* location in file of start of data frames

    bool swabbing;            // DCD being read is not in native format...
  
    std::vector<dcd_real> xcrds, ycrds, zcrds;
  
  };

}



#endif
