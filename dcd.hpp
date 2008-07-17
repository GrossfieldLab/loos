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



#if !defined(DCD_HPP)
#define DCD_HPP


#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <vector>
#include <boost/utility.hpp>

#include <loos_defs.hpp>

#include <AtomicGroup.hpp>
#include <StreamWrapper.hpp>
#include <Trajectory.hpp>

using namespace std;

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
 *  - Does NOT handle ENDIAN issues (i.e. auto-swab)
 *
 *  - Does NOT support fixed atoms
 *
 *  - Does NOT support velocity format
 *
 *  - Reorders the crystal parameters (if present) so they are in
 *    a more sensible order (i.e. a, b, c, alpha, beta, gamma)
 *    
 *  - [Almost] everything returned is a copy
 */
class DCD : public Trajectory {

  // Use a union to convert data to appropriate type...
  typedef union { unsigned int ui; int i; char c[4]; float f; } DataOverlay; 


public:

  // Error classes that we may or may not return...

  //! Error while reading the F77 guard data
  struct RecordError : public exception { virtual const char *what() const throw() { return("Error while reading F77 record"); }; } record_error;
  //! General error while parsing the DCD header
  struct HeaderError : public exception { virtual const char *what() const throw() { return("Error while reading DCD header"); }; } header_error;
  //! General error while reading an F77 data record
  struct LineError : public exception { virtual const char *what() const throw() { return("Error while reading F77 data line"); }; } line_error;
  //! Unexpected EOF
  struct EndOfFile : public exception { virtual const char *what() const throw() { return("Unexpected end of file"); }; } end_of_file;
  //! Endianness of file doesn't match local architecture...
  struct EndianMismatch : public exception { virtual const char *what() const throw() { return("Endianness of file does not match local architecture"); }; } endian_mismatch;
  //! General error...
  struct GeneralError : public exception {
    GeneralError(const char *s) : msg(s) { };
    virtual const char *what() const throw() { return(msg); };
    const char *msg;
  };


  explicit DCD() : Trajectory(), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { }

  //! Begin reading from the file named s
  explicit DCD(const string s) :  Trajectory(s), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { readHeader(); }

  //! Begin reading from the file named s
  explicit DCD(const char* s) :  Trajectory(s), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { readHeader(); }

  //! Begin reading from the stream ifs
  explicit DCD(fstream& fs) : Trajectory(fs), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { readHeader(); };

  //! Read in the header from the stored stream
  void readHeader(void);
  //! Read in the header from the specified stream
  void readHeader(fstream& ifs);

  //! Read the next frame.  Returns false if at EOF
  virtual bool readFrame(void);
  //! Read the ith frame.  Returns false if there is a problem.
  virtual bool readFrame(const unsigned int i);

  //! Rewind the file to the first DCD frame.
  virtual void rewind(void);

  // Accessor methods...

  virtual uint natoms(void) const { return(_natoms); }
  virtual bool hasPeriodicBox(void) const { return(_icntrl[10] == 1); }
  virtual GCoord periodicBox(void) const { return(GCoord(qcrys[0], qcrys[1], qcrys[2])); }

  vector<string> titles(void) const { return(_titles); }

  int icntrl(const int i) const { assert(i>=0 && i<20); return(_icntrl[i]); }
  void icntrl(const int i, const int val) { assert(i>=0 && i<20); _icntrl[i] = val; }

  // * legacy *
  vector<double> crystalParams(void) const { return(qcrys); }
  bool hasCrystalParams(void) const { return(_icntrl[10] == 1); }

  virtual float timestep(void) const { return(_delta); }
  virtual uint nframes(void) const { return(_icntrl[3]); }

  //! Return the raw coords...
  vector<dcd_real> xcoords(void) const { return(xcrds); }
  //! Return the raw coords...
  vector<dcd_real> ycoords(void) const { return(ycrds); }
  //! Return the raw coords...
  vector<dcd_real> zcoords(void) const { return(zcrds); }

  // The following track CHARMm names (more or less...)
  unsigned int nsteps(void) const { return(_icntrl[3]); }
  float delta(void) const { return(_delta); }
  int nsavc(void) const { return(_icntrl[2]); }
  int nfile(void) const { return(_icntrl[0]); }
  int nfixed(void) const { return(_icntrl[8]); }

  //! Auto-interleave the coords into a vector of GCoord()'s.
  /*!  This can be a pretty slow operation, so be careful. */
  virtual vector<GCoord> coords(void);

  //! Interleave coords, selecting entries indexed by map
  // This is slated to go away...
  vector<GCoord> mappedCoords(const vector<int>& map);


  //! Update an AtomicGroup coordinates with the currently-read frame.
  /** This assumes that that atomid's of the
   *AtomicGroup are indices into the DCD frame and are indexed +1,
   *i.e. atomid 7 refers to DCD coords at index 6...
   *
   *There is support for pediodic boundary conditions.  If the DCD has
   *xtal data, then the a, b, and c values are used to update the
   *AtomicGroup::periodicBox().
   */
  virtual void updateGroupCoords(AtomicGroup& g);

private:
  void allocateSpace(const int n);
  void readCrystalParams(void);
  void readCoordLine(vector<float>& v);

  bool endianMatch(StreamWrapper& fsw);

  // For reading F77 I/O
  unsigned int readRecordLen(StreamWrapper& fsw);
  DataOverlay* readF77Line(StreamWrapper& fsw, unsigned int *len);



private:
  int _icntrl[20];          // DCD header data
  int _natoms;              // # of atoms
  vector<string> _titles;   // Vector of title lines from DCD
  vector<double> qcrys;     // Crystal params
  float _delta;             // Timestep (extracted from _icntrl)

  unsigned long frame_size;          // *Internal* size (in bytes) of each frame
  unsigned long first_frame_pos;     // *Internal* location in file of start of data frames
  
  vector<dcd_real> xcrds, ycrds, zcrds;
  
};




#endif
