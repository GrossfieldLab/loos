/*
  dcd.h
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/

#if !defined(DCD_HPP)
#define DCD_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <vector>
#include <tr1/memory>
#include <boost/utility.hpp>
#include <assert.h>

#include <loos.hpp>
#include <AtomicGroup.hpp>

using namespace std;

//! Class for reading DCD files
/*!

  Instantiating a DCD object with either a filename or an ifstream
  only reads the header from the file, not any frames.  When a frame
  is read, the x,y,z coordinates are stored internally in a vector.
  This must be copied out to the caller or used to update the
  coordinates for an AtomicGroup.

  Notes:

    - Does NOT handle ENDIAN issues (i.e. auto-swab)

    - Does NOT support fixed atoms

    - Does NOT support velocity format

    - Reorders the crystal parameters (if present) so they are in
      a more sensible order (i.e. a, b, c, alpha, beta, gamma)
      
    - [Almost] everything returned is a copy
*/
class DCD : public boost::noncopyable {

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
  //! General error...
  struct GeneralError : public exception {
    GeneralError(const char *s) : msg(s) { };
    virtual const char *what() const throw() { return(msg); };
    const char *msg;
  };


  DCD() : del_stream(false), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { };

  //! Begin reading from the file named s
  DCD(const string s) :  del_stream(true), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) {
    _ifs = new ifstream(s.c_str());
    if (!_ifs)
      throw(runtime_error("Cannot open DCD file " + s));
    readHeader();
  }

  //! Begin reading from the file named s
  DCD(const char* s) :  del_stream(true), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) {
    _ifs = new ifstream(s);
    if (!_ifs)
      throw(runtime_error("Cannot open DCD file " + string(s)));
    readHeader();
  }

  //! Begin reading from the stream ifs
  DCD(ifstream& ifs) : del_stream(false), _ifs(&ifs), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { readHeader(); };

  ~DCD() {
    if (del_stream)
      delete _ifs;
  }


  //! Read in the header from the stored stream
  void readHeader(void);
  //! Read in the header from the specified stream
  void readHeader(ifstream& ifs);

  //! Read the next frame.  Returns false if at EOF
  bool readFrame(void);
  //! Read the ith frame.  Returns false if there is a problem.
  bool readFrame(const unsigned int i);

  // Accessor methods...

  int natoms(void) const { return(_natoms); }

  vector<string> titles(void) const { return(_titles); }

  int icntrl(const int i) const { assert(i>=0 && i<20); return(_icntrl[i]); }
  void icntrl(const int i, const int val) { assert(i>=0 && i<20); _icntrl[i] = val; }

  vector<double> crystalParams(void) const { return(qcrys); }

  float delta(void) const { return(_delta); }

  //! Return the raw coords...
  vector<dcd_real> xcoords(void) const { return(xcrds); }
  //! Return the raw coords...
  vector<dcd_real> ycoords(void) const { return(ycrds); }
  //! Return the raw coords...
  vector<dcd_real> zcoords(void) const { return(zcrds); }

  bool hasCrystalParams(void) const { return(_icntrl[10] == 1); }

  // The following track CHARMm names (more or less...)
  unsigned int nsteps(void) const { return((unsigned int)_icntrl[3]); }
  int nsavc(void) const { return(_icntrl[2]); }
  int nfile(void) const { return(_icntrl[0]); }
  int nfixed(void) const { return(_icntrl[8]); }

  //! Auto-interleave the coords into a vector of GCoord()'s.
  /*!  This can be a pretty slow operation, so be careful. */
  vector<GCoord> coords(void);

  //! Interlieve coords, selecting entries indexed by map
  vector<GCoord> mappedCoords(const vector<int>& map);


  //! Update an AtomicGroup coordinates with the currently-read frame.
  /*! This assumes that that atomid's of the
    AtomicGroup are indices into the DCD frame and are indexed +1,
    i.e. atomid 7 refers to DCD coords at index 6... */
  void updateGroupCoords(AtomicGroup& g);

private:
  void allocateSpace(const int n);
  void readCrystalParams(void);
  void readCoordLine(vector<float>& v);

  // For reading F77 I/O
  unsigned int readRecordLen(ifstream* ifs);
  DataOverlay* readF77Line(ifstream* ifs, unsigned int *len);



private:
  bool del_stream;          // Must delete the stream pointer...
  ifstream* _ifs;           // Cached file pointer
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
