/*
  dcd.h
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Notes:

    o Does NOT support fixed atoms

    o Does NOT support velocity format

    o Reorders the crystal parameters (if present) so they are in
      a more sensible order (i.e. a, b, c, alpha, beta, gamma)
      
    o Everything returned is a copy (except for the FILE pointer)

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


class DCD : public boost::noncopyable {

  // Use a union to convert data to appropriate type...
  typedef union { unsigned int ui; int i; char c[4]; float f; } DataOverlay; 

public:

  // Error classes that we may or may not return...
  struct RecordError : public exception { virtual const char *what() const throw() { return("Error while reading F77 record"); }; } record_error;
  struct HeaderError : public exception { virtual const char *what() const throw() { return("Error while reading DCD header"); }; } header_error;
  struct LineError : public exception { virtual const char *what() const throw() { return("Error while reading F77 data line"); }; } line_error;
  struct EndOfFile : public exception { virtual const char *what() const throw() { return("Unexpected end of file"); }; } end_of_file;
  struct GeneralError : public exception {
    GeneralError(const char *s) : msg(s) { };
    virtual const char *what() const throw() { return(msg); };
    const char *msg;
  };


  DCD() : del_stream(false), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { };
  DCD(const string s) :  del_stream(true), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) {
    _ifs = new ifstream(s.c_str());
    if (!_ifs)
      throw(runtime_error("Cannot open DCD file " + s));
    readHeader();
  }

  DCD(const char* s) :  del_stream(true), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) {
    _ifs = new ifstream(s);
    if (!_ifs)
      throw(runtime_error("Cannot open DCD file " + string(s)));
    readHeader();
  }

  DCD(ifstream& ifs) : del_stream(false), _ifs(&ifs), _natoms(0), qcrys(vector<double>(6)), frame_size(0), first_frame_pos(0) { readHeader(); };

  ~DCD() {
    if (del_stream)
      delete _ifs;
  }


  void readHeader(void);
  void readHeader(ifstream& ifs);

  bool readFrame(void);
  bool readFrame(const unsigned int i);

  // Accessor methods...

  int natoms(void) const { return(_natoms); }

  vector<string> titles(void) const { return(_titles); }

  int icntrl(const int i) const { assert(i>=0 && i<20); return(_icntrl[i]); }
  void icntrl(const int i, const int val) { assert(i>=0 && i<20); _icntrl[i] = val; }

  vector<double> crystalParams(void) const { return(qcrys); }

  float delta(void) const { return(_delta); }

  // Return the raw coords...
  vector<dcd_real> xcoords(void) const { return(xcrds); }
  vector<dcd_real> ycoords(void) const { return(ycrds); }
  vector<dcd_real> zcoords(void) const { return(zcrds); }

  bool hasCrystalParams(void) const { return(_icntrl[10] == 1); }

  // The following track CHARMm names (more or less...)
  unsigned int nsteps(void) const { return((unsigned int)_icntrl[3]); }
  int nsavc(void) const { return(_icntrl[2]); }
  int nfile(void) const { return(_icntrl[0]); }
  int nfixed(void) const { return(_icntrl[8]); }

  // Some nice support routines...
  vector<GCoord> coords(void);
  vector<GCoord> mappedCoords(const vector<int>& map);

  // Update the passed AtomicGroup's coords with coords from the
  // currently read frame...  This assumes that that atomid's of the
  // AtomicGroup are indices into the DCD frame and are indexed +1,
  // i.e. atomid 7 refers to DCD coords at index 6...

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
