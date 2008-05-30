/*
  dcd.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  DCD I/O

*/


#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <tr1/memory>

#include <stdio.h>
#include <string.h>
#include <assert.h>

using namespace std;


#include <dcd.hpp>


// Allocate the input buffer and pre-size the coordinate vectors.
// n = number of atoms

void DCD::allocateSpace(const int n) {
  xcrds = vector<dcd_real>(n);
  ycrds = vector<dcd_real>(n);
  zcrds = vector<dcd_real>(n);
}



// Read the F77 record length from the file stream

unsigned int DCD::readRecordLen(ifstream* ifs) {
  DataOverlay o;
  
  ifs->read(o.c, 4);

  if (ifs->eof())
    throw(end_of_file);
  if (ifs->fail())
    throw(record_error);


  return(o.ui);
}


// Read a full line of F77-formatted data.
// Returns a pointer to the read data and puts the # of bytes read into *len

DCD::DataOverlay* DCD::readF77Line(ifstream* ifs, unsigned int *len) {
  DataOverlay* ptr;
  unsigned int n, n2;

  n = readRecordLen(ifs);

  ptr = new DataOverlay[n];

  ifs->read((char *)ptr, n);
  if (ifs->fail())
    throw(line_error);

  n2 = readRecordLen(ifs);
  if (n != n2)
    throw(line_error);

  *len = n;
  return(ptr);
}



// Read in the DCD header...

void DCD::readHeader(void) {
  DataOverlay *ptr;
  unsigned int len;
  char *cp;
  int i;

  assert(_ifs != 0);

  ptr = readF77Line(_ifs, &len);
  if (len != 84)
    throw(header_error);

  // Check for the magic name.  Ignore swabbing for now...
  if (!(ptr->c[0] == 'C' && ptr->c[1] == 'O' && ptr->c[2] == 'R' && ptr->c[3] == 'D'))
    throw(header_error);

  // Copy in the ICNTRL data...
  for (i=0; i<20; i++)
    _icntrl[i] = (ptr[i+1]).i;

  // Extract the delta value...
  _delta = ptr[10].f;

  if (nfixed() != 0)
    throw(GeneralError("Fixed atoms not yet supported"));

  delete[] ptr;

  // Now read in the TITLE info...

  ptr = readF77Line(_ifs, &len);
  char sbuff[81];
  int ntitle = ptr[0].i;
  cp = (char *)(ptr + 1);

  for (i=0; i<ntitle; i++) {
    memcpy(sbuff, cp + 80*i, 80);
    sbuff[80] = '\0';
    string s(sbuff);
    _titles.push_back(s);
  }
  delete[] ptr;

  // get the NATOMS...
  ptr = readF77Line(_ifs, &len);
  if (len != 4)
    throw(header_error);
  _natoms = ptr->i;
  delete[] ptr;


  // Finally, set internal variables and allocate space for a frame...
  first_frame_pos = _ifs->tellg();
  frame_size = 12 * (2 + _natoms);
  if (hasCrystalParams())
    frame_size += 56;

  allocateSpace(_natoms);
}



void DCD::readHeader(ifstream& ifs) {
  if (_ifs && del_stream)
    delete(_ifs);

  _ifs = &ifs;
  del_stream = false;
  readHeader();
}



// Read in and reorder the crystal parameters...
// NOTE: This is already double!


void DCD::readCrystalParams(void) {
  unsigned int len;
  DataOverlay* o;

  o = readF77Line(_ifs, &len);

  if (len != 48)
    throw(GeneralError("Error while reading crystal parameters"));
  double *dp = (double *)o;  // The values are actually doubles, so we
			     // have to use a little trickery... 

  qcrys[0] = dp[0];
  qcrys[1] = dp[2];
  qcrys[2] = dp[5];
  qcrys[3] = dp[1];
  qcrys[4] = dp[3];
  qcrys[5] = dp[4];

  delete[] o;
}



// Read a line of coordinates into the specified vector.

void DCD::readCoordLine(vector<dcd_real>& v) {
  DataOverlay *op;
  int n = _natoms * sizeof(dcd_real);
  unsigned int len;


  op = readF77Line(_ifs, &len);

  if (len != (unsigned int)n)
    throw(GeneralError("Error while reading coordinates"));

  // Recast the values as floats and store them...
  int i;
  for (i=0; i<_natoms; i++)
    v[i] = op[i].f;

  delete[] op;

}


// Read in a frame of data.
// Returns TRUE if success or FALSE if at EOF.
// Throws an exception if there was an error...  (Should we throw EOF
// instead?) 

bool DCD::readFrame(void) {
  if (_ifs->eof())
    return(false);

  try {
    if (hasCrystalParams())
      readCrystalParams();
  }
  catch (EndOfFile& e) { return(false); }
  catch (...) { throw; }

  // Only check for EOF for the first line (in case we didn't try to
  // read crystal params first...
  try { readCoordLine(xcrds); }
  catch (EndOfFile& e) { return(false); }
  catch (...) { throw; }

  readCoordLine(ycrds);
  readCoordLine(zcrds);
  
  return(true);
}


// Read in a specified DCD frame...

bool DCD::readFrame(const unsigned int i) {

  if (first_frame_pos == 0)
    throw(GeneralError("Trying to read a DCD frame without having first read the header"));

  if (i >= nsteps())
    throw(GeneralError("Requested DCD frame is out of range"));

  _ifs->seekg(first_frame_pos + i * frame_size);
  if (_ifs->fail() || _ifs->bad())
    throw(GeneralError("Cannot seek to frame"));
  
  return(readFrame());
}



// ----------------------------------------------------------


vector<GCoord> DCD::coords(void) {
  vector<GCoord> crds(_natoms);
  int i;

  for (i=0; i<_natoms; i++) {
    crds[i].x(xcrds[i]);
    crds[i].y(ycrds[i]);
    crds[i].z(zcrds[i]);
  }

  return(crds);
}

vector<GCoord> DCD::mappedCoords(const vector<int>& indices) {
  vector<int>::const_iterator iter;
  vector<GCoord> crds(indices.size());

  int j = 0;
  for (iter = indices.begin(); iter != indices.end(); iter++, j++) {
    int index = *iter;
    crds[j].x(xcrds[index]);
    crds[j].y(ycrds[index]);
    crds[j].z(zcrds[index]);
  }

  return(crds);
}


void DCD::updateGroupCoords(AtomicGroup& g) {
  AtomicGroup::Iterator iter(g);
  pAtom pa;

  while(pa = iter()) {
    int i = pa->id() - 1;
    if (i < 0 || i >= _natoms)
      throw(runtime_error("Attempting to index a nonexistent atom in DCD::updateStructure()"));
    GCoord c(xcrds[i], ycrds[i], zcrds[i]);
    pa->coords(c);
  }

  // Handle periodic boundary conditions (if present)
  if (hasCrystalParams()) {
    GCoord c(qcrys[0], qcrys[1], qcrys[2]);
    g.periodicBox(c);
  }
}
