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




#include <iostream>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <dcd.hpp>
#include <AtomicGroup.hpp>


namespace loos {

  uint DCD::natoms(void) const { return(_natoms); }
  bool DCD::hasPeriodicBox(void) const { return(_icntrl[10] == 1); }
  GCoord DCD::periodicBox(void) const { return(GCoord(qcrys[0], qcrys[1], qcrys[2])); }


  bool DCD::suppress_warnings = false;
  
  
  std::vector<std::string> DCD::titles(void) const { return(_titles); }

  int DCD::icntrl(const int i) const { assert(i>=0 && i<20); return(_icntrl[i]); }
  void DCD::icntrl(const int i, const int val) { assert(i>=0 && i<20); _icntrl[i] = val; }

  // * legacy *
  std::vector<double> DCD::crystalParams(void) const { return(qcrys); }
  bool DCD::hasCrystalParams(void) const { return(_icntrl[10] == 1); }

  float DCD::timestep(void) const { return(_delta); }
  uint DCD::nframes(void) const { return(_icntrl[0]); }

  std::vector<dcd_real> DCD::xcoords(void) const { return(xcrds); }
  std::vector<dcd_real> DCD::ycoords(void) const { return(ycrds); }
  std::vector<dcd_real> DCD::zcoords(void) const { return(zcrds); }

  // The following track CHARMm names (more or less...)
  unsigned int DCD::nsteps(void) const { return(_icntrl[3]); }
  float DCD::delta(void) const { return(_delta); }
  int DCD::nsavc(void) const { return(_icntrl[2]); }
  int DCD::nfile(void) const { return(_icntrl[0]); }
  int DCD::nfixed(void) const { return(_icntrl[8]); }

  bool DCD::nativeFormat(void) const { return(!swabbing); }


  // Allocate the input buffer and pre-size the coordinate vectors.
  // n = number of atoms

  void DCD::allocateSpace(const int n) {
    xcrds = std::vector<dcd_real>(n);
    ycrds = std::vector<dcd_real>(n);
    zcrds = std::vector<dcd_real>(n);
  }



  // Read the F77 record length from the file stream

  unsigned int DCD::readRecordLen(void) {
    DataOverlay o;
  
    ifs()->read(o.c, 4);

    if (ifs()->eof())
      throw(end_of_file);
    if (ifs()->fail())
      throw(record_error);

    uint data = o.ui;
    if (swabbing)
      data = swab(data);

    return(data);
  }


  // Check for endian-ness
  void DCD::endianMatch(StreamWrapper& fsw) {
    unsigned long curpos = fsw()->tellg();
    unsigned int datum;
    ifs()->read((char *)(&datum), sizeof(datum));
    fsw()->seekg(curpos);

    if (ifs()->eof() || ifs()->fail())
      throw(GeneralError("Unable to read first datum from DCD file"));

    if (datum != 0x54) {
      datum = swab(datum);
      if (datum != 0x54)
        throw(GeneralError("Unable to determine endian-ness of DCD file"));
      swabbing = true;
    } else
      swabbing = false;
  }


  // Read a full line of F77-formatted data.
  // Returns a pointer to the read data and puts the # of bytes read into *len
  // Note:  It is up to the caller to swab individual elements...

  DCD::DataOverlay* DCD::readF77Line(unsigned int *len) {
    DataOverlay* ptr;
    unsigned int n, n2;

    n = readRecordLen();

    ptr = new DataOverlay[n];

    ifs()->read((char *)ptr, n);
    if (ifs()->fail())
      throw(line_error);

    n2 = readRecordLen();
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

    endianMatch(ifs);
    ptr = readF77Line(&len);
    if (len != 84)
      throw(header_error);

    // Check for the magic name.  Ignore swabbing for now...
    DataOverlay o = ptr[0];
    if (!(o.c[0] == 'C' && o.c[1] == 'O' && o.c[2] == 'R' && o.c[3] == 'D'))
      throw(header_error);

    // Copy in the ICNTRL data...
    for (i=0; i<20; i++) {
      if (swabbing)
        _icntrl[i] = swab(ptr[i+1].i);
      else
        _icntrl[i] = ptr[i+1].i;
    }
      
    // Extract the delta value...
    if (swabbing)
      _delta = swab(ptr[10].f);
    else
      _delta = ptr[10].f;

    if (nfixed() != 0)
      throw(GeneralError("Fixed atoms not yet supported"));

    delete[] ptr;

    // Now read in the TITLE info...

    ptr = readF77Line(&len);
    char sbuff[81];
    int ntitle = ptr[0].i;
    if (swabbing)
      ntitle = swab(ntitle);

    cp = (char *)(ptr + 1);

    for (i=0; i<ntitle; i++) {
      memcpy(sbuff, cp + 80*i, 80);
      sbuff[80] = '\0';
      std::string s(sbuff);
      _titles.push_back(s);
    }
    delete[] ptr;

    // get the NATOMS...
    ptr = readF77Line(&len);
    if (len != 4)
      throw(header_error);
    if (swabbing)
      _natoms = swab(ptr->i);
    else
      _natoms = ptr->i;
    delete[] ptr;


    // Finally, set internal variables and allocate space for a frame...
    first_frame_pos = ifs()->tellg();

    frame_size = 12 * (2 + _natoms);
    if (hasCrystalParams())
      frame_size += 56;

    allocateSpace(_natoms);

    // Issue warnings
    if (nframes() == 0 && !suppress_warnings)
        std::cerr << "Warning- DCD '" << _filename << "' appears empty; verify with dcdinfo and fix with fixdcd" << std::endl;
  }



  // Read in and reorder the crystal parameters...
  // NOTE: This is already double!


  void DCD::readCrystalParams(void) {
    unsigned int len;
    DataOverlay* o;

    o = readF77Line(&len);

    if (len != 48)
      throw(GeneralError("Error while reading crystal parameters"));

    double* dp = reinterpret_cast<double*>(o);
    
    qcrys[0] = dp[0];
    qcrys[1] = dp[2];
    qcrys[2] = dp[5];
    qcrys[3] = dp[1];
    qcrys[4] = dp[3];
    qcrys[5] = dp[4];
    
    if (swabbing)
        for (int i=0; i<6; ++i)
            qcrys[i] = swab(qcrys[i]);

    delete[] o;
  }



  // Read a line of coordinates into the specified vector.

  void DCD::readCoordLine(std::vector<dcd_real>& v) {
    DataOverlay *op;
    int n = _natoms * sizeof(dcd_real);
    unsigned int len;


    op = readF77Line(&len);

    if (len != (unsigned int)n)
      throw(GeneralError("Error while reading coordinates"));

    // Recast the values as floats and store them...
    int i;
    for (i=0; i<_natoms; i++)
      if (swabbing)
        v[i] = swab(op[i].f);
      else
        v[i] = op[i].f;

    delete[] op;

  }


  void DCD::seekFrameImpl(const uint i) {
  
    if (first_frame_pos == 0)
      throw(GeneralError("Trying to seek to a DCD frame without having first read the header"));

    if (i >= nframes())
      throw(GeneralError("Requested DCD frame is out of range"));

    ifs()->clear();
    ifs()->seekg(first_frame_pos + i * frame_size);
    if (ifs()->fail() || ifs()->bad()) {
      std::ostringstream s;
      s << "Cannot seek to frame " << i;
      throw(GeneralError(s.str().c_str()));
    }
  }


  // Read in a frame of data.
  // Returns TRUE if success or FALSE if at EOF.
  // Throws an exception if there was an error...  (Should we throw EOF
  // instead?) 

  bool DCD::parseFrame(void) {

    if (first_frame_pos == 0)
      throw(GeneralError("Trying to read a DCD frame without first having read the header."));

    if (ifs()->eof())
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


  void DCD::rewindImpl(void) {
    ifs()->clear();
    ifs()->seekg(first_frame_pos);
    if (ifs()->fail() || ifs()->bad())
      throw(GeneralError("Error rewinding file"));
  }


  // ----------------------------------------------------------


  std::vector<GCoord> DCD::coords(void) {
    std::vector<GCoord> crds(_natoms);
    int i;

    for (i=0; i<_natoms; i++) {
      crds[i].x(xcrds[i]);
      crds[i].y(ycrds[i]);
      crds[i].z(zcrds[i]);
    }

    return(crds);
  }

  std::vector<GCoord> DCD::mappedCoords(const std::vector<int>& indices) {
    std::vector<int>::const_iterator iter;
    std::vector<GCoord> crds(indices.size());

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
    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      int idx = (*i)->id()-1;
      if (idx < 0 || idx >= _natoms)
        throw(LOOSError(**i, "Atom index into the trajectory frame is out of bounds"));
      (*i)->coords(GCoord(xcrds[idx], ycrds[idx], zcrds[idx]));
    }

    // Handle periodic boundary conditions (if present)
    if (hasPeriodicBox()) {
      g.periodicBox(periodicBox());
    }
  }



    void DCD::initTrajectory() 
    {
        readHeader();
        bool b = parseFrame();
        if (!b)
            throw(GeneralError("Cannot read first frame of DCD during initialization"));
        cached_first = true;
    }
    

}
