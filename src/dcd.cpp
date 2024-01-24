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
  uint DCD::nframes(void) const { return(_nframes); }

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
  // Returns 0 at EOF

  unsigned int DCD::readRecordLen(void) {
    DataOverlay o;
  
    ifs->read(o.c, 4);

    if (ifs->eof())
      return(0);
    if (ifs->fail())
      throw(FileReadError(_filename, "Unable to read DCD record length"));

    uint data = o.ui;
    if (swabbing)
      data = swab(data);

    return(data);
  }


  // Check for endian-ness
  void DCD::endianMatch(pStream& fsw) {
    unsigned long curpos = fsw->tellg();
    unsigned int datum;
    fsw->read((char *)(&datum), sizeof(datum));
    if (fsw->eof() || fsw->fail())
      throw(FileReadError(_filename, "Unable to read first datum from DCD file"));

    fsw->seekg(curpos);


    if (datum != 0x54) {
      datum = swab(datum);
      if (datum != 0x54)
        throw(FileReadError(_filename, "Unable to determine endian-ness of DCD file"));
      swabbing = true;
    } else
      swabbing = false;
  }


  // Read a full line of F77-formatted data.
  // Returns a pointer to the read data and puts the # of bytes read into *len
  // Returns a null pointer and 0 length at EOF
  // Note:  It is up to the caller to swab individual elements...

  DCD::DataOverlay* DCD::readF77Line(unsigned int *len) {
    DataOverlay* ptr;
    unsigned int n, n2;

    *len = 0;
    n = readRecordLen();
    if (n == 0)
      return(0);
    
    ptr = new DataOverlay[n];

    ifs->read((char *)ptr, n);
    if (ifs->fail())
      throw(FileReadError(_filename, "Error reading data record from DCD"));

    n2 = readRecordLen();
    if (n != n2)
      throw(FileReadError(_filename, "Mismatch in record length while reading from DCD"));

    *len = n;
    return(ptr);
  }


  // Determine number of frames in trajectory assuming that all frames have the
  // save size.  If the size is not an integral multiple of the frame size, then
  // punt and return 0.
  //
  // NOTE: Does not preserve current file position
  uint DCD::calculateNumberOfFrames() {
    ifs->seekg(0, std::ios::end);
    std::streampos endpos = ifs->tellg();
    long datasize = endpos - first_frame_pos;
    uint n = datasize / frame_size;
    if (datasize % frame_size)
      return(0);
    return(n);
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
      throw(FileReadError(_filename, "Error while reading DCD header"));

    // Check for the magic name.  Ignore swabbing for now...
    DataOverlay o = ptr[0];
    if (!(o.c[0] == 'C' && o.c[1] == 'O' && o.c[2] == 'R' && o.c[3] == 'D'))
      throw(FileReadError(_filename, "DCD is missing CORD magic marker"));

    // Copy in the ICNTRL data...
    for (i=0; i<20; i++) {
      if (swabbing)
        _icntrl[i] = swab(ptr[i+1].i);
      else
        _icntrl[i] = ptr[i+1].i;
    }

    _nframes = _icntrl[0];
    
    // Extract the delta value...
    if (swabbing)
      _delta = swab(ptr[10].f);
    else
      _delta = ptr[10].f;

    if (nfixed() != 0)
      throw(LOOSError("Fixed atoms not yet supported by LOOS DCD reader"));

    delete[] ptr;

    // Now read in the TITLE info...

    ptr = readF77Line(&len);
    if (!ptr)
      throw(FileReadError(_filename, "Unexpected EOF reading DCD TITLE"));
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
      throw(FileReadError(_filename, "Error reading number of atoms from DCD"));
    if (swabbing)
      _natoms = swab(ptr->i);
    else
      _natoms = ptr->i;
    delete[] ptr;


    // Finally, set internal variables and allocate space for a frame...
    first_frame_pos = ifs->tellg();

    frame_size = 12 * (2 + _natoms);
    if (hasCrystalParams())
      frame_size += 56;

    allocateSpace(_natoms);

    // See if the header is missing frame # information...
    if (_nframes == 0) {
      _nframes = calculateNumberOfFrames();
      if (_nframes == 0 && !suppress_warnings)
        std::cerr << "Warning- DCD '" << _filename << "' appears empty and could not determine size; verify with dcdinfo and/or fix with fixdcd" << std::endl;
      
      ifs->seekg(first_frame_pos);
    }
    
  }



  // Read in and reorder the crystal parameters...
  // NOTE: This is already double!


  bool DCD::readCrystalParams(void) {
    unsigned int len;
    DataOverlay* o;

    o = readF77Line(&len);
    if (!o)
      return(false);

    if (len != 48)
      throw(FileReadError(_filename, "Cannot read crystal parameters"));

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

    return(true);
  }



  // Read a line of coordinates into the specified vector.

  bool DCD::readCoordLine(std::vector<dcd_real>& v) {
    DataOverlay *op;
    int n = _natoms * sizeof(dcd_real);
    unsigned int len;


    op = readF77Line(&len);
    if (!op)
      return(false);
    
    if (len != (unsigned int)n)
      throw(FileReadError(_filename, "Size of coords stored in frame does not match model size"));

    // Recast the values as floats and store them...
    uint i;
    for (i=0; i<_natoms; i++)
      if (swabbing)
        v[i] = swab(op[i].f);
      else
        v[i] = op[i].f;

    delete[] op;

    return(true);
  }


  void DCD::seekFrameImpl(const uint i) {
  
    if (first_frame_pos == 0)
      throw(FileError(_filename, "Trying to seek to a DCD frame without having first read the header"));

    if (i >= nframes())
      throw(FileError(_filename, "Requested DCD frame is out of range"));

    ifs->clear();
    ifs->seekg(first_frame_pos + i * frame_size);
    if (ifs->fail() || ifs->bad())
      throw(FileError(_filename, "Cannot seek to requested frame"));
  }


  // Read in a frame of data.
  // Returns TRUE if success or FALSE if at EOF.
  // Throws an exception if there was an error...  (Should we throw EOF
  // instead?) 

  bool DCD::parseFrame(void)  {

    if (first_frame_pos == 0)
      throw(FileReadError(_filename, "Trying to read a DCD frame without first having read the header."));

    // This will not catch most cases of reading to the end of the file...
    if (ifs->eof())
      return(false);

    if (hasCrystalParams())
      if (!readCrystalParams())
	return(false);

    if (!readCoordLine(xcrds))
      return(false);
    
    if (!readCoordLine(ycrds))
      throw(FileReadError(_filename, "Unexpected EOF reading Y-coordinates from DCD"));
    if (!readCoordLine(zcrds))
      throw(FileReadError(_filename, "Unexepcted EOF reading Z-coordinates from DCD"));
  
    return(true);
  }


  void DCD::rewindImpl(void) {
    ifs->clear();
    ifs->seekg(first_frame_pos);
    if (ifs->fail() || ifs->bad())
      throw(FileError(_filename, "Cannot rewind DCD trajectory"));
  }


  // ----------------------------------------------------------


  std::vector<GCoord> DCD::coords(void) const {
    std::vector<GCoord> crds(_natoms);

    for (uint i=0; i<_natoms; i++) {
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


  void DCD::updateGroupCoordsImpl(AtomicGroup& g) {
    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx >= _natoms)
        throw(TrajectoryError("updating group coords", _filename, "Atom index into trajectory frame is out of bounds"));
      (*i)->coords(GCoord(xcrds[idx], ycrds[idx], zcrds[idx]));
    }

    // Handle periodic boundary conditions (if present)
    if (hasPeriodicBox()) {
      g.periodicBox(periodicBox());
    }
  }



  void DCD::initTrajectory() {
        readHeader();
        bool b = parseFrame();
        if (!b)
            throw(TrajectoryError("reading first frame of DCD during initialization"));
        cached_first = true;
    }
    

}
