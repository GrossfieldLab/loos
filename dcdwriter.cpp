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




#include <loos.hpp>
#include <dcdwriter.hpp>

void DCDWriter::writeF77Line(StreamWrapper& ofs, const char* const data, const unsigned int len) {
  DataOverlay d;

  d.ui = len;

  ofs()->write((char *)&d, sizeof(len));
  ofs()->write(data, len);
  ofs()->write((char *)&d, sizeof(len));
}


string DCDWriter::fixStringSize(const string& s, const unsigned int n) {
  string result(s);

  if (s.size() < n) {
    int i = n - s.size();
    result += string (i, ' ');
  } else if (s.size() > n)
    result = s.substr(0, n);

  return(result);
}


void DCDWriter::writeHeader(void) {
  unsigned int icntrl[21];
  DataOverlay *dop = (DataOverlay *)(&icntrl[0]);
  unsigned int i;
  for (i=0; i<21; i++)
    icntrl[i] = 0;

  icntrl[1] = _nsteps;
  icntrl[2] = 1;
  icntrl[3] = 1;
  icntrl[4] = _nsteps;
  icntrl[8] = _natoms * 3 - 6;
  icntrl[9] = 0;
  dop[10].f = _timestep;
  icntrl[11] = _has_box;
  icntrl[20] = 27;
  dop[0].c[0] = 'C'; dop[0].c[1] = 'O'; dop[0].c[2] = 'R'; dop[0].c[3] = 'D';

  writeF77Line(_ofs, (char *)dop, 21 * sizeof(unsigned int));

  unsigned int size = 4 + 80 * _titles.size();
  char *ptr = new char[size];
  unsigned int *iptr = (unsigned int *)ptr;
  *iptr = _titles.size();
  for (i=0; i<_titles.size(); i++) {
    string s = fixStringSize(_titles[i], 80);
    memcpy(ptr + sizeof(unsigned int) + 80*i, s.c_str(), 80);
  }
  writeF77Line(_ofs, ptr, size);
  delete[] ptr;
  

  dop[0].ui = _natoms;
  writeF77Line(_ofs, (char *)dop, 1 * sizeof(unsigned int));
}



void DCDWriter::writeBox(const GCoord& box) {
  double xtal[6] = { 1.0, box[0], 1.0, box[1], box[2], 1.0 };

  writeF77Line(_ofs, (char *)xtal, 6*sizeof(double));
}



void DCDWriter::writeFrame(const AtomicGroup& grp) {

  if (_current > _nsteps)
    throw(runtime_error("Attempting to write more frames than requested."));

  if (grp.size() != _natoms)
    throw(runtime_error("Frame group atom count mismatch"));

  if (!_has_box && grp.isPeriodic())
    throw(runtime_error("Frame has periodic info but none was requested to be written to the DCD."));

  if (_has_box)
    writeBox(grp.periodicBox());

  float *data = new float[_natoms];
  int i;
  for (i=0; i<_natoms; i++)
    data[i] = grp[i]->coords().x();
  writeF77Line(_ofs, (char *)data, _natoms * sizeof(float));

  for (i=0; i<_natoms; i++)
    data[i] = grp[i]->coords().y();
  writeF77Line(_ofs, (char *)data, _natoms * sizeof(float));

  for (i=0; i<_natoms; i++)
    data[i] = grp[i]->coords().z();
  writeF77Line(_ofs, (char *)data, _natoms * sizeof(float));

  delete[] data;

  ++_current;
}


void DCDWriter::writeFrame(const vector<AtomicGroup>& grps) {
  vector<AtomicGroup>::const_iterator i;

  for (i= grps.begin(); i != grps.end(); i++)
    writeFrame(*i);
}

