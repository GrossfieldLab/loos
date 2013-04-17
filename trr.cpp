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


#include <trr.hpp>


namespace loos {

  // GROMACS magic constants...
  const int TRR::magic = 1993;
  const int TRR::DIM = 3;


  // GROMACS voodoo to determine precision of stored data...
  int TRR::floatSize(Header& hdr) {
    int i;

    if (hdr.box_size)
      i = hdr.box_size/(DIM*DIM);
    else if (hdr.x_size)
      i = hdr.x_size/(hdr.natoms * DIM);
    else if (hdr.v_size)
      i = hdr.v_size/(hdr.natoms * DIM);
    else if (hdr.f_size)
      i = hdr.f_size/(hdr.natoms * DIM);
    else
      throw(std::runtime_error("Cannot determine float size"));

    if (i != sizeof(float) && i != sizeof(double))
      throw(std::runtime_error("Float size does not match native sizes"));

    return(i);
  }


  bool TRR::readHeader(Header& hdr) {
    int magic_no;

    xdr_file.read(magic_no);
    if ((xdr_file.get())->eof())
      return(false);

    if (magic_no != magic)
      throw(std::runtime_error("Invalid TRR magic number"));

    std::string version;
    xdr_file.read(version);
    xdr_file.read(hdr.ir_size);
    xdr_file.read(hdr.e_size);
    xdr_file.read(hdr.box_size);
    xdr_file.read(hdr.vir_size);
    xdr_file.read(hdr.pres_size);
    xdr_file.read(hdr.top_size);
    xdr_file.read(hdr.sym_size);
    xdr_file.read(hdr.x_size);
    xdr_file.read(hdr.v_size);
    xdr_file.read(hdr.f_size);
    xdr_file.read(hdr.natoms);
    
    int flsz = floatSize(hdr);
    hdr.bDouble = (flsz == sizeof(double));

    xdr_file.read(hdr.step);
    xdr_file.read(hdr.nre);

    if (hdr.bDouble) {
      xdr_file.read(hdr.td);
      hdr.tf = hdr.td;
      xdr_file.read(hdr.lambdad);
      hdr.lambdaf = hdr.lambdad;
    } else {
      xdr_file.read(hdr.tf);
      hdr.td = hdr.tf;
      xdr_file.read(hdr.lambdaf);
      hdr.lambdad = hdr.lambdaf;
    }

    if ((xdr_file.get())->fail())
      throw(std::runtime_error("Error reading TRR header"));

    return(true);
  }


  // Initialize the object, along with scanning file for frames to
  // build the frame index, and finally caches the first frame.
  void TRR::init(void) {
    Header h;
    h.natoms = 0;

    // Since the trajectory can have different number of atoms per
    // frame, we track the max and then reserve that space...
    int maxatoms = 0;

    rewindImpl();
    frame_indices.clear();
    
    size_t frame_start = (xdr_file.get())->tellg();
    while (readHeader(h)) {
      frame_indices.push_back(frame_start);
      if (h.natoms > maxatoms)
        maxatoms = h.natoms;

      uint b = sizeof(internal::XDRReader::block_type);

      size_t offset = (h.box_size ? DIM*DIM*b : 0) + (h.vir_size ? DIM*DIM*b : 0)
        + (h.pres_size ? DIM*DIM*b : 0) + (h.x_size ? h.natoms*DIM*b : 0)
        + (h.v_size ? h.natoms*DIM*b : 0) + (h.f_size ? h.natoms*DIM*b : 0);

      (xdr_file.get())->seekg(offset, std::ios_base::cur);
      frame_start = (xdr_file.get())->tellg();
    }

    coords_.reserve(maxatoms);
    velo_.reserve(maxatoms);
    forc_.reserve(maxatoms);

    rewindImpl();
    parseFrame();
    cached_first = true;
    hdr_ = h;
  }

  void TRR::updateGroupCoords(AtomicGroup& g) {

    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      int idx = (*i)->id()-1;
      if (idx < 0 || static_cast<uint>(idx) > natoms())
        throw(std::runtime_error("atom index into trajectory frame is out of range"));
      (*i)->coords(coords_[idx]);
    }

    if (hdr_.box_size)
      g.periodicBox(box);
  }

  void TRR::seekFrameImpl(uint i) {
    if (i >= frame_indices.size())
      throw(std::runtime_error("Trying to seek past the end of the file"));
    
    ifs()->clear();
    ifs()->seekg(frame_indices[i], std::ios_base::beg);
  }


}
