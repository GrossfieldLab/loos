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



#include <mdtrajtraj.hpp>
#include <AtomicGroup.hpp>
#include <iomanip>
#include <sstream>

namespace loos {

  // Scan the trajectory file to determine frame sizes and box
  void MDTrajTraj::init(void) {
    file.openFile(_filename, H5F_ACC_RDONLY);

    hsize_t box_dims[3];
    if (H5Lexists(file.getId(), "cell_lengths", H5P_DEFAULT)) {
      periodic = true;
      box_dataset = file.openDataSet("cell_lengths");
      box_dataspace = box_dataset.getSpace();
      box_datatype = box_dataset.getDataType();
      int box_ndims = box_dataspace.getSimpleExtentDims(box_dims, NULL);
    } else {
      periodic = false;
    }

    coords_dataset = file.openDataSet("coordinates");
    coords_dataspace = coords_dataset.getSpace();
    coords_datatype = coords_dataset.getDataType();
    hsize_t coords_dims[3];
    int coords_ndims = coords_dataspace.getSimpleExtentDims(coords_dims, NULL);

    if (periodic) {
      if (coords_dims[0] != box_dims[0]) {
        throw(FileError(_filename, "Number of frames in box and coords datasets in HDF5 do not match"));
      }
    }

    _nframes = coords_dims[0];
    if (_natoms != coords_dims[1]) {
      throw(FileError(_filename, "Number of atoms in HDF5 does not match the AtomicGroup"));
    }

    // Allocate space to store the coordinates
    frame.resize(_natoms);
    one_frame = new float[_natoms][3];

    // Now cache the first frame...
		readRawFrame(0);
		cached_first = true;
  }

  bool MDTrajTraj::parseFrame(void) {
    if (atEnd()) {
      return(false);
    }
    readRawFrame(_current_frame);
    return(true);
  }

  void MDTrajTraj::readRawFrame(const uint i) {
    
    // Read the periodic box
    if (periodic) {
      hsize_t offset[2] = {i, 0};
      hsize_t count[2] = {1, 3};
      hsize_t count_out[1] = {3};
      float one_box[3];
      H5::DataSpace memspace(1, count_out);        
      box_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
      box_dataset.read(one_box, box_datatype, memspace, box_dataspace);
      // copy into periodic box and convert from nm to Angstroms
      for (int j = 0; j < 3; ++j) {
        box[j] = 10.0*one_box[j];
      }
    }

    // Read the coordinates

    hsize_t offset_coord[3] = {i, 0, 0};
    hsize_t count_coord[3] = {1, _natoms, 3};
    hsize_t count_coord_out[2] = {_natoms, 3};
    H5::DataSpace memspace_coord(2, count_coord_out);
    coords_dataspace.selectHyperslab(H5S_SELECT_SET, count_coord, offset_coord);
    coords_dataset.read(one_frame, coords_datatype, memspace_coord, coords_dataspace);
    
    // copy coords into frame and convert from nm to Angstroms
    for (int j=0; j < _natoms; ++j) {
      for (int k=0; k < 3; ++k) {
        frame[j][k] = 10.0*one_frame[j][k];
      }
    }
  }

  void MDTrajTraj::seekFrameImpl(const uint i) {

    /*
    if (i >= _nframes)
      throw(FileError(_filename, "Attempting to seek frame beyond end of trajectory"));

    _current_frame = i;
    */
  }


  void MDTrajTraj::updateGroupCoordsImpl(AtomicGroup& g) {

    for (AtomicGroup::iterator i = g.begin(); i != g.end(); ++i) {
      uint idx = (*i)->index();
      if (idx >= _natoms)
        throw(LOOSError(**i, "Atom index into trajectory is out of bounds"));
      (*i)->coords(frame[idx]);
    }
    
    if (periodic)
      g.periodicBox(box);
  }
}
