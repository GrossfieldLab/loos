/*
  Testing code for HDF5 support 

  (c) 2023 Alan Grossfield
           Department of Biochemistry and Biophysics
           University of Rochester School of Medicine and Dentistry


*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2011 Tod D. Romo
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
#include <H5Cpp.h>
#include <boost/json.hpp>
//using namespace boost::json;

std::string getTopology(H5::H5File &file) {
  H5::DataSet dataset = file.openDataSet("topology");
  H5::DataSpace dataspace = dataset.getSpace();
  H5::DataType datatype = dataset.getDataType();
  H5std_string topology;
  dataset.read(topology, datatype, dataspace, dataspace);
  return std::string(topology);
}


int main(int argc, char *argv[]) {

    std::string filename = argv[1];

    // Turn off the auto-printing when failure occurs so that we can
  	// handle the errors appropriately
    H5::Exception::dontPrint();

    H5::H5File file(filename, H5F_ACC_RDONLY);
    std::string topology_json = getTopology(file);
    boost::json::value topology = boost::json::parse(topology_json);

    loos::AtomicGroup ag;
    uint index = 0;

    boost::json::array chains = topology.at("chains").as_array();
    for (auto chain: chains) {
      boost::json::array residues = chain.at("residues").as_array();
      for (auto residue: residues) {
        int resnum = residue.at("resSeq").as_int64();
        std::string resname = residue.at("name").as_string().data();

        boost::json::array atoms = residue.at("atoms").as_array();
        for (auto atom: atoms) {
          std::string atom_name = atom.at("name").as_string().data();
          int id = atom.at("index").as_int64();
          id++; // HDF5 is 0 based, LOOS is 1 based
          std::string element = atom.at("element").as_string().data();

          loos::pAtom pa(new loos::Atom);
          pa->name(atom_name);
          pa->id(id);
          pa->index(index);
          pa->resid(resnum);
          pa->resname(resname);
          pa->PDBelement(element);

          ag.append(pa);

          index++;
        }
      }
    }

    // Add the bonds
    // We're assuming the atoms are in order
    boost::json::array bonds = topology.at("bonds").as_array();
    for (auto bond: bonds) {
      int atom1 = bond.at(0).as_int64();
      int atom2 = bond.at(1).as_int64();
      ag[atom1]->addBond(ag[atom2]);
      ag[atom2]->addBond(ag[atom1]);
    }

    // TODO: I need an example that has constraints
    // If there are constraints, we should treat them like bonds

    // Read the coordinates and box size
    H5::DataSet box_dataset = file.openDataSet("cell_lengths");
    H5::DataSpace box_dataspace = box_dataset.getSpace();
    H5::DataType box_datatype = box_dataset.getDataType();
    hsize_t box_dims[3];
    int ndims = box_dataspace.getSimpleExtentDims(box_dims, NULL);
    std::cout << ndims << std::endl;
    std::cout << box_dims[0] << std::endl;
    std::cout << box_dims[1] << std::endl;

    int num_elements = box_dims[0]* box_dims[1];
    auto box_lengths = new float[num_elements];

    // TODO: this reads the whole traj at once, but I should learn how to read just
    //       one frame at a time to do the traj right
    //box_dataset.read(box_lengths, box_datatype, box_dataspace, box_dataspace);

    // Read the nth frame
    hsize_t n = 32;
    hsize_t offset[2] = {n, 0};
    hsize_t count[2] = {1, box_dims[1]};
    float one_box[box_dims[1]];
    hsize_t offset_out[1] = {0};
    hsize_t count_out[1] = {box_dims[1]};
    H5::DataSpace memspace(1, count_out);

    box_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
    box_dataset.read(one_box, box_datatype, memspace, box_dataspace);


    // Set the periodic box to the first frame of the traj and fix the units
    loos::GCoord box(one_box[0], one_box[1], one_box[2]);
    box *= 10.0; // convert to Angstroms
    ag.periodicBox(box);

    std::cout << box << std::endl;
    //std::cout << loos::PDB::fromAtomicGroup(ag) << std::endl;
  
    return 0;
}