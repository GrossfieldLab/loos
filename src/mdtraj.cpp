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

#include <mdtraj.hpp>
#include <utils.hpp>
#include <boost/algorithm/string.hpp>


namespace loos {


  MDTraj* MDTraj::clone(void) const {
    return(new MDTraj(*this));
  }

  MDTraj MDTraj::copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    MDTraj p(grp);
    
    return(p);
  }

  void MDTraj::read() {
    // turn off printed exceptions
    //H5::Exception::dontPrint();
    // Open the hdf5 file
    H5::H5File file(_filename, H5F_ACC_RDONLY);

    // Retrieve the topology from the HDF5 file
    std::string topology_json = getTopology(file);

    // Parse the json into a property tree
    boost::json::value topology = boost::json::parse(topology_json);
    
    // get the atoms
    topologyToAtoms(topology);

    // get the bonds
    topologyToBonds(topology);

    // TODO: read the first frame of coordinates from the HDF5 file

  }

  std::string MDTraj::getTopology(H5::H5File &file) {
    H5::DataSet dataset = file.openDataSet("topology");
    H5::DataSpace dataspace = dataset.getSpace();
    H5::DataType datatype = dataset.getDataType();
    H5std_string topology;
    dataset.read(topology, datatype, dataspace, dataspace);
    return std::string(topology);   
  }     

  void MDTraj::topologyToAtoms(const boost::json::value& topology) {
    
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
          id++; // HDF5 is 0-based, most files expect id to be 1-based
          std::string element = atom.at("element").as_string().data();

          pAtom pa(new loos::Atom);
          pa->name(atom_name);
          pa->id(id);
          pa->index(index);
          pa->resid(resnum);
          pa->resname(resname);
          pa->PDBelement(element);

          append(pa);

          _max_index = index;
          index++;
        }
      }
    }
  }

void MDTraj::topologyToBonds(const boost::json::value& topology) {
  // We're assuming the atoms are in order
  boost::json::array bonds = topology.at("bonds").as_array();
  for (auto bond: bonds) {
    int atom1 = bond.at(0).as_int64();
    int atom2 = bond.at(1).as_int64();
    atoms[atom1]->addBond(atoms[atom2]);
    atoms[atom2]->addBond(atoms[atom1]);
  } 
}

}
