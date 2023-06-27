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


    std::cout << loos::PDB::fromAtomicGroup(ag) << std::endl;
  
    return 0;
}