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
    //H5:Exception::dontPrint();

    H5::H5File file(filename, H5F_ACC_RDONLY);
#if 0
    H5::DataSet dataset = file.openDataSet("topology");

    H5::DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    std::cout << rank << std::endl;    
    int length = dataspace.getSimpleExtentNpoints();
    std::cout << length << std::endl;    
    const int ndims = dataspace.getSimpleExtentNdims();
    std::cout << ndims << std::endl;
    hsize_t * dims = new hsize_t[1]; 
    dataspace.getSimpleExtentDims(dims, NULL);
    std::cout << dims[0] << std::endl;    

    H5::DataType datatype = dataset.getDataType();
    H5std_string topology;
    dataset.read(topology, datatype, dataspace, dataspace);
    std::cout << topology << std::endl;


    delete[] dims;
#endif

    std::string topology_json = getTopology(file);
    //std::cout << topology << std::endl;

    // Need to put this into try catch
    boost::json::value topology = boost::json::parse(topology_json);

    loos::AtomicGroup ag;
    uint index = 0;

    boost::json::array chains = topology.at("chains").as_array();
    for (auto chain: chains) {
      boost::json::array residues = chain.at("residues").as_array();
      for (auto residue: residues) {
        uint resnum = residue.at("resSeq").as_uint64();
        std::string resname = std::string_view(residue.at("name").as_string());

        boost::json::array atoms = residue.at("atoms").as_array();
        for (auto atom: atoms) {
          std::stringview atom_name = boost::json::value_to<std::string>(atom.at("name").as_string());
          uint id = atom.at("index").as_uint64();
          std::stringview element = boost::json::value_to<std::string>(atom.at("element").as_string());

          loos::pAtom pa(new loos::Atom);
          pa->name(atom_name);
          pa->id(id);
          pa->index(index);
          pa->resid(resnum);
          pa->resname(resname);
          pa->PDBelement(element);

          //ag._atomid_to_patom[pa->id()] = pa;
          ag.append(pa);

          index++;
        }
      }
    }

  
    return 0;
}