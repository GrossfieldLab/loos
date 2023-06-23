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
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

int main(int argc, char *argv[]) {

    std::string filename = argv[1];

    HighFive::File file(filename, HighFive::File::ReadOnly);
    auto dataset = file.getDataSet("topology");
    std::vector<std::size_t> dims = dataset.getSpace().getDimensions();
    std::size_t n = dataset.getElementCount();
    std::cout << "n = " << n << std::endl;
    std::cout << "storage size = " << dataset.getStorageSize() << std::endl;

    auto dataspace = dataset.getSpace();
    std::cout << "elements: " << dataspace.getElementCount() << std::endl;

    auto datatype = dataset.getDataType();
    std::cout << "datatype: " << datatype.string() << std::endl;
    std::cout << "fixed len:" << datatype.isFixedLenStr() << std::endl;

    std::cout << H5Easy::getShape(file, "topology")[0] << std::endl;
    std::cout << H5Easy::getShape(file, "topology")[1] << std::endl;
    HighFive::FixedLenStringArray<200000000> str;

    std::cerr << "got here" << std::endl;
    dataset.read(str);
    std::cerr << "got here" << std::endl;
    std::cout << str.getString(0) << std::endl;

    for (auto i : dims) {
        std::cout << i << "\t" << dims[i] << std::endl;
    } 


}