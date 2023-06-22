/*
  test_xml

  code to play with reading openmm xml files

*/


/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2023 Alan Grossfield
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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/detail/file_parser_error.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {

    fstream xml_file("/home/alan/Downloads/system.xml");
    boost::property_tree::ptree system_tree;

    try {
        read_xml(xml_file, system_tree);
    } catch (boost::property_tree::xml_parser_error &e) {
        std :: cout << "Failed to parse the xml string." << e.what();
    } catch (...) {
        std :: cout << "Failed !!!";
    }

    for (auto& p : system_tree.get_child("System.Particles")) {
        std :: cout << "[" << p.first << "]" << std :: endl;
        for (auto& c : p.second) {
            std :: cout << "Tag : [" << c.first.data() << "], ";
            double mass = c.second.get<double>("mass");
            std :: cout << "Value : [" << mass << "]" << std :: endl;
        }   
    }

return 0;
}

