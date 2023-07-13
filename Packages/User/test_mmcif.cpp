
#include <loos.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>    // cif::Document -> Structure
#include <gemmi/gz.hpp> 

namespace cif = gemmi::cif;


int main(int argc, char *argv[]) {
    std::string filename = std::string(argv[1]);
    
    auto structure = gemmi::read_structure_file(filename, gemmi::CoorFormat::Mmcif);
}