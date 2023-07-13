
#include <loos.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>    // cif::Document -> Structure
#include <gemmi/gz.hpp> 

namespace cif = gemmi::cif;


int main(int argc, char *argv[]) {
    std::string filename = std::string(argv[1]);
    
    //cif::Document doc = cif::read(gemmi::MaybeGzipped(filename));
    auto structure = gemmi::read_structure_file(filename);
}