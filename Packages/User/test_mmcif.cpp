
#include <loos.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>    // cif::Document -> Structure
#include <gemmi/gz.hpp> 

namespace cif = gemmi::cif;


int main(int argc, char *argv[]) {
    std::string filename = std::string(argv[1]);
    
    auto structure = gemmi::read_structure_file(filename, gemmi::CoorFormat::Mmcif);
    // TODO: hard-wired to read first model, but there should probably be a way to read others
    auto model = structure.first_model();
    int atom_index = 0;
    int residue_number = 1;
    for (auto chain:model.chains) {
        std::string chain_name = chain.name;
        for (auto residue:chain.residues) {
            std::string residue_name = residue.name;
            std::string label_seq = residue.label_seq.str();
            std::string res_entity_id = residue.entity_id;
            std::cout << "Residue: "
                      << residue_number << "\t"
                      << label_seq << "\t"
                      << res_entity_id << std::endl;
            for (auto atom:residue.atoms) {
                std::cout << atom.name << "\t" 
                          << atom.element.name() << "\t" 
                          << atom.pos.x << "\t" 
                          << atom.pos.y << "\t" << atom.pos.z 
                          << std::endl;
            

            atom_index++;
            }
        residue_number++;
        }
    }
}