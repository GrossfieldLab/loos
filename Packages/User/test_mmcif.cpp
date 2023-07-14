
#include <loos.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>    // cif::Document -> Structure
#include <gemmi/gz.hpp> 

namespace cif = gemmi::cif;


int main(int argc, char *argv[]) {
    std::string filename = std::string(argv[1]);
    
    loos::AtomicGroup ag;
    
    auto structure = gemmi::read_structure_file(filename, gemmi::CoorFormat::Mmcif);
    auto unit_cell = structure.cell;
    auto box = loos::GCoord(unit_cell.a, unit_cell.b, unit_cell.c);
    ag.periodicBox(box);

    // TODO: hard-wired to read first model, but there should probably be a way to read others
    auto model = structure.first_model();
    int atom_index = 0;
    int residue_number = 1;
    for (auto chain:model.chains) {
        std::string chain_name = chain.name;
        for (auto residue:chain.residues) {
            std::string residue_name = residue.name;
            auto label_seq = residue.label_seq;
            std::string res_entity_id = residue.entity_id;
            std::cout << "Residue: "
                      << residue_number << "\t"
                      << label_seq.str() << "\t"
                      << res_entity_id << "\t" 
                      << residue.name << "\t"
                      << residue.seqid.num.value << "\t"
                      << std::endl;
            for (auto atom:residue.atoms) {
                loos::pAtom pa(new loos::Atom);
                pa->index(atom.serial);
                pa->id(atom_index);
                pa->name(atom.name);
                pa->PDBelement(atom.element.name());
                pa->coords().x(atom.pos.x);
                pa->coords().y(atom.pos.y);
                pa->coords().z(atom.pos.z);
                pa->resid(residue.seqid.num.value);
                pa->chainId(chain_name);
                // might as well fill in segid, even though it's not official
                pa->segid(chain_name);
                pa->resname(residue.name);
                std::cout << atom.name << "\t" 
                          << atom.element.name() << "\t" 
                          << atom.pos.x << "\t" 
                          << atom.pos.y << "\t" 
                          << atom.pos.z << std::endl;
                // TODO: charge is a char*, looks like it's usually "?" in actual
                //       mmCIF files. Perhaps a try/catch block to convert to float?

                // TODO: since I've got the element, in principle I can look up the 
                //       mass and atomic number, but doing so will require some changes to
                //       AtomicNumberDeducer
                ag.append(pa); 

            atom_index++;
            }
        residue_number++;
        }
    }
    std::cout << ag << std::endl;
}