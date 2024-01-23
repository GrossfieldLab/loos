/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2023, Tod D. Romo, Alan Grossfield
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

#include <mmcif.hpp>

namespace loos {

  MMCIF* MMCIF::clone(void) const {
    return(new MMCIF(*this));
  }

  void MMCIF::read(const std::string& filename) {
        auto structure = gemmi::read_structure_file(filename, gemmi::CoorFormat::Mmcif);
        auto unit_cell = structure.cell;
        auto box = loos::GCoord(unit_cell.a, unit_cell.b, unit_cell.c);

        periodicBox(box);
        
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
                    // TODO: charge is a char*, looks like it's usually "?" in actual
                    //       mmCIF files. Perhaps a try/catch block to convert to float?

                    pa->atomic_number(atom.element.atomic_number());
                    append(pa); 

                    atom_index++;
                }
            }
        }
    
    // assign the masses
    uint n_assigned = deduceMassFromAtomicNumber();
    // TODO: need to add bonds
    }

}


