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

#include <mmcif.hpp>

namespace loos {

void MMCIF::read(std::istream& is) {
    OpenBabel::OBConversion obconversion(&is);
    obconversion.SetInFormat("cif");
    OpenBabel::OBMol mol;
    double *coords = new double[3];

    if (!obconversion.Read(&mol)) {
        throw(FileReadError(_fname), std::string("Error reading mmcif file"));
    }

    uint natoms = mol.NumAtoms();

    // TODO: Should probably loop over residues, then atoms, so we don't
    //       have to create a new residue object for each atom
    OpenBabel::OBAtom *a;
    uint index = 0;
    FOR_ATOMS_OF_MOL(a, mol) {
        coords = a->GetCoordinate();
        pAtom pa(new Atom);

        pa->index(index);
        pa->id(a->GetId());
        pa->coords(GCoord(coords[0], coords[1], coords[2]));

        pa->mass(a->GetAtomicMass());
        pa->atomic_number(a->GetAtomicNum());

        OpenBabel::OBResidue *residue = a->GetResidue();
        pa->name(residue->GetAtomID( &(*a) ));
        pa->resname(residue->GetName());
        pa->resid(residue->GetNum());
        pa->chainId(std::string(1, residue->GetChain()));

        // TODO: get bond info, if it exists
        // probably could do it more efficiently by separately looping over bonds
        /*
        for (auto iter = a->BeginBonds();
                  iter != a->EndBonds();
                  iter = a->NextBond()) {

                  }
        */

        append(pa);

        index++;
    }

    // TODO: this causes a segfault, but it looks like a leak to me
    //delete [] coords;
}



}
