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

    // Turn off warnings from OpenBabel, just do errors
    OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obError);
    // Read the file
    if (!obconversion.Read(&mol)) {
        throw(FileReadError(_fname), std::string("Error reading mmcif file"));
    }

    // TODO: Should probably loop over residues, then atoms, so we don't
    //       have to create a new residue object for each atom
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

        append(pa);
        index++;

        // Record which pAtom belongs to this atomid.
        // NOTE: duplicate atomid's are NOT checked for
        _atomid_to_patom[pa->id()] = pa;
    }

    uint begin_bonded, end_bonded;
    bool has_bonds = false;
    FOR_BONDS_OF_MOL(b, mol) {
        has_bonds = true;
        begin_bonded = b->GetBeginAtomIdx();
        end_bonded = b->GetEndAtomIdx();

        // These will throw if the atoms aren't found
        pAtom first = findAtom(begin_bonded);
        pAtom second = findAtom(end_bonded);
        first->addBond(second);
        second->addBond(first);
    }
    if (has_bonds) {
        setGroupConnectivity();
    }

    // Don't need this info anymore
    _atomid_to_patom.clear();

    // Get unit cell information, if it's there
    if (mol.HasData("UnitCell")) {
        OpenBabel::OBUnitCell *unitcell = static_cast<OpenBabel::OBUnitCell *>
                                          (mol.GetData(OpenBabel::OBGenericDataType::UnitCell));
        cell.a(unitcell->GetA());
        cell.b(unitcell->GetB());
        cell.c(unitcell->GetC());
        cell.alpha(unitcell->GetAlpha());
        cell.beta(unitcell->GetBeta());
        cell.gamma(unitcell->GetGamma());
        cell.spaceGroup(unitcell->GetSpaceGroupName());
        _has_cryst = true;

        // Use the cell to set the periodic box
        // Note: this is wrong if it's not a simple orthonormal space group.
        //
        GCoord c(cell.a(), cell.b(), cell.c());
        periodicBox(c);
    }

    // TODO: this causes a segfault, but it looks like a leak to me
    //delete[] coords;
}


MMCIF MMCIF::copy(void) const {
  AtomicGroup grp = this->AtomicGroup::copy();
  MMCIF p(grp);

  p.cell = cell;

  return(p);
}

MMCIF* MMCIF::clone(void) const {
  return(new MMCIF(*this));
}


MMCIF MMCIF::fromAtomicGroup(const AtomicGroup& g) {
  MMCIF p(g);

  if (p.isPeriodic())
    p.unitCell(UnitCell(p.periodicBox()));

  return(p);
}

// Private function to search the map of atomid's -> pAtoms
// Throws an error if the atom is not found
pAtom MMCIF::findAtom(const uint id) {
  std::map<uint, pAtom>::iterator i = _atomid_to_patom.find(id);
  if (i == _atomid_to_patom.end()) {
    std::ostringstream oss;
    oss << "Cannot find atom corresponding to atomid " << id << " for making a bond.";
    throw(LOOSError(oss.str()));
  }
  return(i->second);
}

const UnitCell& MMCIF::unitCell(void) { return(cell); }
void MMCIF::unitCell(const UnitCell& c) { _has_cryst = true; cell = c; }

OpenBabel::OBMol MMCIF::toOpenBabel(void) const {
    /*
    Looking at https://github.com/openbabel/openbabel/blob/master/src/mol.cpp,
    line 1227 (the operator= method) to understand what it expects
    */

    OpenBabel::OBMol obmol;
    obmol.BeginModify();

    // loop over atoms
    for (auto a = begin(); a != end(); ++a) {
    //for (uint i=0; i<size(); ++i) {
        OpenBabel::OBAtom atom;

        atom.SetIdx((*a)->index());
        atom.SetId(static_cast<unsigned long>((*a)->id()));

        double x = (*a)->coords()[0];
        double y = (*a)->coords()[1];
        double z = (*a)->coords()[2];
        atom.SetVector(x, y, z);

        atom.SetType((*a)->name());
        atom.SetAtomicNum((*a)->atomic_number());

        //atom.SetChain((*a)->chainId().c_str());


        // set partial charges?
        atom.SetPartialCharge((*a)->charge());
        if (!obmol.AddAtom(atom)) {
            // Adding atom failed
            throw(LOOSError("Error adding atom to MMCIF"));
        }

    }

    // loop over residues

    // loop over bonds?

    // unit cell
    OpenBabel::OBUnitCell unitcell;
    unitcell.SetData(cell.a(), cell.b(), cell.c(),
                     cell.alpha(), cell.beta(), cell.gamma());
    obmol.SetData(&unitcell);


    // total Charge
    obmol.SetTotalCharge(this->totalCharge());

    return(obmol);
}

    //! Output the group as an MMCIF
    std::ostream& operator<<(std::ostream os, const MMCIF& m) {
        OpenBabel::OBMol om = m.toOpenBabel();
        OpenBabel::OBConversion conv(NULL, &os);
        conv.SetOutFormat("CIF");
        conv.Write(&om);

        return(os);
    }

}
