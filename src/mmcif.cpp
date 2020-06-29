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
    obconversion.SetInFormat("MMCIF");
    OpenBabel::OBMol mol;
    double *coords = new double[3];

    // Turn off warnings from OpenBabel, just do errors
    OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obError);
    // Read the file
    if (!obconversion.Read(&mol)) {
        throw(FileReadError(_fname), std::string("Error reading mmcif file"));
    }

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

#if 0
OpenBabel::OBMol * MMCIF::toOpenBabel(void) const {
    /*
    Looking at https://github.com/openbabel/openbabel/blob/master/src/mol.cpp,
    line 1227 (the operator= method) to understand what it expects
    */

    OpenBabel::OBMol *obmol = new OpenBabel::OBMol;
    obmol->BeginModify();

    bool has_charges = false;

    // loop over atoms
    for (auto a = begin(); a != end(); ++a) {
        //OpenBabel::OBAtom * atom = new OpenBabel::OBAtom;
        OpenBabel::OBAtom * atom = obmol->NewAtom();

        //atom->SetIdx((*a)->index());
        atom->SetId(static_cast<unsigned long>((*a)->id()));

        double x = (*a)->coords()[0];
        double y = (*a)->coords()[1];
        double z = (*a)->coords()[2];
        atom->SetVector(x, y, z);

        //atom.SetChain((*a)->chainId().c_str());


        // set partial charges?
        if ((*a)->checkProperty(Atom::chargebit)) {
            atom->SetPartialCharge((*a)->charge());
            has_charges = true;
        }

        /*
        if (!obmol->AddAtom(*atom)) {
            // Adding atom failed
            throw(LOOSError("Error adding atom to MMCIF"));
        }
        */

    }

    // loop over residues

    uint atom_index = 0;
    std::vector<AtomicGroup> residues = splitByResidue();
    for (uint i=0; i < residues.size(); ++i) {
        pAtom a = residues[i][0];
        OpenBabel::OBResidue *res = obmol->NewResidue();
        res->SetIdx(i);
        res->SetNum(a->resid());
        res->SetName(a->resname());

        for (uint j=0; j<residues[i].size(); ++j) {
            OpenBabel::OBAtom * babel_atom = obmol->GetAtom(atom_index);
            pAtom loos_atom = residues[i][j];
            res->AddAtom(babel_atom);
            //res->SetHetAtom(babel_atom, true);

            // weird, but I think I need to set the metadata here
            res->SetAtomID(babel_atom, std::string(loos_atom->name()));

            // Only set the atomic number if we know it
            if (loos_atom->atomic_number() > 0) {
                babel_atom->SetAtomicNum(loos_atom->atomic_number());
            }

            atom_index++;
        }
    }

    // loop over bonds?


    // unit cell
    OpenBabel::OBUnitCell *unitcell = new OpenBabel::OBUnitCell;
    if (_has_cryst) {
        unitcell->SetData(cell.a(), cell.b(), cell.c(),
                          cell.alpha(), cell.beta(), cell.gamma());
        unitcell->SetSpaceGroup(cell.spaceGroup());
        obmol->SetData(unitcell);
    } else if (isPeriodic()) {
        GCoord c = periodicBox();
        unitcell->SetData(c.x(), c.y(), c.z(), 90., 90., 90.);
        unitcell->SetSpaceGroup(std::string("P1"));
        obmol->SetData(unitcell);
    }
    else {
        delete(unitcell);
    }

    // total Charge
    if (has_charges) {
        obmol->SetTotalCharge(this->totalCharge());
    }

    // TODO:: state of debugging
    /*
        When I loop over residues, the metadata is there -- I get correct
        atom and residue names.
        However, in the immediately following loop, residue DOES NOT have
        names, ID, etc, and I have no friggin' idea why.
     */

    FOR_RESIDUES_OF_MOL(r, obmol) {
        std::cerr << r->GetName() << "\t"
                  << r->GetNum() << std::endl;
        OpenBabel::OBResidue *residue = &(*r);
        FOR_ATOMS_OF_RESIDUE(a2, residue) {
            std::cerr << r->GetAtomID(&(*a2)) << std::endl;
        }
    }

    FOR_ATOMS_OF_MOL(a, obmol) {

        std::cerr << a->GetId() << std::endl;
        std::cerr << a->GetAtomicMass() << std::endl;
        std::cerr << a->GetAtomicNum() << std::endl;

        OpenBabel::OBResidue *residue = a->GetResidue();
        std::cerr << residue->GetAtomID( &(*a) ) << std::endl;
        std::cerr << residue->GetName() << std::endl;
        std::cerr << residue->GetNum() << std::endl;

        std::cerr << std::endl;
    }


    obmol->EndModify();

    return(obmol);
#endif

#if 0
    //! Output the group as an MMCIF
    std::ostream& operator<<(std::ostream& os, const MMCIF& m) {
        OpenBabel::OBMol * om = m.toOpenBabel();
        std::cerr << "In operator << " << std::endl;
        FOR_ATOMS_OF_MOL(a, om) {

            std::cerr << a->GetId() << std::endl;
            std::cerr << a->GetAtomicMass() << std::endl;
            std::cerr << a->GetAtomicNum() << std::endl;

            OpenBabel::OBResidue *residue = a->GetResidue();
            std::cerr << residue->GetAtomID( &(*a) ) << std::endl;
            std::cerr << residue->GetName() << std::endl;
            std::cerr << residue->GetNum() << std::endl;

        }


        std::ifstream dummy("/dev/null");
        OpenBabel::OBConversion conv(&dummy, &os);

        if (!conv.SetOutFormat("MMCIF")) {
            std::cerr << "Alan is an idiot" << std::endl;
        }
        conv.Write(om, &os);

        return os;
    }
#endif

    //! Output the group as an MMCIF
    /**
        Code is derived from the OpenBabel MMCIF writing code.
    */
    std::ostream& operator<<(std::ostream& os, const MMCIF& m) {

        os << "# - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        os << "# " << std::endl;
        os << "# MMCIF file generated by LOOS " << std::endl;
        os << "# " << std::endl;
        os << "# - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        os << std::endl;

        // Identify of the system -- for now, fill with junk. All instances
        // where this is true marked with TODO: foobar
        std::string id = std::string("foobar");

        os << "data_" << id << std::endl; // TODO: foobar
        os << std::endl;

        // Entry ID: again, for now filled with junk
        os << "###########" << std::endl;
        os << "## ENTRY ##" << std::endl;
        os << "###########" << std::endl;
        os << std::endl;
        os << "_entry.id\t" << id << std::endl; // TODO: foobar
        os << std::endl;

        // I'm skipping the CHEMICAL and CHEMICAL FORMULA entries, since they
        // appear to be optional.

        os << "###############" << std::endl;
        os << "## ATOM_SITE ##" << std::endl;
        os << "###############" << std::endl;
        os << std::endl;
        os << "loop_" << std::endl;
        os << "_atom_site.id" << std::endl;
        os << "_atom_site.type_symbol" << std::endl;
        os << "_atom_site.label_atom_id" << std::endl;
        os << "_atom_site.label_comp_id" << std::endl;
        os << "_atom_site.label_entity_id" << std::endl;
        os << "_atom_site.label_seq_id" << std::endl;
        os << "_atom_site.Cartn_x" << std::endl;
        os << "_atom_site.Cartn_y" << std::endl;
        os << "_atom_site.Cartn_z" << std::endl;

        for (auto atom = m.begin(); atom != m.end(); ++atom) {
            pAtom a = (*atom);
            os << "\t" << a->id()
               << "\t" << a->name()[0]  // really, this should be the element name
               << "\t " << a->name()  // the space after the tab helps pymol
               << "\t" << a->resname()
               << "\t" << a->segid()  // really, this is the chain number
               << "\t" << a->resid()
               << "\t" << a->coords().x()
               << "\t" << a->coords().y()
               << "\t" << a->coords().z()
               << std::endl;
        }

    if (m.isPeriodic()) {
        os << std::endl;
        os << "##########" << std::endl;
        os << "## CELL ##" << std::endl;
        os << "##########" << std::endl;
        os << std::endl;
        os << "_cell.entry_id\t" << id << std::endl;  // TODO: foobar
        os << "_cell.length_a\t" << m.cell.a() << std::endl;
        os << "_cell.length_b\t" << m.cell.b() << std::endl;
        os << "_cell.length_c\t" << m.cell.c() << std::endl;
        os << "_cell.angle_alpha\t" << m.cell.alpha() << std::endl;
        os << "_cell.angle_beta\t"  << m.cell.beta() << std::endl;
        os << "_cell.angle_gamma\t" << m.cell.gamma() << std::endl;
        os << std::endl;
    }

    return(os);
    }

}
