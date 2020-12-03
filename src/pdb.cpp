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


#include <pdb.hpp>
#include <utils.hpp>
#include <Fmt.hpp>

#include <iomanip>
#include <boost/unordered_set.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>


namespace loos {

  // Assume we're only going to find spaces in a PDB file...
  bool PDB::emptyString(const std::string& s) {
    std::string::const_iterator i;

    for (i = s.begin(); i != s.end(); ++i)
      if (*i != ' ')
        return(false);

    return(true);
  }


  // Special handling for REMARKs to ignore the line code, if
  // present...

  void PDB::parseRemark(const std::string& s) {
    std::string t;

    if (s[6] == ' ' && isdigit(s[7]))
      t = s.substr(11, 58);
    else
      t = s.substr(7, 62);

    _remarks.add(t);
  }


  // Parse an ATOM or HETATM record...
  // Note: ParseErrors can come from parseStringAs

  void PDB::parseAtomRecord(const std::string& s) {
    greal r;
    gint i;
    std::string t;
    GCoord c;
    pAtom pa(new Atom);

    pa->index(_max_index++);

    t = parseStringAs<std::string>(s, 0, 6);
    pa->recordName(t);

    i = parseStringAsHybrid36(s, 6, 5);
    pa->id(i);

    t = parseStringAs<std::string>(s, 12, 4);
    pa->name(t);

    t = parseStringAs<std::string>(s, 16, 1);
    pa->altLoc(t);

    t = parseStringAs<std::string>(s, 17, 4);
    pa->resname(t);

    t = parseStringAs<std::string>(s, 21, 1);
    pa->chainId(t);

    i = parseStringAsHybrid36(s, 22, 4);
    pa->resid(i);

    t = parseStringAs<std::string>(s, 26, 1);

    // Special handling of resid field since it may be frame-shifted by
    // 1 col in some cases...
    if (strictness_policy) {
      char c = t[0];

      if (c != ' ' && !isalpha(c))
        throw(ParseError("Non-alpha character in iCode column of PDB"));
    } else {
      char c = t[0];

      // Assume that if we see this variant, then we're not using hybrid-36
      if (c != ' ' && isdigit(c)) {
        i = parseStringAs<int>(s, 22, 5);
        pa->resid(i);
        t = " ";
      }

    }
    pa->iCode(t);

    c[0] = parseStringAs<float>(s, 30, 8);
    c[1] = parseStringAs<float>(s, 38, 8);
    c[2] = parseStringAs<float>(s, 46, 8);
    pa->coords(c);

    if (s.size() > 54) {
      r = parseStringAs<float>(s, 54, 6);
      pa->occupancy(r);

      if (s.size() > 60) {
	r = parseStringAs<float>(s, 60, 6);
	pa->bfactor(r);

	if (s.size() > 72) {
	  t = parseStringAs<std::string>(s, 72, 4);
	  pa->segid(t);

	  if (s.size() > 76) {
	    t = parseStringAs<std::string>(s, 76, 2);
	    pa->PDBelement(t);

	    // Charge is not currently handled...
	    // t = parseStringAs<std::string>(s, 78, 2);
	  }
	} else { // segid
	  _missing_segid = true;
	}
      } else { // b-factor
	_missing_b = _missing_segid = true;
      }
    } else { // occupancies
      _missing_q = _missing_b = _missing_segid = true;
    }
    append(pa);

    // Record which pAtom belongs to this atomid.
    // NOTE: duplicate atomid's are NOT checked for
    _atomid_to_patom[pa->id()] = pa;
  }



  // Convert an Atom to a string with a PDB format...

  std::string PDB::atomAsString(const pAtom p) const {
    std::ostringstream s;

    // Float formatter for coords
    Fmt crdfmt(3);
    crdfmt.width(8);
    crdfmt.right();
    crdfmt.trailingZeros(true);
    crdfmt.fixed();

    // Float formatter for B's and Q's...
    Fmt bqfmt(2);
    bqfmt.width(6);
    bqfmt.right();
    bqfmt.trailingZeros(true);
    bqfmt.fixed();

    // We don't worry about strings exceeding field-widths (yet),
    // but do check for numeric overflows...
    s << std::setw(6) << std::left << p->recordName();
    s << hybrid36AsString(p->id(), 5) << " ";
    s << std::setw(4) << std::left << p->name();

    s << std::setw(1) << p->altLoc();
    s << std::setw(4) << std::left << p->resname();

    s << std::setw(1) << std::right << p->chainId();
    s << hybrid36AsString(p->resid(), 4);
    s << std::setw(2) << p->iCode();
    s << "  ";

    s << crdfmt(p->coords().x());
    s << crdfmt(p->coords().y());
    s << crdfmt(p->coords().z());
    s << bqfmt(p->occupancy());
    s << bqfmt(p->bfactor());
    s << "      ";
    s << std::setw(4) << std::left << p->segid();
    s << std::setw(2) << std::right << p->PDBelement();
    if (_show_charge)
      s << std::setw(2) << p->charge();
    else
      s << "  ";

    return(s.str());
  }


  // Private function to search the map of atomid's -> pAtoms
  // Throws an error if the atom is not found
  pAtom PDB::findAtom(const int id) {
    std::map<int, pAtom>::iterator i = _atomid_to_patom.find(id);
    if (i == _atomid_to_patom.end()) {
      std::ostringstream oss;
      oss << "Cannot find atom corresponding to atomid " << id << " for making a bond.";
      throw(LOOSError(oss.str()));
    }
    return(i->second);
  }




  // If the PDB already had symmetrical CONECT records, then it will
  // have duplicates thanks to how parseConectRecord works.  We sort
  // and filter out repeated bonds to fix this.
  void PDB::uniqueBonds() {

    for (iterator atom = begin(); atom != end(); ++atom) {
      std::vector<int> bonds = (*atom)->getBonds();

      if (!bonds.empty()) {
        std::sort(bonds.begin(), bonds.end());
        std::vector<int> unique_bonds;

        unique_bonds.push_back(bonds[0]);
        for (std::vector<int>::const_iterator i = bonds.begin()+1; i != bonds.end(); ++i)
          if (*i != *(i-1))
            unique_bonds.push_back(*i);

        (*atom)->setBonds(unique_bonds);
      }
    }
  }


  // Parse CONECT records, updating the referenced atoms...
  // Couple of issues:
  //
  //    Will accept up to 8 bound atoms and considers them all equal,
  // but the PDB standard says some are h-bonded and some are
  // salt-bridged...
  //
  //    No check is made for overflow of fields...
  //
  //    Adding a bond requires finding the bound atoms, which can
  // be expensive, so we build-up a map of atomid to pAtoms as the
  // ATOM records are being parsed.  This is then searched to find
  // the pAtom that matches the CONECT atomid (rather than using
  // findById()).

  void PDB::parseConectRecord(const std::string& s) {
    int bound_id = parseStringAsHybrid36(s, 6, 5);


    // Rely on findAtom to throw if bound_id not found...
    pAtom bound = findAtom(bound_id);

    // This currently includes fields designated as H-bond indices...
    // Should we do this? or separate them out?  Hmmm...
    for (int i=0; i<8; ++i) {
      int j = i * 5 + 11;
      if (static_cast<uint>(j) >= s.length())
        break;
      std::string t = s.substr(j, 5);
      if (emptyString(t))
        break;
      int id = parseStringAsHybrid36(t);

      // findAtom will throw if id is not found
      pAtom boundee = findAtom(id);
      bound->addBond(boundee);
      boundee->addBond(bound);
    }
  }


  void PDB::parseCryst1Record(const std::string& s) {
    greal r;
    gint i;
    std::string t;

    UnitCell newcell;

    r = parseStringAs<float>(s, 6, 9);
    newcell.a(r);
    r = parseStringAs<float>(s, 15, 9);
    newcell.b(r);
    r = parseStringAs<float>(s, 24, 9);
    newcell.c(r);

    r = parseStringAs<float>(s, 33, 7);
    newcell.alpha(r);
    r = parseStringAs<float>(s, 40, 7);
    newcell.beta(r);
    r = parseStringAs<float>(s, 47, 7);
    newcell.gamma(r);

    // Special handling in case of mangled CRYST1 record...
    if (s.length() < 66) {
      t = s.substr(55);

      newcell.spaceGroup(t);
      newcell.z(-1);   // ???
    } else {
      t = parseStringAs<std::string>(s, 55, 11);
      newcell.spaceGroup(t);
      i = parseStringAs<int>(s, 66, 4);
      newcell.z(i);
    }

    cell = newcell;
    _has_cryst = true;
  }


  //! Top level parser...
  //! Reads a PDB from an input stream
  /*
   * Will transform any caught exceptions into a FileReadError
   */
  void PDB::read(std::istream& is) {
    std::string input;
    bool has_cryst = false;
    bool has_bonds = false;
    boost::unordered_set<std::string> seen;

    while (getline(is, input)) {
      try {
	if (input.substr(0, 4) == "ATOM" || input.substr(0,6) == "HETATM")
	  parseAtomRecord(input);
	else if (input.substr(0, 6) == "REMARK")
	  parseRemark(input);
	else if (input.substr(0,6) == "CONECT") {
	  has_bonds = true;
	  parseConectRecord(input);
	} else if (input.substr(0, 6) == "CRYST1") {
	  parseCryst1Record(input);
	  has_cryst = true;
	} else if (input.substr(0,3) == "TER")
	  ;
	else if (input.substr(0,3) == "END")
	  break;
	else {
	  int space = input.find_first_of(' ');
	  std::string record = input.substr(0, space);
	  if (seen.find(record) == seen.end()) {
	    std::cerr << "Warning - unknown PDB record '" << record << "'" << std::endl;
	    seen.insert(record);
	  }
	}
      }
      catch(LOOSError& e) {
	throw(FileReadError(_fname, e.what()));
      }
      catch(...) {
	throw(FileReadError(_fname, "Unknown exception"));
      }
    }
    if (isMissingFields())
      std::cerr << "Warning- PDB is missing fields.  Default values will be used.\n";

    // Clean-up temporary storage...
    _atomid_to_patom.clear();

    // Do some post-extraction...
    if (loos::remarksHasBox(_remarks)) {
      GCoord c;
      try {
	c = loos::boxFromRemarks(_remarks);
      }
      catch(ParseError& e) {
	throw(FileReadError(_fname, e.what()));
      }
      periodicBox(c);
    } else if (has_cryst) {
      GCoord c(cell.a(), cell.b(), cell.c());
      periodicBox(c);
    }

    // Force atom id's to be monotonic if there was an overflow event...
    if (atoms.size() >= 100000)
      renumber();

    // Set bonds state...
    if (has_bonds) {
      setGroupConnectivity();
      uniqueBonds();
    }

  }


  std::ostream& FormattedUnitCell(std::ostream& os, const UnitCell& u) {
    os << "CRYST1";
    Fmt dists(3);
    dists.width(9).right().trailingZeros(true).fixed();
    Fmt angles(2);
    angles.width(7).right().trailingZeros(true).fixed();

    os << dists(u.a()) << dists(u.b()) << dists(u.c());
    os << angles(u.alpha()) << angles(u.beta()) << angles(u.gamma());
    os << " " << std::setw(11) << std::left << u.spaceGroup() << std::setw(4) << u.z();

    return(os);
  }

  std::ostream& XTALLine(std::ostream& os, const GCoord& box) {
    os << "REMARK  XTAL "
       << box.x() << " "
       << box.y() << " "
       << box.z();
    return(os);
  }


  std::ostream& FormatConectRecords(std::ostream& os, const PDB& p) {
    AtomicGroup::iterator ci;

    // Copy and sort
    PDB sorted = p;
    sorted.sort();

    for (ci = sorted.begin(); ci != sorted.end(); ++ci) {
      if ((*ci)->checkProperty(Atom::bondsbit)) {
        int id = (*ci)->id();

        std::vector<int> bonds = (*ci)->getBonds();
        if (bonds.empty())
          continue;


        // Filter bonds...  Any bond to an atomid smaller than the current atom
        // will have already been written, so omit.
        std::sort(bonds.begin(), bonds.end());
        std::vector<int> filtered_bonds;
        for (std::vector<int>::const_iterator i = bonds.begin(); i != bonds.end(); ++i)
          if (*i >= id)
            filtered_bonds.push_back(*i);

        if (filtered_bonds.empty())
          continue;

        os << "CONECT" << hybrid36AsString(id, 5);
        int i = 0;

        std::vector<int>::const_iterator cj;
        for (cj = filtered_bonds.begin(); cj != filtered_bonds.end(); ++cj) {
          if (++i > 4) {
            i = 1;
            os << "\nCONECT" << hybrid36AsString(id, 5);
          }
          int bound_id = *cj;
          pAtom pa = sorted.findById(bound_id);
          if (pa != 0)
            os << hybrid36AsString(bound_id, 5);
        }
        os << std::endl;
      }
    }

    return(os);
  }


  //! Output the group as a PDB...
  /**
   * There are some formatting changes that occur when the group has a
   * large number of atoms or resids.  The most significant is when
   * you have 100,000 or more, in which case you lose the altloc and
   * chainid fields on output.  However, the output PDB will load into
   * pymol...
   *
   */
  std::ostream& operator<<(std::ostream& os, const PDB& p) {
    AtomicGroup::const_iterator i;
    //std::vector<AtomicGroup>::const_iterator m;

    os << p._remarks;
    if (p.isPeriodic())
      XTALLine(os, p.periodicBox()) << std::endl;
    if (p._has_cryst)
      FormattedUnitCell(os, p.cell) << std::endl;

    if (!p.hasBonds()) {
      for (i = p.atoms.begin(); i != p.atoms.end(); ++i) {
        os << p.atomAsString(*i) << std::endl;
      }
    }
    else if (p._auto_ter){
      std::vector<AtomicGroup> molecules = p.splitByMolecule();
      for (auto m = molecules.begin(); m != molecules.end(); ++m) {
        PDB tmppdb = PDB::fromAtomicGroup(*m);
        for (i = tmppdb.atoms.begin(); i != tmppdb.atoms.end(); ++i) {
          os << p.atomAsString(*i) << std::endl;
        }
        os << "TER  " << std::endl;
      }
    }


    if (p.hasBonds()) {
      int maxid = 0;
      for (i = p.atoms.begin(); i != p.atoms.end(); ++i)
        if ((*i)->id() > maxid)
          maxid = (*i)->id();

      FormatConectRecords(os, p);
    }

    return(os);
  }

  PDB PDB::copy(void) const {
    AtomicGroup grp = this->AtomicGroup::copy();
    PDB p(grp);

    p._show_charge = _show_charge;
    p._auto_ter = _auto_ter;
    p._has_cryst = _has_cryst;
    p._remarks = _remarks;
    p.cell = cell;

    return(p);
  }

  PDB* PDB::clone(void) const {
    return(new PDB(*this));
  }

  PDB PDB::fromAtomicGroup(const AtomicGroup& g) {
    PDB p(g);

    if (p.isPeriodic())
      p.unitCell(UnitCell(p.periodicBox()));

    return(p);
  }


  bool PDB::showCharge(void) const { return(_show_charge); }
  void PDB::showCharge(bool b) { _show_charge = b; }
  bool PDB::strict(void) const { return(strictness_policy); }
  void PDB::strict(const bool b) { strictness_policy = b; }

  bool PDB::autoTerminate(void) const { return(_auto_ter); }

  void PDB::autoTerminate(bool b) { _auto_ter = b; }

  Remarks& PDB::remarks(void) { return(_remarks); }
  void PDB::remarks(const Remarks& r) { _remarks = r; }

  const UnitCell& PDB::unitCell(void) { return(cell); }
  void PDB::unitCell(const UnitCell& c) { _has_cryst = true; cell = c; }



}
