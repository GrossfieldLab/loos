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


#include <amber.hpp>
#include <exceptions.hpp>
#include <stdlib.h>
#include <ctype.h>

#include <iomanip>

namespace loos {


  // Parse simple Fortran format specifications, extracted from a %FORMAT tag...
  // Takes a string of expected format types (characters) that the extracted format
  // is compared against.  For example, to parse floats, expected could be any of "FEG".
  // If the format type does not match (for example, say "%FORMAT (20I8)" when expected
  // is "FEG"), then an error will be thrown.
  //
  // Format specs are converted to upper-case for validation

  Amber::FormatSpec Amber::parseFormat(const std::string& expected_types, const std::string& where) {

    reader.getNext();
    char period;
    FormatSpec fmt;

    // Verify line has a %FORMAT tag
    if (reader.line().compare(0, 7, "%FORMAT") != 0)
      throw(FileReadErrorWithLine(reader.name(), "Expected format for " + where, reader.lineNumber()));

    // Extract format spec between parens...
    boost::char_separator<char> sep("()");
    std::string input_line(reader.line());     // Ubuntu 16.04 apparently requires that the string in tokenizer
    tokenizer tokens(input_line, sep);         // be non-const
    tokenizer::iterator toks = tokens.begin();

    ++toks;
    std::istringstream iss(*toks);

    // try nXw.d first
    if (! (iss >> fmt.repeat >> fmt.type >> fmt.width >> period >> fmt.precision) ) {
      iss.clear();
      iss.seekg(0);
      // Now try nXw
      if (! (iss >> fmt.repeat >> fmt.type >> fmt.width) ) {
        iss.clear();
        iss.seekg(0);
        // Now try Xw
        if (! (iss >> fmt.type >> fmt.width) ) {
          iss.clear();
          iss.seekg(0);
          // And finally just try X
          if (! (iss >> fmt.type) )
            throw(FileReadErrorWithLine(reader.name(), "Cannot parse format for " + where, reader.lineNumber()));
        }
      }
    }

    // Compare against expected types
    std::string expected_types_UC = boost::to_upper_copy(expected_types);

    if (expected_types_UC.find_first_of(toupper(fmt.type)) == std::string::npos)
      throw(FileReadErrorWithLine(reader.name(), "Invalid format type for " + where, reader.lineNumber()));

    return(fmt);

  }


  void Amber::parseCharges() {
    FormatSpec fmt = parseFormat("EFG", "charges");

    std::vector<double> charges = readBlock<double>(fmt.width);
    if (charges.size() != atoms.size())
      throw(FileReadErrorWithLine(reader.name(), "Error parsing charges from amber file", reader.lineNumber()));

    for (uint i=0; i<charges.size(); ++i)
      atoms[i]->charge(charges[i]);
  }


  void Amber::parseMasses()  {
    FormatSpec fmt = parseFormat("EFG", "masses");

    std::vector<double> masses = readBlock<double>(fmt.width);
    if (masses.size() != atoms.size())
      throw(FileReadErrorWithLine(reader.name(), "Error parsing masses from amber file", reader.lineNumber()));

    for (uint i=0; i<masses.size(); ++i)
      atoms[i]->mass(masses[i]);

  }



  void Amber::parseResidueLabels() {
    FormatSpec fmt = parseFormat("a", "residue labels");

    std::vector<std::string> labels = readBlock<std::string>(fmt.width);
    if (labels.size() != nres)
      throw(FileReadErrorWithLine(reader.name(), "Error parsing residue labels from amber file", reader.lineNumber()));

    residue_labels = labels;
  }


  void Amber::parseResiduePointers() {
    FormatSpec fmt = parseFormat("I", "residue pointers");

    std::vector<uint> pointers = readBlock<uint>(fmt.width);
    if (pointers.size() != nres)
      throw(FileReadErrorWithLine(reader.name(), "Error parsing residue pointers from amber file", reader.lineNumber()));
    residue_pointers = pointers;
  }


  void Amber::assignResidues(void) {
    if (!(residue_pointers.size() == nres && residue_labels.size() == nres))
      throw(std::runtime_error("Unable to assign residues."));

    int curresid = 0;
    uint i = 0;
    std::string curresname;

    for (i=0; i<nres-1; i++) {
      ++curresid;
      curresname = residue_labels[i];
      for (uint j = residue_pointers[i]; j<residue_pointers[i+1]; j++) {
        atoms[j-1]->resid(curresid);
        atoms[j-1]->resname(curresname);
      }
    }

    // Fix the end...
    ++curresid;
    curresname = residue_labels[i];
    for (uint j = residue_pointers[i]-1; j<natoms; j++) {
      atoms[j]->resid(curresid);
      atoms[j]->resname(curresname);
    }
  }




  void Amber::parseBonds(const uint n) {
    FormatSpec fmt = parseFormat("I", "bonds");


    std::vector<int> bond_list = readBlock<int>(fmt.width);

    if (bond_list.size() != 3*n)
      throw(FileReadErrorWithLine(reader.name(), "Error parsing bonds in amber file", reader.lineNumber()));

    for (uint i=0; i<bond_list.size(); i += 3) {
      if (bond_list[i] == bond_list[i+1])
        continue;

      pAtom aatom = atoms[bond_list[i]/3];
      pAtom batom = atoms[bond_list[i+1]/3];

      // Amber bond lists are not symmetric, so make sure we add both pairs...
      if (!(aatom->isBoundTo(batom)))
        aatom->addBond(batom);
      if (!(batom->isBoundTo(aatom)))
        batom->addBond(aatom);

    }

  }


  void Amber::parsePointers() {
    FormatSpec fmt = parseFormat("I", "pointers");

    std::vector<uint> pointers = readBlock<uint>(fmt.width);

    // Now build up the atomic-group...
    if (atoms.size() != 0)
      throw(std::logic_error("Internal error: trying to read in an amber parmtop into a non-empty group!"));


    natoms = pointers[0];
    nbonh = pointers[2];
    mbona = pointers[3];
    nres = pointers[11];

    for (uint i=0; i<natoms; i++) {
      pAtom pa(new Atom);
      pa->id(i+1);
      pa->index(i);
      atoms.push_back(pa);
    }

  }


  // Simply slurp up the title (for now)
  void Amber::parseTitle() {

    FormatSpec fmt = parseFormat("a", "title");
    std::vector<std::string> titles = readBlock<std::string>(fmt.width);
    for (std::vector<std::string>::const_iterator i = titles.begin(); i != titles.end(); ++i)
      _title += *i;
  }


  void Amber::parseAtomNames() {
    FormatSpec fmt = parseFormat("a", "atom names");

    std::vector<std::string> names = readBlock<std::string>(fmt.width);
    if (names.size() != natoms)
      throw(FileReadErrorWithLine(reader.name(), "Error parsing atom names", reader.lineNumber()));
    for (uint i=0; i<names.size(); ++i)
      atoms[i]->name(names[i]);
  }

  void Amber::parseBoxDimensions() {
    FormatSpec fmt = parseFormat("E", "box dimensions");

    std::vector<double> dimensions = readBlock<double>(fmt.width);
    const double epsilon = 1e-8;

    double angle = dimensions[0];
    double diff = fabs(angle - 90.0);
    if (diff > epsilon) {
      std::cerr << "prmtop specifies non-rectangular box: ignoring"
                << std::endl;
      return;
    }

    GCoord box = GCoord(dimensions[1], dimensions[2], dimensions[3]);
    periodicBox(box);
  }

  void Amber::parseAmoebaRegularBondNumList() {
    FormatSpec fmt = parseFormat("I", "amoeba_regular_num_bond_list");
    reader.getNext();
    std::istringstream iss(reader.line());

    if (! (iss >> std::setw(fmt.width) >> _amoeba_regular_bond_num_list))
      throw(FileReadErrorWithLine(reader.name(), "Error parsing amoeba_regular_bond_num_list", reader.lineNumber()));
  }



  void Amber::parseAmoebaRegularBondList(const uint n) {
    FormatSpec fmt = parseFormat("I", "amoeba_regular_bond_list");


    std::vector<int> bond_list = readBlock<int>(fmt.width);

    if (bond_list.size() != 3*n)
      throw(FileReadErrorWithLine(reader.name(), "Error parsing amoeba bonds in amber file", reader.lineNumber()));

    for (uint i=0; i<bond_list.size(); i += 3) {
      if (bond_list[i] == bond_list[i+1])
        continue;

      // Amoeba bond indices appear not to be /3 as regular amber bonds are...
      // Are we sure???
      pAtom aatom = atoms[bond_list[i]-1];
      pAtom batom = atoms[bond_list[i+1]-1];

      // Amber bond lists are not symmetric, so make sure we add both pairs...
      if (!(aatom->isBoundTo(batom)))
        aatom->addBond(batom);
      if (!(batom->isBoundTo(aatom)))
        batom->addBond(aatom);

    }

  }



  void Amber::read(std::istream& ifs) {
    reader.stream(ifs);

    while (reader.getNext()) {

      boost::char_separator<char> sep(" \t");
      std::string input_line(reader.line());      // See parseFormat() for explanation...required by Ubuntu 16
      tokenizer tokens(input_line, sep);
      tokenizer::iterator toks = tokens.begin();

      if (*toks != "%FLAG")
        continue;

      ++toks;
      if (*toks == "TITLE")
        parseTitle();
      else if (*toks == "POINTERS")
        parsePointers();
      else if (*toks == "ATOM_NAME")
        parseAtomNames();
      else if (*toks == "CHARGE")
        parseCharges();
      else if (*toks == "MASS")
        parseMasses();
      else if (*toks == "RESIDUE_LABEL")
        parseResidueLabels();
      else if (*toks == "RESIDUE_POINTER")
        parseResiduePointers();
      else if (*toks == "BONDS_INC_HYDROGEN")
        parseBonds(nbonh);
      else if (*toks == "BONDS_WITHOUT_HYDROGEN")
        parseBonds(mbona);
      else if (*toks == "BOX_DIMENSIONS")
        parseBoxDimensions();
      else if (*toks == "AMOEBA_REGULAR_BOND_NUM_LIST")
        parseAmoebaRegularBondNumList();
      else if (*toks == "AMOEBA_REGULAR_BOND_LIST")
        parseAmoebaRegularBondList(_amoeba_regular_bond_num_list);

    }

    assignResidues();
    deduceAtomicNumberFromMass();
    setGroupConnectivity();
  }


}
