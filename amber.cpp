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

  void Amber::getNextLine(std::istream& is) {

    if (_unget) {
      _unget = false;
      return;
    }

    while (! getline(is, _current_line).eof() ) {
      ++_lineno;

      if (_current_line.compare(0, 8, "%COMMENT") != 0)
        break;

    }
  }


  Amber::FormatSpec Amber::parseFormat(std::istream& is, const std::string& expected_types, const std::string& where) {
    getNextLine(is);
    char period;
    FormatSpec fmt;

    // Verify line has a %FORMAT tag
    if (_current_line.compare(0, 7, "%FORMAT") != 0)
      throw(FileParseError("Expected format for " + where, _lineno));

    // Extract format spec between parens...
    boost::char_separator<char> sep("()");
    tokenizer tokens(_current_line, sep);
    tokenizer::iterator toks = tokens.begin();

    ++toks;
    std::istringstream iss(*toks);

    // Attempt to parse...
    if (! (iss >> fmt.repeat >> fmt.type >> fmt.width >> period >> fmt.precision) ) {
      iss.clear();
      iss.seekg(0);
      if (! (iss >> fmt.repeat >> fmt.type >> fmt.width) ) {
        iss.clear();
        iss.seekg(0);
        if (! (iss >> fmt.type >> fmt.width) ) {
          iss.clear();
          iss.seekg(0);
          if (! (iss >> fmt.type) )
            throw(FileParseError("Cannot parse format for " + where, _lineno));
        }
      }
    }
    
    if (expected_types.find_first_of(fmt.type) == std::string::npos)
      throw(FileParseError("Invalid format type for " + where, _lineno));

    return(fmt);

  }


  void Amber::parseCharges(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "EFG", "charges");

    std::vector<double> charges = readBlock<double>(is, fmt.width);
    if (charges.size() != atoms.size())
      throw(FileParseError("Error parsing charges from amber file", _lineno));
    
    for (uint i=0; i<charges.size(); ++i)
      atoms[i]->charge(charges[i]);
  }


  void Amber::parseMasses(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "EFG", "masses");

    std::vector<double> masses = readBlock<double>(is, fmt.width);
    if (masses.size() != atoms.size())
      throw(FileParseError("Error parsing masses from amber file", _lineno));
    
    for (uint i=0; i<masses.size(); ++i)
      atoms[i]->mass(masses[i]);

  }



  void Amber::parseResidueLabels(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "a", "residue labels");

    residue_labels = readBlock<std::string>(is, fmt.width);
    if (residue_labels.size() != nres)
      throw(FileParseError("Error parsing residue labels from amber file", _lineno));

  }


  void Amber::parseResiduePointers(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "I", "residue pointers");

    residue_pointers = readBlock<uint>(is, fmt.width);
    if (residue_pointers.size() != nres)
      throw(FileParseError("Error parsing residue pointers from amber file", _lineno));
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




  void Amber::parseBonds(std::istream& is, const uint n) {
    FormatSpec fmt = parseFormat(is, "I", "bonds");

  
    std::vector<int> bond_list = readBlock<int>(is, fmt.width);

    if (bond_list.size() != 3*n)
      throw(FileParseError("Error parsing bonds in amber file", _lineno));
    
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


  void Amber::parsePointers(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "I", "pointers");

    std::vector<uint> pointers = readBlock<uint>(is, fmt.width);

    natoms = pointers[0];
    nbonh = pointers[2];
    mbona = pointers[3];
    nres = pointers[11];

    // Now build up the atomic-group...
    if (atoms.size() != 0)
      throw(std::logic_error("Internal error: trying to read in an amber parmtop into a non-empty group!"));

    for (uint i=0; i<natoms; i++) {
      pAtom pa(new Atom);
      pa->id(i+1);
      atoms.push_back(pa);
    }

  }


  // Simply slurp up the title (for now)
  void Amber::parseTitle(std::istream& is) {

    FormatSpec fmt = parseFormat(is, "a", "title");
    std::vector<std::string> titles = readBlock<std::string>(is, fmt.width);
    for (std::vector<std::string>::const_iterator i = titles.begin(); i != titles.end(); ++i)
      _title += *i;
  }


  void Amber::parseAtomNames(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "a", "atom names");
    
    std::vector<std::string> names = readBlock<std::string>(is, fmt.width);
    if (names.size() != natoms)
      throw(FileParseError("Error parsing atom names", _lineno));
    for (uint i=0; i<names.size(); ++i)
      atoms[i]->name(names[i]);
  }

  void Amber::parseAmoebaRegularBondNumList(std::istream& is) {
    FormatSpec fmt = parseFormat(is, "I", "amoeba_regular_num_bond_list");
    getNextLine(is);
    std::istringstream iss(_current_line);

    if (! (iss >> std::setw(fmt.width) >> _amoeba_regular_bond_num_list))
      throw(FileParseError("Error parsing amoeba_regular_bond_num_list", _lineno));
  }

  void Amber::parseAmoebaRegularBondList(std::istream& is, const uint n) {
    FormatSpec fmt = parseFormat(is, "I", "amoeba_regular_bond_list");

  
    std::vector<int> bond_list = readBlock<int>(is, fmt.width);

    if (bond_list.size() != 3*n)
      throw(FileParseError("Error parsing amoeba bonds in amber file", _lineno));
    
    for (uint i=0; i<bond_list.size(); i += 3) {
      if (bond_list[i] == bond_list[i+1])
        continue;

      pAtom aatom = atoms[bond_list[i]-1];
      pAtom batom = atoms[bond_list[i+1]-1];
    
      // Amber bond lists are not symmetric, so make sure we add both pairs...
      if (!(aatom->isBoundTo(batom)))
        aatom->addBond(batom);
      if (!(batom->isBoundTo(aatom)))
        batom->addBond(aatom);
      
    }
    
  }



  void Amber::read(std::istream& is) {

    while (true) {
      getNextLine(is);
      if (is.eof())
        break;
      if (is.fail())
        throw(FileParseError("Error reading amber file ", _lineno));

      boost::char_separator<char> sep(" \t");
      tokenizer tokens(_current_line, sep);
      tokenizer::iterator toks = tokens.begin();
      
      if (*toks != "%FLAG")
        continue;

      ++toks;
      if (*toks == "TITLE")
        parseTitle(is);
      else if (*toks == "POINTERS")
        parsePointers(is);
      else if (*toks == "ATOM_NAME")
        parseAtomNames(is);
      else if (*toks == "CHARGE")
        parseCharges(is);
      else if (*toks == "MASS")
        parseMasses(is);
      else if (*toks == "RESIDUE_LABEL")
        parseResidueLabels(is);
      else if (*toks == "RESIDUE_POINTER")
        parseResiduePointers(is);
      else if (*toks == "BONDS_INC_HYDROGEN")
        parseBonds(is, nbonh);
      else if (*toks == "BONDS_WITHOUT_HYDROGEN")
        parseBonds(is, mbona);
      else if (*toks == "AMOEBA_REGULAR_BOND_NUM_LIST")
        parseAmoebaRegularBondNumList(is);
      else if (*toks == "AMOEBA_REGULAR_BOND_LIST")
        parseAmoebaRegularBondList(is, _amoeba_regular_bond_num_list);
      
    }

    assignResidues();
    deduceAtomicNumberFromMass();
    setGroupConnectivity();
  }


}

