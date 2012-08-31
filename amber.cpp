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


  void Amber::verifyFormat(std::istream& is, const std::string& fmt, const std::string& where) {

    getNextLine(is);
    boost::char_separator<char> sep("()");
    tokenizer tokens(_current_line, sep);
    tokenizer::iterator toks = tokens.begin();

    ++toks;
    if (*toks != fmt)
      throw(FileParseError("Bad format spec while parsing amber file - " + where, _lineno));
  }


  void Amber::parseCharges(std::istream& is) {
    verifyFormat(is, "5E16.8", "charges");

    std::vector<double> charges = readBlock<double>(is, 16);
    if (charges.size() != atoms.size())
      throw(FileParseError("Error parsing charges from amber file", _lineno));
    
    for (uint i=0; i<charges.size(); ++i)
      atoms[i]->charge(charges[i]);
  }


  void Amber::parseMasses(std::istream& is) {
    verifyFormat(is, "5E16.8", "masses");

    std::vector<double> masses = readBlock<double>(is, 16);
    if (masses.size() != atoms.size())
      throw(FileParseError("Error parsing masses from amber file", _lineno));
    
    for (uint i=0; i<masses.size(); ++i)
      atoms[i]->mass(masses[i]);

  }



  void Amber::parseResidueLabels(std::istream& is) {
    verifyFormat(is, "20a4", "residue labels");

    residue_labels = readBlock<std::string>(is, 4);
    if (residue_labels.size() != nres)
      throw(FileParseError("Error parsing residue labels from amber file", _lineno));

  }


  void Amber::parseResiduePointers(std::istream& is) {
    verifyFormat(is, "10I8", "residue pointers");

    residue_pointers = readBlock<uint>(is, 8);
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
    verifyFormat(is, "10I8", "bonds");

  
    std::vector<int> bond_list = readBlock<int>(is, 8);

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
    verifyFormat(is, "10I8", "pointers");

    std::vector<uint> pointers = readBlock<uint>(is, 8);

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


    // The Amber spec says the title format should be 20a4, but at least one
    // variant we've seen uses just a.  So we have to handle both...
    try {
      verifyFormat(is, "20a4", "title");
    }
    catch (FileParseError& e) { _unget = true; verifyFormat(is, "a", "title"); }
    catch (...) { throw; }

    // This reads the actual title...
    getline(is, _title);
  }


  void Amber::parseAtomNames(std::istream& is) {
    verifyFormat(is, "20a4", "atom names");
    
    std::vector<std::string> names = readBlock<std::string>(is, 4);
    if (names.size() != natoms)
      throw(FileParseError("Error parsing atom names", _lineno));
    for (uint i=0; i<names.size(); ++i)
      atoms[i]->name(names[i]);
  }

  void Amber::parseAmoebaRegularBondNumList(std::istream& is) {
    verifyFormat(is, "I8", "amoeba_regular_num_bond_list");
    getNextLine(is);
    std::istringstream iss(_current_line);

    if (! (iss >> std::setw(8) >> _amoeba_regular_bond_num_list))
      throw(FileParseError("Error parsing amoeba_regular_bond_num_list", _lineno));
  }

  void Amber::parseAmoebaRegularBondList(std::istream& is, const uint n) {
    verifyFormat(is, "10I8", "amoeba_regular_bond_list");

  
    std::vector<int> bond_list = readBlock<int>(is, 8);

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

