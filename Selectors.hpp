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






#if !defined(SELECTORS_HPP)
#define SELECTORS_HPP

#include <loos.hpp>
#include <AtomicGroup.hpp>
#include <Kernel.hpp>


//! Predicate for selecting CA atoms
struct CAlphaSelector : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    return(pa->name() == "CA");
  }
};

//! Predicate for selecting backbone
struct BackboneSelector : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    string s = pa->name();
    return(s == "C" || s == "CA" || s == "O" || s == "N");
  }
};


//! Predicate for selecting atoms based on the passed segid string
struct SegidSelector : public AtomSelector {
  explicit SegidSelector(const string s) : str(s) { }
  bool operator()(const pAtom& pa) const {
    return(pa->segid() == str);
  }

  string str;
};

//! Predicate for selecting atoms from a range of resid's
struct ResidRangeSelector : public AtomSelector {
  ResidRangeSelector(const int low, const int high) : _low(low), _high(high) { }
  bool operator()(const pAtom& pa) const {
    return(pa->resid() >= _low && pa->resid() <= _high);
  }

  int _low, _high;
};

//! Predicate for selecting atoms in a specific range of z values
struct ZSliceSelector : public AtomSelector {
    ZSliceSelector(const greal min, const greal max) : _min(min), _max(max) { }
    bool operator()(const pAtom& pa) const {
        greal z = (pa->coords()).z();
        return ( (z>=_min) && (z<_max) );
    }

    greal _min, _max;
};


//! Negates a selection predicate
/*!
  Example:

  \verbatim
  SegidSelector solvsel("SOLV");
  NotSelector notsolvsel(solvsel);
  \endverbatim

  This will select all atoms that are NOT solvent
*/

struct NotSelector : public AtomSelector {
  explicit NotSelector(const AtomSelector& s) : sel(s) { }
  bool operator()(const pAtom& pa) const {
    return(!(sel(pa)));
  }

  const AtomSelector& sel;
};


//! Select hydrogen atoms
struct HydrogenSelector  : public AtomSelector {
  bool operator()(const pAtom& pa) const {
    string n = pa->name();
    char c = n[0];
    return( (c == 'H') && (pa->mass() < 1.1) );
  }
};

//! Select non-hydrogen atoms
struct HeavyAtomSelector : public AtomSelector {
    HydrogenSelector hsel;
    NotSelector not_heavy;
    HeavyAtomSelector() : not_heavy(hsel) { }
    bool operator()(const pAtom& pa) const {
        return (not_heavy(pa));
    }
};


//! Combines two selectors with a logical "and"
/** Example:
 *  \verbatim
 *  SegidSelector prot("PROT");
 *  MainChainSelector main_chain;
 *  AndSelector main_chain_protein(main_chain, prot);
 *  \endverbatim
 *
 *  The main_chain_protein selector will select for all atoms that are
 *  both main chain and have a segid of "PROT".
*/

struct AndSelector : public AtomSelector {
  AndSelector(const AtomSelector& x, const AtomSelector& y) : lhs(x), rhs(y) { }
  bool operator()(const pAtom& pa) const {
    return(lhs(pa) && rhs(pa));
  }

  const AtomSelector& lhs;
  const AtomSelector& rhs;
};



//! Combines two selectors with a logical "or"
/** Example:
 *  \verbatim
 *  SegidSelector prot("PROT");
 *  SegidSelector heme("HEME");
 *  OrSelector prot_with_heme(prot, heme);
 *  \endverbatim
 *
 *  This selector will pick any atom that has a segid of either "PROT"
 *  or "HEME".
*/

struct OrSelector : public AtomSelector {
  OrSelector(const AtomSelector& x, const AtomSelector& y) : lhs(x), rhs(y) { }
  bool operator()(const pAtom& pa) const {
    return(lhs(pa) || rhs(pa));
  }

  const AtomSelector& lhs;
  const AtomSelector& rhs;
};







//! Predicate for selecting solvent based on common solvent SEGIDs
struct SolventSelector : public AtomSelector {
  SolventSelector() : s1("SOLV"), s2("BULK"), osel(s1, s2) {}
  bool operator()(const pAtom& pa) const {
    return(osel(pa));
  }

  SegidSelector s1, s2;
  OrSelector osel;
};




//! Selection predicate that executes a compiled Kernel
/**
 * This predicate takes a compiled Kernel and executes it once for each
 * Atom.  This is primarily for use in conjunction with the Parser for
 * handling selections based on user input.
 *
 * Example:
 * \verbatim
 * Parser parsed(selection_string);
 * KernelSelector sel(parsed.kernel());
 * \endverbatim
 *
 */
class KernelSelector : public AtomSelector {
public:
  explicit KernelSelector(loos::Kernel& k) : krnl(k) { }

  bool operator()(const pAtom& pa) const {
    krnl.execute(pa);
    if (krnl.stack().size() != 1) {
      throw(runtime_error("Execution error - unexpected values on stack"));
    }

    loos::Value results = krnl.stack().pop();
    if (results.type != loos::Value::INT)
      throw(runtime_error("Execution error - unexpected value on top of stack"));

    return(results.itg);
  }


private:
  loos::Kernel& krnl;

};




#endif

