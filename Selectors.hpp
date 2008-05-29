/*
  Selectors.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  A selector library (of sorts...)
*/




#if !defined(SELECTORS_HPP)
#define SELECTORS_HPP


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
  SegidSelector(const string s) : str(s) { }
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
  NotSelector(const AtomSelector& s) : sel(s) { }
  bool operator()(const pAtom& pa) const {
    return(!(sel(pa)));
  }

  const AtomSelector& sel;
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
  KernelSelector(loos::Kernel& k) : krnl(k) { }

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

