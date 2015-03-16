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






#if !defined(LOOS_SELECTORS_HPP)
#define LOOS_SELECTORS_HPP

#include <loos_defs.hpp>
#include <AtomicGroup.hpp>
#include <Kernel.hpp>


namespace loos {

  //! Predicate for selecting CA atoms
  struct CAlphaSelector : public AtomSelector {
    bool operator()(const pAtom&) const;
  };

  //! Predicate for selecting backbone
  class BackboneSelector : public AtomSelector {
    static const uint nresnames = 35;
    static std::string residue_names[nresnames];

    static const uint natomnames = 33;
    static std::string atom_names[natomnames];

  public:
    bool operator()(const pAtom&) const;
  };


  //! Predicate for selecting atoms based on the passed segid string
  struct SegidSelector : public AtomSelector {
    explicit SegidSelector(const std::string s) : str(s) { }
    bool operator()(const pAtom&) const;

    std::string str;
  };


  //! Predicate for selecting atoms based on explicit name matching
  struct AtomNameSelector : public AtomSelector {
    explicit AtomNameSelector(const std::string& s) : str(s) { }
    bool operator()(const pAtom&) const;

    std::string str;
  };


  //! Predicate for selecting atoms from a range of resid's
  struct ResidRangeSelector : public AtomSelector {
    ResidRangeSelector(const int low, const int high) : _low(low), _high(high) { }
    bool operator()(const pAtom&) const;

    int _low, _high;
  };

  //! Predicate for selecting atoms in a specific range of z values
  struct ZSliceSelector : public AtomSelector {
    ZSliceSelector(const greal min, const greal max) : _min(min), _max(max) { }
    bool operator()(const pAtom& pa) const;

    greal _min, _max;
  };


  //! Negates a selection predicate
  /*!
    Example:

    \code
    SegidSelector solvsel("SOLV");
    NotSelector notsolvsel(solvsel);
    \endcode

    This will select all atoms that are NOT solvent
  */

  struct NotSelector : public AtomSelector {
    explicit NotSelector(const AtomSelector& s) : sel(s) { }
    bool operator()(const pAtom&) const;

    const AtomSelector& sel;
  };


  //! Select hydrogen atoms
  struct HydrogenSelector  : public AtomSelector {
    bool operator()(const pAtom&) const;
  };

  //! Select non-hydrogen atoms
  struct HeavyAtomSelector : public AtomSelector {
    HydrogenSelector hsel;
    NotSelector not_heavy;
    HeavyAtomSelector() : not_heavy(hsel) { }
    bool operator()(const pAtom& pa) const;
  };


  //! Combines two selectors with a logical "and"
  /** Example:
   *  \code
   *  SegidSelector prot("PROT");
   *  MainChainSelector main_chain;
   *  AndSelector main_chain_protein(main_chain, prot);
   *  \endcode
   *
   *  The main_chain_protein selector will select for all atoms that are
   *  both main chain and have a segid of "PROT".
   */

  struct AndSelector : public AtomSelector {
    AndSelector(const AtomSelector& x, const AtomSelector& y) : lhs(x), rhs(y) { }
    bool operator()(const pAtom& pa) const;

    const AtomSelector& lhs;
    const AtomSelector& rhs;
  };



  //! Combines two selectors with a logical "or"
  /** Example:
   *  \code
   *  SegidSelector prot("PROT");
   *  SegidSelector heme("HEME");
   *  OrSelector prot_with_heme(prot, heme);
   *  \endcode
   *
   *  This selector will pick any atom that has a segid of either "PROT"
   *  or "HEME".
   */

  struct OrSelector : public AtomSelector {
    OrSelector(const AtomSelector& x, const AtomSelector& y) : lhs(x), rhs(y) { }
    bool operator()(const pAtom& pa) const;

    const AtomSelector& lhs;
    const AtomSelector& rhs;
  };







  //! Predicate for selecting solvent based on common solvent SEGIDs
  struct SolventSelector : public AtomSelector {
    SolventSelector() : s1("SOLV"), s2("BULK"), osel(s1, s2) {}
    bool operator()(const pAtom& pa) const;

    SegidSelector s1, s2;
    OrSelector osel;
  };



  //! Select only heavy solvent atoms
  struct HeavySolventSelector : public AtomSelector {
    HeavySolventSelector() : sel(s1, s2) { }
    bool operator()(const pAtom& pa) const;

    SolventSelector s1;
    HeavyAtomSelector s2;
    AndSelector sel;
  };




  //! Selection predicate that executes a compiled Kernel
  /**
   * This predicate takes a compiled Kernel and executes it once for each
   * Atom.  This is primarily for use in conjunction with the Parser for
   * handling selections based on user input.
   *
   * Example:
   * \code
   * Parser parsed(selection_string);
   * KernelSelector sel(parsed.kernel());
   * \endcode
   *
   */
  class KernelSelector : public AtomSelector {
  public:
    explicit KernelSelector(Kernel& k) : krnl(k) { }

    bool operator()(const pAtom& pa) const;

  private:
    Kernel& krnl;

  };


}

#endif

