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



#include <ios>
#include <sstream>
#include <iomanip>

#include <assert.h>
#include <vector>
#include <list>
#include <algorithm>
#include <unordered_set>
#include <functional>

#include <boost/random.hpp>

#include <AtomicGroup.hpp>
#include <AtomicNumberDeducer.hpp>
#include <Selectors.hpp>
#include <Parser.hpp>

#include <boost/unordered_map.hpp>

namespace loos {


  typedef boost::unordered_map<int,int>    IMap;



  const double AtomicGroup::superposition_zero_singular_value  =  1e-10;


  AtomicGroup* AtomicGroup::clone(void) const {
    return(new AtomicGroup(*this));
  }


  AtomicGroup AtomicGroup::copy(void) const {
    const_iterator i;
    AtomicGroup res;

    for (i = atoms.begin(); i != atoms.end(); i++) {
      pAtom pa(new Atom(**i));
      res.append(pa);
    }
    res._sorted = _sorted;
    res.box = box.copy();

    return(res);
  }


  // Internal: verify the index into the atom array...
  int AtomicGroup::rangeCheck(int i) const {
    if (i < 0)
      i = atoms.size() + i;
    if ((unsigned int)i >= atoms.size())
      throw(std::out_of_range("Bad index for an atom"));

    return(i);
  }


  // Returns the ith atom...
  pAtom AtomicGroup::getAtom(const int i) const {
    int j = rangeCheck(i);

    return(atoms[j]);
  }

  // Should these invalidate sort status?
  pAtom& AtomicGroup::operator[](const int i) {
    int j = rangeCheck(i);
    return(atoms[j]);
  }

  // For const objects...
  const pAtom& AtomicGroup::operator[](const int i) const {
    int j = rangeCheck(i);
    return(atoms[j]);
  }

  // Internal: removes an atom from this group based on the address of the shared pointer...
  void AtomicGroup::deleteAtom(pAtom pa) {
    std::vector<pAtom>::iterator iter;

    iter = find(atoms.begin(), atoms.end(), pa);
    if (iter == atoms.end())
      throw(LOOSError(*pa, "Attempting to delete an atom that is not in the passed AtomicGroup"));

    atoms.erase(iter);
    _sorted = false;
  }


  // Append each atom from the passed vector onto this group...
  AtomicGroup& AtomicGroup::append(std::vector<pAtom> pas) {
    std::vector<pAtom>::iterator i;

    for (i=pas.begin(); i != pas.end(); i++)
      atoms.push_back(*i);

    _sorted = false;
    return(*this);
  }


  // Append all atoms from the passed group onto this one
  AtomicGroup& AtomicGroup::append(const AtomicGroup& grp) {
    std::vector<pAtom>::const_iterator i;

    for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
      addAtom(*i);

    _sorted = false;
    return(*this);
  }


  // Remove all atoms in the passed vector
  AtomicGroup& AtomicGroup::remove(std::vector<pAtom> pas) {
    std::vector<pAtom>::iterator i;

    for (i=pas.begin(); i != pas.end(); i++)
      deleteAtom(*i);

    _sorted = false;
    return(*this);
  }


  // Removes all atoms contained in the passed group from this one...
  AtomicGroup& AtomicGroup::remove(const AtomicGroup& grp) {


    if (&grp == this)
      atoms.clear();      // Assume caller meant to clean out AtomicGroup
    else {
      std::vector<pAtom>::const_iterator i;

      for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
        deleteAtom(*i);

      _sorted = false;
      return(*this);
    }

    return(*this);
  }


  // Concatenation operations...
  AtomicGroup& AtomicGroup::operator+=(const AtomicGroup& rhs) {
    append(rhs);
    return(*this);
  }

  AtomicGroup& AtomicGroup::operator+=(const pAtom& rhs) {
    atoms.push_back(rhs);
    _sorted = false;
    return(*this);
  }

  AtomicGroup AtomicGroup::operator+(const AtomicGroup& rhs) {
    AtomicGroup res(*this);
    res += rhs;
    return(res);
  }

  AtomicGroup AtomicGroup::operator+(const pAtom& rhs) {
    AtomicGroup res(*this);
    res += rhs;
    return(res);
  }


  bool AtomicGroup::allHaveProperty(const Atom::bits& property) const {
    for (const_iterator ci = atoms.begin(); ci != atoms.end(); ++ci)
      if (!(*ci)->checkProperty(property))
        return(false);

    return(true);
  }

  bool AtomicGroup::anyHaveProperty(const Atom::bits& property) const {
    for (const_iterator ci = atoms.begin(); ci != atoms.end(); ++ci)
      if ((*ci)->checkProperty(property))
        return(true);

    return(false);
  }

  bool AtomicGroup::hasBonds(void) const {
    const_iterator ci;

    for (ci = atoms.begin(); ci != atoms.end(); ++ci)
      if ((*ci)->checkProperty(Atom::bondsbit))
        return(true);

    return(false);
  }

  // returns a vector of pairs of atom IDs corresponding to bonded atoms
  std::vector<std::pair<int, int>> AtomicGroup::getBondsIDs() const {
    std::unordered_set<std::pair<int, int>, unordered_pair_hash<int>, unordered_pair_eq<int>> bond_set;
    const_iterator ci;

    for (ci = atoms.begin(); ci != atoms.end(); ++ci){
      std::vector<int> bonds = (*ci)->getBonds();
      int id = (*ci)->id();
      for(auto b : bonds){
        bond_set.emplace(std::make_pair(id, b));
      }
    }
    std::vector<std::pair<int, int>> bond_list(bond_set.begin(), bond_set.end());
    return bond_list;
  }
  
  // returns a vector of pairs of atom IDs corresponding to bonded atoms
  std::vector<AtomicGroup> AtomicGroup::getBondsAGs() const {
    std::unordered_set<std::pair<int, int>, unordered_pair_hash<int>, unordered_pair_eq<int>> bond_set;
    const_iterator ci;

    for (ci = atoms.begin(); ci != atoms.end(); ++ci){
      std::vector<int> bonds = (*ci)->getBonds();
      int id = (*ci)->id();
      for(auto b : bonds){
        bond_set.emplace(std::make_pair(id, b));
      }
    }
    std::vector<AtomicGroup> bond_list;
    for (const auto& bond : bond_set) {
      bond_list.emplace_back(groupFromID(bond));
    }
    return bond_list;
  }

  void AtomicGroup::clearBonds(void) {
    const_iterator ci;

    for (ci = atoms.begin(); ci != atoms.end(); ++ci)
      (*ci)->clearBonds();
  }


  bool AtomicGroup::hasCoords(void) const {
    const_iterator ci;

    for (ci = atoms.begin(); ci != atoms.end(); ++ci)
      if (!((*ci)->checkProperty(Atom::coordsbit)))
        return(false);

    return(true);
  }


  void AtomicGroup::sort(void) {
    CmpById comp;

    if (! _sorted)
      std::sort(atoms.begin(), atoms.end(), comp);

    _sorted = true;
  }



  // Internal: calculates the start and stop iterators given offset and len args
  // as in PERL's substr()...

  boost::tuple<AtomicGroup::iterator, AtomicGroup::iterator> AtomicGroup::calcSubsetIterators(const int offset, const int len) {
    unsigned int a, b;

    if (offset < 0) {
      b = atoms.size() + offset + 1;
      a = (len == 0) ? 0 : b - len;
    } else {
      a = offset;
      b = (len == 0) ? atoms.size() : a + len;
    }

    if (b-a >= atoms.size())
      throw(std::range_error("Indices out of bounds for subsetting"));

    boost::tuple<iterator, iterator> res(atoms.begin() + a, atoms.begin() + b);

    return(res);
  }



  AtomicGroup AtomicGroup::subset(const int offset, const int len) {
    AtomicGroup res;

    boost::tuple<iterator, iterator> iters = calcSubsetIterators(offset, len);
    res.atoms.insert(res.atoms.begin(), boost::get<0>(iters), boost::get<1>(iters));

    res.box = box;
    return(res);
  }


  AtomicGroup AtomicGroup::excise(const int offset, const int len) {
    AtomicGroup res;

    boost::tuple<iterator, iterator> iters = calcSubsetIterators(offset, len);

    res.atoms.insert(res.atoms.begin(), boost::get<0>(iters), boost::get<1>(iters));
    atoms.erase(boost::get<0>(iters), boost::get<1>(iters));

    _sorted = false;

    res.box = box;
    return(res);
  }

  // Select atoms from the current group, adding them to a new group
  // based on the sel functor/predicate...

  AtomicGroup AtomicGroup::select(const AtomSelector& sel) const {
    AtomicGroup res;

    std::vector<pAtom>::const_iterator i;
    for (i=atoms.begin(); i != atoms.end(); i++)
      if (sel(*i))
        res.addAtom(*i);

    res.box = box;
    return(res);
  }


  // Split up a group into a vector of groups based on unique segids...
  std::vector<AtomicGroup> AtomicGroup::splitByUniqueSegid(void) const {
    const_iterator i;
    std::list<std::string> segids;

    for (i = atoms.begin(); i != atoms.end(); i++)
      if (find(segids.begin(), segids.end(), (*i)->segid()) == segids.end())
        segids.push_front((*i)->segid());

    // Reverse the order so that the chunks are in the same order they appear
    // in the file...
    uint n = segids.size();
    std::vector<AtomicGroup> results(n);
    uint j = n-1;
    for (std::list<std::string>::const_iterator i = segids.begin(); i != segids.end(); ++i) {
      SegidSelector segid(*i);
      results[j--] = select(segid);
    }

    std::vector<AtomicGroup>::iterator g;
    for (g=results.begin(); g!=results.end(); g++) {
      g->box = box;
    }

    return(results);
  }

  std::map<std::string, AtomicGroup> AtomicGroup::splitByName(void) const {
    const_iterator i;
    std::map<std::string, AtomicGroup> groups;

    // Loop over atoms, adding them to groups based on their name, creating the
    // map entry for each new name as we find it
    std::map<std::string, AtomicGroup>::iterator g;
    for (i = atoms.begin(); i != atoms.end(); ++i) {
        g = groups.find((*i)->name());
        if (g == groups.end()) { // not found, need to create a new AG
            AtomicGroup ag;
            ag.append(*i);
            ag.box = box; // copy the current groups periodic box
            groups[(*i)->name()] = ag;
        }  else {              // found group for that atom name,
                               // so add the atom to it
            g->second.append(*i);
        }

    }

    return(groups);
  }

  /** The idea is that we iterate over the list of contained atoms.  For
   * each atom, we recurse through the list of bonded atoms.  Each time
   * we visit an atom, we mark it as having been seen via the hash_set.
   * In the recursive function, every time we find a new unseen atom, we
   * append it to the current group and mark it as seen, then recurse
   * through all of its bonded atoms.
   *
   * Since splitByMolecule() recurses through the connectivity list,
   * the ordering of the atoms returned will not necessarily be the
   * same as the input group.  We therefore sort each group returned
   * to hopefully maintain correct relative ordering.
   *
   * If we find a bond that goes to an atom that does not exist in the
   * current group, a std::runtime_error is thrown.
   *
   * Note that for efficiency of atom lookup, the group gets sorted.
   * This is why the public function splitByMolecule() makes a copy of
   * itself and calls sortingSplitByMolecule() on that...
   */

  std::vector<AtomicGroup> AtomicGroup::sortingSplitByMolecule(void) {
    std::vector<AtomicGroup> molecules;

    // If no connectivity, just return the entire group...
    if (!hasBonds()) {
      sort();
      molecules.push_back(*this);
    } else {
      HashInt seen;                      // Track what atoms we've already processed
      AtomicGroup current;               // The molecule we're currently building...
      AtomicGroup working(*this);        // Copy, so we can sort without mucking up original order
      working.sort();

      int n = size();
      for (int i=0; i<n; i++) {
        HashInt::iterator it = seen.find(working.atoms[i]->id());
        if (it != seen.end())
          continue;

        walkBonds(current, seen, working, working.atoms[i]);
        if (current.size() != 0) {       // Just in case...
          current.sort();
          molecules.push_back(current);
          current = AtomicGroup();
        }

      }
    }
    

    // copy the box over
    std::vector<AtomicGroup>::iterator m;
    for (m=molecules.begin(); m!=molecules.end(); m++) {
      m->box = box;
    }
    return(molecules);
  }

  std::vector<AtomicGroup> AtomicGroup::sortingSplitByMolecule(const std::string& selection) {
    std::vector<AtomicGroup> molecules;
    Parser parser(selection);
    KernelSelector parsed_sel(parser.kernel());


    // If no connectivity, just return the entire group...
    if (!hasBonds()) {
      sort();
      molecules.push_back(this->select(parsed_sel));
    } else {
      HashInt seen;                      // Track what atoms we've already processed
      AtomicGroup current;               // The molecule we're currently building...
      AtomicGroup working(*this);        // Copy, so we can sort without mucking up original order
      working.sort();

      int n = size();
      for (int i=0; i<n; i++) {
        HashInt::iterator it = seen.find(working.atoms[i]->id());
        if (it != seen.end())
          continue;

        walkBonds(current, seen, working, working.atoms[i]);
        if (current.size() != 0) {       // Just in case...
          current.sort();
          molecules.push_back(current.select(parsed_sel));
          current = AtomicGroup();
        }

      }
    }
    

    // copy the box over
    std::vector<AtomicGroup>::iterator m;
    for (m=molecules.begin(); m!=molecules.end(); m++) {
      m->box = box;
    }
    return(molecules);
  } 


  void AtomicGroup::walkBonds(AtomicGroup& current, HashInt& seen, AtomicGroup& working, pAtom& moi) {
    int myid = moi->id();
    HashInt::iterator it = seen.find(myid);

    // If we've touched this atom before, stop recursing and return.
    if (it != seen.end())
      return;

    // This is a novel atom, so append it into the group we're currently building.
    seen.insert(myid);
    current.addAtom(moi);

    // Just in case it's a solo-atom...  This probably should indicate
    // some kind of error...?
    if (!(moi->hasBonds()))
      return;

    // Now find atoms that are bound to the current atom and recurse
    // through them...

    std::vector<int> bonds = moi->getBonds();
    std::vector<int>::const_iterator citer;
    for (citer = bonds.begin(); citer != bonds.end(); citer++) {
      pAtom toi = working.findById(*citer);
      if (toi != 0)
        walkBonds(current, seen, working, toi);
    }
  }


  /**
   * Splits an AtomicGroup into individual residues.  The residue
   * boundary is marked by either a change in the resid or in the
   * segid.
   */
  std::vector<AtomicGroup> AtomicGroup::splitByResidue(void) const {
    std::vector<AtomicGroup> residues;

    int curr_resid = atoms[0]->resid();
    std::string curr_segid = atoms[0]->segid();

    AtomicGroup residue;
    AtomicGroup::const_iterator ci;
    for (ci = atoms.begin(); ci != atoms.end(); ++ci) {
      if (curr_resid != (*ci)->resid() || (*ci)->segid() != curr_segid) {
        residues.push_back(residue);
        residue = AtomicGroup();
        curr_resid = (*ci)->resid();
        curr_segid = (*ci)->segid();
      }
      residue.append(*ci);
    }

    if (residue.size() != 0)
      residues.push_back(residue);


    // Copy the box information
    for (std::vector<AtomicGroup>::iterator i = residues.begin(); i != residues.end(); ++i)
      i->box = box;

    return(residues);
  }


  // Find an atom based on atomid
  // Returns 0 (null shared_ptr) if not found...
  pAtom AtomicGroup::findById_binarySearch(const int id) const {
    int bottom = 0, top = size()-1, middle;

    while (top > bottom) {
      middle = bottom + (top - bottom) / 2;
      if (atoms[middle]->id() < id)
        bottom = middle + 1;
      else
        top = middle;
    }

    if (atoms[bottom]->id() == id)
      return(atoms[bottom]);

    return(pAtom());
  }


  pAtom AtomicGroup::findById_linearSearch(const int id) const {
    for (AtomicGroup::const_iterator i = begin(); i != end(); ++i)
      if ((*i)->id() == id)
        return(*i);

    return(pAtom());
  }


  pAtom AtomicGroup::findById(const int id) const {
    if (sorted())
      return(findById_binarySearch(id));

    return(findById_linearSearch(id));
  }


  /**
   * Note: when calling this, you'll want to make sure you use the
   * outermost group (eg the psf or pdb you used to create things, rather than
   * using a subselection, unless you're sure the subsection contains these
   * atoms as well.  The main use of this routine is to create a group of atoms
   * bound to another atom.
   *
   * Any missing atoms are ignored...  This is in contrast with the previous
   * behavior where missing atoms would throw an exception
   */

  AtomicGroup AtomicGroup::groupFromID(const std::vector<int> &id_list) const {
    AtomicGroup result;

    result.box = box;

    for (unsigned int i=0; i<id_list.size(); i++) {
      pAtom pa = findById(id_list[i]);
      if (pa)
        result.addAtom(pa);
    }
    return(result);
  }
  
  // Note this specialization can also produce an empty AG, or an AG with one element.
  AtomicGroup AtomicGroup::groupFromID(const std::pair<int, int> &id_pair) const {
    AtomicGroup result;

    result.box = box;

    pAtom pa1 = findById(id_pair.first);
    if (pa1)
      result.addAtom(pa1);
    pAtom pa2 = findById(id_pair.second);
    if (pa2)
      result.addAtom(pa2);
    return(result);
  }

  // Get all atoms associated with the residue that contains the
  // passed atom...  The returned atoms will not be in order.  If
  // you want that, then explicitly sort the group.
  AtomicGroup AtomicGroup::getResidue(pAtom res) {
    iterator i;
    AtomicGroup result;

    result.box = box;
    i = find(atoms.begin(), atoms.end(), res);
    if (i == atoms.end())
      return(result);

    iterator j = i;

    while (j >= atoms.begin()) {
      if ((*j)->resid() == res->resid() && (*j)->segid() == res->segid())
        result.addAtom(*j);
      else
        break;
      --j;
    }

    j = i+1;
    while (j < atoms.end()) {
      if ((*j)->resid() == res->resid() && (*j)->segid() == res->segid())
        result.addAtom(*j);
      else
        break;

      ++j;
    }

    return(result);
  }


  namespace {
    // renumber the atomids of the group...
    void renumberWithoutBonds(AtomicGroup& grp, const int start, const int stride) {
      AtomicGroup::iterator i;
      int id = start;

      for (i=grp.begin(); i != grp.end(); i++, id += stride)
        (*i)->id(id);
    }


    void renumberWithBonds(AtomicGroup& grp, const int start, const int stride) {
      IMap bond_map;

      int id = start;
      for (AtomicGroup::iterator i = grp.begin(); i != grp.end(); ++i, id += stride) {
        bond_map[(*i)->id()] = id;
        (*i)->id(id);
      }

      // Now transform bond list
      for (AtomicGroup::iterator i = grp.begin(); i != grp.end(); ++i, id += stride) {
        if ((*i)->hasBonds()) {
          std::vector<int> bonds = (*i)->getBonds();
          for (std::vector<int>::iterator j = bonds.begin(); j != bonds.end(); ++j) {
            IMap::iterator k = bond_map.find(*j);
            if (k != bond_map.end())
              *j = k->second;
          }
          (*i)->setBonds(bonds);
        }
      }
    }
  };

  void AtomicGroup::renumber(const int start, const int stride) {

    if (hasBonds())
      renumberWithBonds(*this, start, stride);
    else
      renumberWithoutBonds(*this, start, stride);
  }

  // Get the min and max atomid's...
  int AtomicGroup::minId(void) const {
    const_iterator i;

    if (atoms.size() == 0)
      return(-1);

    int min = atoms[0]->id();
    for (i = atoms.begin()+1; i != atoms.end(); i++)
      if ((*i)->id() < min)
        min = (*i)->id();

    return(min);
  }


  int AtomicGroup::maxId(void) const {
    const_iterator i;

    if (atoms.size() == 0)
      return(-1);
    int max = atoms[0]->id();
    for (i=atoms.begin()+1; i != atoms.end(); i++)
      if ((*i)->id() > max)
        max = (*i)->id();

    return(max);
  }


  int AtomicGroup::minResid(void) const {
    const_iterator i;

    if (atoms.size() == 0)
      return(-1);

    int min = atoms[0]->resid();
    for (i = atoms.begin()+1; i != atoms.end(); i++)
      if ((*i)->resid() < min)
        min = (*i)->resid();

    return(min);
  }

  int AtomicGroup::maxResid(void) const {
    const_iterator i;

    if (atoms.size() == 0)
      return(-1);

    int max = atoms[0]->resid();
    for (i = atoms.begin()+1; i != atoms.end(); i++)
      if ((*i)->resid() > max)
        max = (*i)->resid();

    return(max);
  }


  // Count the number of higher structural elements...
  int AtomicGroup::numberOfResidues(void) const {

    if (atoms.size() == 0)
      return(0);

    const_iterator i;
    int n = 1;
    int curr_resid = atoms[0]->resid();
    std::string curr_segid = atoms[0]->segid();

    for (i=atoms.begin()+1; i !=atoms.end(); i++)
      if (((*i)->resid() != curr_resid) ||
          ((*i)->segid() != curr_segid)) {
        ++n;
        curr_resid = (*i)->resid();
        curr_segid = (*i)->segid();
      }

    return(n);
  }


  int AtomicGroup::numberOfSegids(void) const {

    if (atoms.size() == 0)
      return(0);

    const_iterator i;
    int n = 1;
    std::string curr_segid = atoms[0]->segid();

    for (i=atoms.begin()+1; i !=atoms.end(); i++)
      if ((*i)->segid() != curr_segid) {
        ++n;
        curr_segid = (*i)->segid();
      }

    return(n);
  }



  AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs) {
    AtomicGroup res;
    res.append(lhs);
    res.append(rhs);
    return(res);
  }


  AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs) {
    AtomicGroup res(rhs);
    res += lhs;
    return(res);
  }



  bool AtomicGroup::operator==(AtomicGroup& rhs) {

    if (size() != rhs.size())
      return(false);

    if (this == &rhs)
      return(true);

    sort();
    rhs.sort();

    int n = size();
    for (int i = 0; i < n; i++)
      if (atoms[i] != rhs.atoms[i])
        return(false);

    return(true);
  }


  bool AtomicGroup::operator==(const AtomicGroup& rhs) const {
    if (size() != rhs.size())
      return(false);

    if (this == &rhs)
      return(true);

    const std::vector<pAtom> *lp;
    const std::vector<pAtom> *rp;
    std::vector<pAtom> lhs_atoms, rhs_atoms;
    CmpById comp;
    if (!sorted()) {
      lhs_atoms = atoms;
      std::sort(lhs_atoms.begin(), lhs_atoms.end(), comp);
      lp = &lhs_atoms;
    } else
      lp = &atoms;

    if (!rhs.sorted()) {
      rhs_atoms = rhs.atoms;
      std::sort(rhs_atoms.begin(), rhs_atoms.end(), comp);
      rp = &rhs_atoms;
    } else
      rp = &rhs.atoms;

    int n = size();
    for (int i = 0; i<n; i++)
      if ((*lp)[i] != (*rp)[i])
        return(false);

    return(true);
  }


  // this function takes a bond pair offset for the whole AG.
  // Returns the unit vector projection across all such bond pairs.
  const greal AtomicGroup::ocf(uint offset){
    greal part_ocf = 0;
    GCoord bv1, bv2; 
    for (auto i = 0; i < atoms.size() - offset - 1; i++){
      // compute bond vector between atom and its next neighbor.
      bv1 = atoms[i]->coords() - atoms[i+1]->coords();
      // compute bond vector between offset atom and its next neighbor.
      bv2 = atoms[i + offset]->coords() - atoms[i + offset + 1]->coords();
      // accumulate dot product of unit vectors 
      part_ocf += bv1.uvdot(bv2);
    }
    return(part_ocf / (atoms.size() - offset - 1));
  }

  void AtomicGroup::reimage() {
    if (!(isPeriodic()))
      throw(LOOSError("trying to reimage a non-periodic group"));
    GCoord com = centroid();
    GCoord reimaged = com;
    reimaged.reimage(periodicBox());
    GCoord trans = reimaged - com;
    const_iterator a;
    for (a=atoms.begin(); a!=atoms.end(); a++) {
      (*a)->coords() += trans;
    }
  }

  void AtomicGroup::reimageByAtom () {
    if (!(isPeriodic()))
      throw(LOOSError("trying to reimage a non-periodic group"));
    const_iterator a;
    GCoord box = periodicBox();
    for (a=atoms.begin(); a!=atoms.end(); a++) {
      (*a)->coords().reimage(box);
    }
  }

  /** Works by translating the system so one atom is in the center of the
   *  box, reimaging by atom (so now the group is all in the middle of the box),
   *  and then translating back.
   *
   *  If you don't want to give it a reference atom, call the version
   *  that takes no argument; it uses the first atom in the
   *  AtomicGroup.
   *
   */
  void AtomicGroup::mergeImage(pAtom &p ) {
      GCoord ref = p->coords();

      translate(-ref);
      reimageByAtom();
      translate(ref);
  }

  /** Does the same as the other mergeImage, only using the first atom in the
   *  AtomicGroup as the reference atom.
   */
  void AtomicGroup::mergeImage() {

    // Ignore reimaging if the group is empty...
    if (!empty()) {
      mergeImage(atoms[0]);
    }
  }






  /**
   * The connectivity list is searched for each atom and if a bond is
   * not found in the current group, then it is removed from the
   * bond-list for that atom.  Note that this means that any
   * AtomicGroup sharing the atom in question will also now have the
   * modified bond list.  It's therefore recommended that this
   * function be called on a copy (AtomicGroup::copy()).  Also note
   * that FindById() does not implicitly sort the atoms for more
   * efficient searching.  You may want to call AtomicGroup::sort()
   * prior to AtomicGroup::pruneBonds() if the exact atom order does
   * not matter.
   */

  void AtomicGroup::pruneBonds() {

    for (AtomicGroup::iterator j = begin(); j != end(); ++j)
      if ((*j)->hasBonds()) {
        std::vector<int> bonds = (*j)->getBonds();
        std::vector<int> pruned_bonds;
        for (std::vector<int>::const_iterator i = bonds.begin(); i != bonds.end(); ++i)
          if (findById(*i) != 0)
            pruned_bonds.push_back(*i);
        (*j)->setBonds(pruned_bonds);
      }
  }


  /**
   * The Atom index is the original ordering of atoms from whatever
   * file format the model came from.  This is used as an index
   * into each frame of the trajectory for corresponding atom
   * properties (such as coordinates).  If an AtomicGroup is a subset,
   * then it may be necessary to reset the indices when working with
   * a subsetted trajectory as well.  This function will reset the
   * atom indices to be sequential, beginning with 0.
   */
  void AtomicGroup::resetAtomIndices()
  {
    for (uint i=0; i<size(); ++i)
      atoms[i]->index(i);
  }


  /**
   * If an atom has a mass, then this is used to look up it's atomic
   * number.  Note that LOOS only has the first 96 elements in its
   * tables.  If a mass is not found in the LOOS table, then the
   * atomic number is not modified (or set), otherwise any existing
   * atomic number is overwritten.
   */
  uint AtomicGroup::deduceAtomicNumberFromMass(const double tol) {
    uint n = 0;

    for (AtomicGroup::iterator i = begin(); i != end(); ++i)
      if ((*i)->checkProperty(Atom::massbit)) {
        uint an = loos::deduceAtomicNumberFromMass((*i)->mass(), tol);
        if (an) {
          (*i)->atomic_number(an);
          ++n;
        }
      }

    return(n);
  }


  void AtomicGroup::copyCoordinatesWithIndex(const std::vector<GCoord> &coords) {
    if (! atoms.empty())
      if (! atoms[0]->checkProperty(Atom::indexbit))
        throw(LOOSError(*(atoms[0]), "Cannot use copyCoordinatesWithIndex() on an atom that does not have an index set"));

    for (uint i=0; i<atoms.size(); ++i)
    {
      uint index = atoms[i]->index();
      atoms[i]->coords( coords.at(index) );
    }
  }

  void AtomicGroup::copyVelocitiesWithIndex(const std::vector<GCoord> &velocities) {
    if (! atoms.empty())
      if (! atoms[0]->checkProperty(Atom::indexbit))
        throw(LOOSError(*(atoms[0]), "Cannot use copyVelocitiesWithIndex() on an atom that does not have an index set"));

    for (uint i=0; i<atoms.size(); ++i)
    {
      uint index = atoms[i]->index();
      atoms[i]->velocities( velocities.at(index) );
    }
  }


  void AtomicGroup::copyCoordinatesFrom(const AtomicGroup& g, const uint offset, const uint length) {
    uint n = (length == 0 || length > g.size()) ? g.size() : length;


    for (uint i=0; i<n && i+offset<atoms.size(); ++i)
      atoms[i+offset]->coords(g[i]->coords());
  }


  std::vector<uint> AtomicGroup::atomOrderMapFrom(const AtomicGroup& g) {
    if (g.size() != size())
      throw(LOOSError("Cannot map atom order between groups of different sizes"));

    std::vector<AtomicGroup> other = g.splitByResidue();
    std::vector<AtomicGroup> self = splitByResidue();

    std::vector<uint> order;
    uint idx = 0;
    for (uint k=0; k<self.size(); ++k) {
      if (self[k][0]->resname() != other[k][0]->resname())
        throw(LOOSError("Mismatched residues while trying to map atom order between groups"));
      for (uint j=0; j<self[k].size(); ++j) {
        uint i;
        for (i=0; i<other[k].size(); ++i)
          if (self[k][j]->name() == other[k][i]->name())
            break;
        if (i < other[k].size())
          order.push_back(idx + i);
        else
          throw(LOOSError(*(self[k][j]), "Cannot find match while constructing atom order map"));
      }
      idx += other[k].size();
    }

    return(order);
  }

  void AtomicGroup::copyMappedCoordinatesFrom(const AtomicGroup& g, const std::vector<uint>& map) {
    if (g.size() != map.size())
      throw(LOOSError("Atom order map is of incorrect size to copy coordinates"));
    if (g.size() != size())
      throw(LOOSError("Cannot copy coordinates (with atom ordering) from an AtomicGroup of a different size"));

    for (uint i=0; i<map.size(); ++i)
      atoms[i]->coords(g[map[i]]->coords());
  }

  void AtomicGroup::copyMappedCoordinatesFrom(const AtomicGroup& g) {
    std::vector<uint> map = atomOrderMapFrom(g);
    copyMappedCoordinatesFrom(g, map);
  }


  // XMLish output...
  std::ostream& operator<<(std::ostream& os, const AtomicGroup& grp) {
    AtomicGroup::const_iterator i;
    if (grp.isPeriodic())
      os << "<GROUP PERIODIC='" << grp.box.box() << "'>\n";
    else
      os << "<GROUP>\n";
    for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
      os << "   " << **i << std::endl;
    os << "</GROUP>";

    return(os);
  }


  void AtomicGroup::setGroupConnectivity() {
    for (AtomicGroup::iterator i = atoms.begin(); i != atoms.end(); ++i)
      (*i)->setProperty(Atom::bondsbit);
  }



  void AtomicGroup::setCoords(double* seq, int m, int n) {
    if (n != 3 || static_cast<uint>(m) != size())
      throw(LOOSError("Invalid dimensions in AtomicGroup::setCoords()"));

    for (int j=0; j<m; ++j)
      for (int i=0; i<n; ++i)
	atoms[j]->coords()[i] = seq[j*n+i];
  }


  void AtomicGroup::getCoords(double** outseq, int* m, int* n) {
    double* dp = static_cast<double*>(malloc(size() * 3 * sizeof(double)));
    for (uint j=0; j<size(); ++j)
      for (uint i=0; i<3; ++i)
	dp[j*3+i] = atoms[j]->coords()[i];

    *m = size();
    *n = 3;

    *outseq = dp;
  }

  AtomicGroup AtomicGroup::centrifyByMolecule() const {
    std::vector<AtomicGroup> mols = splitByMolecule();
    AtomicGroup centers;
    for (std::vector<AtomicGroup>::const_iterator i = mols.begin(); i != mols.end(); ++i) {
      pAtom orig = (*i)[0];
      pAtom atom(new Atom(*orig));
      atom->name("CEN");
      atom->coords((*i).centerOfMass());
      centers.append(atom);
    }
    centers.box = box;
    return(centers);
  }

  AtomicGroup AtomicGroup::centrifyByResidue() const {
    std::vector<AtomicGroup> residues = splitByResidue();
    AtomicGroup centers;
    for (std::vector<AtomicGroup>::const_iterator i = residues.begin(); i != residues.end(); ++i) {
      pAtom orig = (*i)[0];
      pAtom atom(new Atom(*orig));
      atom->name("CEN");
      atom->coords((*i).centerOfMass());
      centers.append(atom);
    }
    centers.box = box;
    return(centers);
  }



}
