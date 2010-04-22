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
#include <algorithm>

#include <boost/random.hpp>

#include <AtomicGroup.hpp>


namespace loos {

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
      throw(std::runtime_error("Attempting to delete a non-existent atom"));

    atoms.erase(iter);
    _sorted = false;
  }


  // Append each atom from the passed vector onto this group...
  void AtomicGroup::append(std::vector<pAtom> pas) {
    std::vector<pAtom>::iterator i;

    for (i=pas.begin(); i != pas.end(); i++)
      atoms.push_back(*i);

    _sorted = false;
  }


  // Append all atoms from the passed group onto this one
  void AtomicGroup::append(const AtomicGroup& grp) {
    std::vector<pAtom>::const_iterator i;

    for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
      addAtom(*i);

    _sorted = false;
  }


  // Remove all atoms in the passed vector
  void AtomicGroup::remove(std::vector<pAtom> pas) {
    std::vector<pAtom>::iterator i;

    for (i=pas.begin(); i != pas.end(); i++)
      deleteAtom(*i);

    _sorted = false;
  }


  // Removes all atoms contained in the passed group from this one...
  void AtomicGroup::remove(const AtomicGroup& grp) {
    std::vector<pAtom>::const_iterator i;

    for (i=grp.atoms.begin(); i != grp.atoms.end(); i++)
      deleteAtom(*i);

    _sorted = false;
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


  // Internal: sort the atom array by atomid
  void AtomicGroup::sort(void) {
    CmpById comp;

    if (!sorted())
      std::sort(atoms.begin(), atoms.end(), comp);

    sorted(true);
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
    UniqueStrings unique;

    for (i = atoms.begin(); i != atoms.end(); i++)
      unique.add((*i)->segid());

    int n = unique.size();
    int j;
    std::vector<AtomicGroup> results(n);
    for (i = atoms.begin(); i != atoms.end(); i++) {
      j = unique.find((*i)->segid());
      if (j < 0)
        throw(std::runtime_error("Could not find an atom we already found..."));
      results[j].append(*i);
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
   * SUCKS
   */

  std::vector<AtomicGroup> AtomicGroup::splitByMolecule(void) {
    HashInt seen;                      // Track what atoms we've already
    // processed... 
    std::vector<AtomicGroup> molecules;
    AtomicGroup current;               // The molecule we're currently building...

    // If no connectivity, just return the entire group...
    if (!hasBonds()) {
      sort();
      molecules.push_back(*this);
    } else {

      int n = size();
      for (int i=0; i<n; i++) {
        HashInt::iterator it = seen.find(atoms[i]->id());
        if (it != seen.end())
          continue;
      
        walkBonds(current, seen, atoms[i]);
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


  void AtomicGroup::walkBonds(AtomicGroup& current, HashInt& seen, pAtom moi) {
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
      pAtom toi = findById(*citer);
      if (toi == 0)
        throw(std::runtime_error("Missing bonds while trying to walk the connectivity tree."));
      walkBonds(current, seen, toi);
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
    
    return(residues);
  }


  // Find an atom based on atomid
  // Returns 0 (null shared_ptr) if not found...
  pAtom AtomicGroup::findById(const int id) {
    sort();
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

  //! Note: when calling this, you'll want to make sure you use the 
  //! outermost group (eg the psf or pdb you used to create things, rather than
  //! using a subselection, unless you're sure the subsection contains these
  //! atoms as well.  The main use of this routine is to create a group of atoms
  //! bound to another atom.
  AtomicGroup AtomicGroup::groupFromID(const std::vector<int> &id_list) {
    AtomicGroup result;

    result.box = box;

    for (unsigned int i=0; i<id_list.size(); i++) {
      pAtom pa = findById(id_list[i]);
      if (!pa) throw(std::out_of_range("Atom id doesn't exist"));
      result.addAtom(pa);
    }
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


  // renumber the atomids of the group...
  void AtomicGroup::renumber(const int start, const int stride) {
    iterator i;
    int id = start;

    for (i=atoms.begin(); i != atoms.end(); i++, id += stride)
      (*i)->id(id);
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
      if ((*i)->resid() < max)
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


  void AtomicGroup::reimage() {
    if (!(isPeriodic()))
      throw(std::runtime_error("trying to reimage a non-periodic group"));
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
      throw(std::runtime_error("trying to reimage a non-periodic group"));
    const_iterator a;
    GCoord box = periodicBox();
    for (a=atoms.begin(); a!=atoms.end(); a++) {
      (*a)->coords().reimage(box);
    }
  }
    

  /** Uses a not-very-bright algorithm that compares all atoms against
   * all atoms... 
   */

  AtomicGroup AtomicGroup::within(const double dist, AtomicGroup& grp) {
    int na = size();
    int nb = grp.size();
    AtomicGroup res;

    res.box = box;
    double dist2 = dist * dist;
    std::vector<int> ids;

    for (int j=0; j<nb; j++) {
      for (int i=0; i<na; i++) {
        if (atoms[i]->coords().distance2(grp.atoms[j]->coords()) <= dist2)
          ids.push_back(atoms[i]->id());
      }
    }

    // Abort the rest if nothing was found...
    if (ids.size() == 0)
      return(res);

    std::vector<int> unique_ids;
    std::sort(ids.begin(), ids.end());
    std::vector<int>::const_iterator ci;
    int last_id = ids[0];
    unique_ids.push_back(last_id);
    for (ci = ids.begin()+1; ci != ids.end(); ci++)
      if (*ci != last_id) {
        last_id = *ci;
        unique_ids.push_back(last_id);
      }

    for (ci = unique_ids.begin(); ci != unique_ids.end(); ci++) {
      pAtom pa = findById(*ci);
      if (pa == 0)
        throw(std::logic_error("Cannot find a found atom in AtomicGroup::atomsWithin()"));

      res.addAtom(pa);
    }

    return(res);
  }


  void AtomicGroup::findBonds(const double dist) {
    AtomicGroup::iterator ij;
    double dist2 = dist * dist;
    
    for (ij = begin(); ij != end() - 1; ++ij) {
      AtomicGroup::iterator ii;
      GCoord u = (*ij)->coords();

      for (ii = ij + 1; ii != end(); ++ii) {
        if (u.distance2((*ii)->coords()) < dist2) {
          (*ij)->addBond(*ii);
          (*ii)->addBond(*ij);
        }
      }
    }
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


}
