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





#if !defined(ATOMICGROUP_HPP)
#define ATOMICGROUP_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>

#include <tr1/unordered_set>


#include <loos_defs.hpp>


#include <Atom.hpp>
#include <XForm.hpp>
#include <UniqueStrings.hpp>
#include <PeriodicBox.hpp>
#include <utils.hpp>
#include <Matrix.hpp>



namespace loos {


  //! Virtual base-class for selecting atoms from a group

  struct AtomSelector {
    //! Predicate function for selecting atoms.  If true, then the passed
    //! Atom is selected for an operation (or addition to a new group).
    //! If false, then the passed Atom is skipped.
    virtual bool operator()(const pAtom& atom) const =0;
    virtual ~AtomSelector() { }
  };


  class AtomicGroup;
  typedef boost::shared_ptr<AtomicGroup> pAtomicGroup;


  //! Class for handling groups of Atoms (pAtoms, actually)
  /** This class contains a collection of shared pointers to Atoms
   * (i.e. pAtoms).  Copying an AtomicGroup is a light-copy.  You can,
   * however, perform a deep copy by using the AtomicGroup::copy()
   * method.  Note that atomid's are assumed to be unique for any given
   * AtomicGroup.
   *
   * Valid operators are '+' and '+=' and can combine either
   * AtomicGroup objects or pAtom objects.
   *
   * AtomicGroups also support periodic boundary conditions via the
   * periodicBox() method.  If a box has been set, then isPeriodic()
   * will return true.  The periodic box is shared between the parent
   * group and all derived groups.  AtomicGroup copies have non-shared
   * periodic boxes...
   */


  class AtomicGroup {
  public:
    typedef std::vector<pAtom>::iterator       iterator;
    typedef std::vector<pAtom>::const_iterator const_iterator;

  public:
    AtomicGroup() : _sorted(false) { }

    //! Creates a new AtomicGroup with \a n un-initialized atoms.
    /** The atoms will all have ascending atomid's beginning with 1, but
     *  otherwise no other properties will be set.
     */
    AtomicGroup(const int n) : _sorted(true) {
      assert(n >= 1 && "Invalid size in AtomicGroup(n)");
      for (int i=1; i<n; i++) {
        pAtom pa(new Atom);
        pa->id(i);
        atoms.push_back(pa);
      }
    }

    virtual ~AtomicGroup() { }

    //! Creates a deep copy of this group
    /** This creates a non-polymorphic deep copy of an AtomicGroup.  The
     * additional catch is that it may end up involving extra
     * data-movement as the copy is constructed and then copied back out
     * to wherever you're putting it.
     */

    AtomicGroup copy(void) const;

    //! Creates a lightweight clone of this group (for polymorphism)
    /** Despite the name, this is meant for polymorphic use.  It is
     *  <I>not</I> a deep copy.  If you don't understand what any of
     *  this means, then you almost certainly want to be using the
     *  copy() method instead.
     */
    virtual AtomicGroup* clone(void) const;


    int length(void) const { return(atoms.size()); }
    int size(void) const { return(atoms.size()); }
    bool empty(void) const { return(atoms.empty()); }

    //! Get the ith atom from this group.
    pAtom getAtom(const int i) const;

    //! Same as getAtom(i)
    pAtom& operator[](const int i);
    const pAtom& operator[](const int i) const;

    //! Append the atom onto the group
    void append(pAtom pa) { atoms.push_back(pa); _sorted = false; }
    //! Append a vector of atoms
    void append(std::vector<pAtom> pas);
    //! Append an entire AtomicGroup onto this one (concatenation)
    void append(const AtomicGroup& grp);

    //! Delete a single atom
    void remove(pAtom pa) { deleteAtom(pa); }
    //! Deletes a set of atoms
    void remove(std::vector<pAtom> pas);
    //! Deletes all atoms in the passed grp that are also in the current group.
    void remove(const AtomicGroup& grp);

    // Concatenation of groups and/or atoms
    AtomicGroup& operator+=(const AtomicGroup& rhs);
    AtomicGroup& operator+=(const pAtom& rhs);
    AtomicGroup operator+(const AtomicGroup& rhs);
    AtomicGroup operator+(const pAtom& rhs);

    //! Equality test for two groups
    /**The test for equality is based on whether or not the contained
     *atom pointers are the same.  This operator will also force both
     *sides of the equation to be sorted.
     */
    bool operator==(AtomicGroup& rhs);

    //! Equality test for const groups
    /**Similar to the non-const version, but it will sort <I>copies</I>
     *of the atom lists if they are not already sorted...
     */
    bool operator==(const AtomicGroup& rhs) const;

    //! Inequality test for two groups
    bool operator!=(AtomicGroup& rhs) {
      return(!(operator==(rhs)));
    }

    //! Inequality test for two groups
    bool operator!=(const AtomicGroup& rhs) const {
      return(!(operator==(rhs)));
    }

    //! subset() and excise() args are patterned after perl's substr...
    /** If offset is negative, then it's relative to the end of the
     * group.  If length is 0, then everything from offset to the
     * appropriate end is used...
     */
    AtomicGroup subset(const int offset, const int len = 0);

    //! excise returns the excised atoms as a group...
    AtomicGroup excise(const int offset, const int len = 0);


    //! Determines if a pAtom is contained in this group using the EqualsOp atom-equality policy
    /**
     * The problem with determining containment/intersection/etc is
     * how to define when two atoms are equal...  This is done by
     * specifying a comparison functor (the EqualsOp) as a policy.
     * There are two comparison policies currently in LOOS: AtomEquals
     * and AtomCoordsEquals.  The default behavior is to use
     * AtomEquals which only compares a subset of the available Atom
     * metadata.  You can specify the more restrictive policy (or an
     * user-defined policy) like:
     * \code
     * bool b = group.contains(an_atom, loos::AtomCoordsEquals());
     * \endcode
     *
     * Or as another example, comparing only residue numbers...
     * \code
     * struct ResidEquals : public std::binary_function<pAtom, pAtom, bool> {
     *   bool operator()(const pAtom& a, const pAtom& b) { return(a.resid() == b.resid()); }
     * };
     *
     * bool b = group.contains(an_atom, ResidEquals());
     * \endcode
     */

    template<class EqualsOp> bool contains(const pAtom& p, const EqualsOp& op) {
      const_iterator ci = std::find_if(begin(), end(), bind2nd(op, p));
      return(ci != end());
    }
    
    //! Determines if a pAtom is contained in this group using the AtomEquals policy (ie the default comparison policy)
    bool contains(const pAtom& p) { return(contains(p, AtomEquals())); }


    //! Determines if the passed group is a subset of the current group using the EqualsOp atom-equality policy
    template<class EqualsOp> bool contains(const AtomicGroup& g, const EqualsOp& op) {
      for (const_iterator cj = g.begin(); cj != g.end(); ++cj)
        if (std::find_if(begin(), end(), bind2nd(op, *cj)) == end())
          return(false);
      return(true);
    }
    
    //! Determines if a group is a subset of the current group using the default AtomEquals policy
    bool contains(const AtomicGroup& g) { return(contains(g, AtomEquals())); }
    
    //! Computes the intersection of two groups using the EqualsOp atom-equality policy
    /**
     * See AtomicGroup::contains(const pAtom&, const EqualsOp&) for more details
     */
    template<class EqualsOp> AtomicGroup intersect(const AtomicGroup& g, const EqualsOp& op) {
      AtomicGroup result;

      for (const_iterator cj = begin(); cj != end(); ++ cj)
        if (std::find_if(g.begin(), g.end(), bind2nd(op, *cj)) != g.end())
          result.addAtom(*cj);

      result.box = box;
      return(result);
    }

    //! Intersection of two groups
    AtomicGroup intersect(const AtomicGroup& g) { return(intersect(g, AtomEquals())); }

    //! Union of two groups using the specified atom-equality policy
    /**
     * Note that the periodic box of the current group is unchanged by this operation
     */
    template<class EqualsOp> AtomicGroup merge(const AtomicGroup& g, const EqualsOp& op) {
      AtomicGroup result = copy();

      for (const_iterator ci = g.begin(); ci != g.end(); ++ci)
        if (std::find_if(begin(), end(), bind2nd(op, *ci)) == end())
          result.addAtom(*ci);

      return(result);
    }


    //! Union of two groups using the default AtomEquals atom-equality policy
    AtomicGroup merge(const AtomicGroup& g) { return(merge(g, AtomEquals())); }


    //! Return a group consisting of atoms for which sel predicate returns true...
    AtomicGroup select(const AtomSelector& sel) const;

    //! Returns a vector of AtomicGroups split from the current group based on segid
    std::vector<AtomicGroup> splitByUniqueSegid(void) const;

    //! Returns a vector of AtomicGroups split based on bond connectivity
    std::vector<AtomicGroup> splitByMolecule(void);

    //! Returns a vector of AtomicGroups, each comprising a single residue
    std::vector<AtomicGroup> splitByResidue(void) const;

    //! Returns a vector of AtomicGroups, each containing atoms with the same name
    std::map<std::string, AtomicGroup> splitByName(void) const;

    //! Find a contained atom by its atomid
    pAtom findById(const int id);

    //! Create a new group from a vector of atomids
    AtomicGroup groupFromID(const std::vector<int> &id_list);

    //! Given an Atom, return a group of all the atoms contained by its
    //! containing residue 
    AtomicGroup getResidue(pAtom res);

    //! Output the group in pseudo-XML format...
    friend std::ostream& operator<<(std::ostream& os, const AtomicGroup& grp);
  
    // Some misc support routines...

    //! Renumber the atomid's of the contained atoms...
    void renumber(const int start = 1, const int stride = 1);
    int minId(void) const;
    int maxId(void) const;
    int minResid(void) const;
    int maxResid(void) const;
    int numberOfResidues(void) const;
    int numberOfSegids(void) const;

    //! True if all atoms in the group have the passed property(ies)
    bool allHaveProperty(const Atom::bits& property) const;
    
    //! True if any atom in the group have the passed property(ies)
    bool anyHaveProperty(const Atom::bits& property) const;

    // These are now deprecated in favor of the above functions...
    //! Does any atom in the group have bond information???
    bool hasBonds(void) const;

    //! Does all the atoms in the group have coordinates?
    bool hasCoords(void) const;

    //! Remove any bonding information present in contained atoms
    void clearBonds(void);

    //! Is the array of atoms already sorted???
    bool sorted(void) const { return(_sorted); }

    //! Sort based on atomid
    void sort(void);

    //! Test whether or not periodic boundary conditions are set
    bool isPeriodic(void) const { return(box.isPeriodic()); }
  
    //! Fetch the periodic boundary conditions.
    GCoord periodicBox(void) const { return(box.box()); }

    //! Set the periodic boundary conditions.  
    void periodicBox(const GCoord& c) { box.box(c); }

    //! Set the periodic boundary conditions
    void periodicBox(const greal x, const greal y, const greal z) { 
      box.box(GCoord(x,y,z));
    }

    //! Translate the entire group so that the centroid is in the 
    //! primary cell
    void reimage();
  
    //! Reimage atoms individually into the primary cell
    void reimageByAtom();
  
    //! Find atoms in \a grp that are within \a dist angstroms of atoms
    //! in the current group.
    AtomicGroup within(const double dist, AtomicGroup& grp);

    //! Apply a functor or a function to each atom in the group.
    /** apply() let's you apply a functor or a function pointer to each
     * atom in the group.  The functor is passed a pAtom.  The functor
     * object is also returned (in case it retained state).  For
     * example, the following code snippet shows how to calculate the
     * centroid of a group using apply and a functor...
     \code
     struct Functor {
     Functor() : C(GCoord(0,0,0)), n(0) { }
     void operator()(pAtom& p) { C += p->coords(); ++n; }
     GCoord center(void) const { return(C/n); }

     GCoord C;
     int n;
     };

     Functor f = group.apply(Functor());
     GCoord centroid = f.center();
     \endcode
    */


    //! Distance-based search for bonds
    /** Searches for bonds within an AtomicGroup based on distance.
     *  does NOT clear the existing bond list prior to building new
     *  bonds.  The default distance cutoff is 1.25
     */

    // Larger distances cause problems with hydrogens...
    void findBonds(const double dist = 1.65);


    template<class T> T apply(T func) {
      for (iterator i = atoms.begin(); i != atoms.end(); ++i)
        func(*i);
      return(func);
    }

    // *** Helper classes...

    //! Our own simple iterator for stepping over all managed atoms.
    /** Example:
     *  \code
     *  AtomicGroup::Iterator iter(an_atomic_group);
     *  pAtom p;
     *
     *   while (p = iter())
     *    do_something(p);
     *  \endcode
     *
     *  Note that the shared atom returned is a copy of the shared-atom
     *  pointer stored, rather than a ref to the shared atom pointer...
     *  You should exercise GREAT care in modifying the atom while
     *  iterating, or performing any operations that modify the group
     *  you're iterating over.  In fact, don't do it, unless you are
     *  sure you know what you're doing.
     */
    class Iterator {
    public:
      explicit Iterator(const AtomicGroup& grp) : iter(grp.atoms.begin()), final(grp.atoms.end()) { }
      pAtom operator()(void) {
        if (iter >= final)
          return(pAtom());
        return(*iter++);
      }
    private:
      std::vector<pAtom>::const_iterator iter, final;
    };

    // STL-iterator access
    // Should these reset sort status?
    iterator begin(void) { return(atoms.begin()); }
    const_iterator begin(void) const { return(atoms.begin()); }

    iterator end(void) { return(atoms.end()); }
    const_iterator end(void) const { return(atoms.end()); }


    // Statistical routines...
    //! Bounding box for the group...
    std::vector<GCoord> boundingBox(void) const;

    //! Translates the group so that the centroid is at the origin.
    /**
     * Returns the old centroid of the group
     */
    GCoord centerAtOrigin(void);

    //! Centroid of atoms (ignores mass, operates in group coordinates)
    GCoord centroid(void) const;

    //! Maximum radius from centroid of all atoms (not gyration)
    greal radius(void) const;

    //! Center of mass of the group (in group coordinates)
    GCoord centerOfMass(void) const;
    GCoord centerOfCharge(void) const;
    GCoord centerOfElectrons(void) const;
    GCoord dipoleMoment(void) const;
    greal totalCharge(void) const;
    greal totalMass(void) const;
    greal radiusOfGyration(void) const;

    //! Compute the RMSD between two groups
    /**Sorts both groups (if necessary), then assumes a 1:1
     *correspondence between ith atoms.
     */
    greal rmsd(AtomicGroup&);

    // Geometric transformations...
  
    //! Returns a vector of coordinates transformed by the passed XForm
    /**
     * Does not alter the group's coordinates...
     */
    std::vector<GCoord> getTransformedCoords(const XForm&) const;
  
    void translate(const GCoord & v);

    //! Apply the given transform to the group's coordinates...
    void applyTransform(const XForm&);


    //! Copy coordinates from one group into another...
    /** Requires that the groups be the same size and that the ith atom
     * in group g matches the ith atom in the current group.
     */
  
    void copyCoordinates(AtomicGroup& g) {
      iterator i, j;

      for (i = atoms.begin(), j = g.atoms.begin(); i != atoms.end(); i++, j++)
        (*i)->coords((*j)->coords());
    }

    //! Each atom is moved in a random direction by a vector of the passed size
    void perturbCoords(const greal);

    //! Compute the principal axes of a group
    /** Calculates the eigendecomposition of AA' where A is column-wise
     * concatenation of coordinates from all atoms in the group.  The mean
     * coordinate is automatically subtracted from A...  Returns a vector
     * of GCoord's in order of decreasing magnitude of the corresponding
     * eigenvalue.  The eigenvalues are returned as a GCoord after the
     * eigenvectors.
     *
     * Example
     * \code
     *     vector<GCoord> V = group_of_atoms.principalAxes();
     *     GCoord eigenvalues = V[3];
     *     GCoord first_eigenvector = V[0];   // Most significant
     *     GCoord second_eigenvector = V[1];
     *     GCoord third_eigenvector = V[2];   // Least significant
     * \endcode
     *
     * Notes
     *  - Any errors encountered in the BLAS/LAPACK routines cause
     *    a runtime exception to be thrown...
     * 
     *  - Coord type of contained atoms will always be upcast to double.
     *
     *  - Potential issue with f77int under linux when not on a 64-bit
     *    architecture. 
     */
    std::vector<GCoord> principalAxes(void) const;


    //! Computes the moments of inertia for a group
    std::vector<GCoord> momentsOfInertia(void) const;

    //! Calculates the transformation matrix for superposition of groups.
    /**
     * Uses the Kabsch alignment method (via SVD) to calculate the
     * transformation matrix that superimposes the current group onto
     * the passed group.  Returns the matrix.
     */
    GMatrix superposition(const AtomicGroup&);

    //! Superimposes the current group onto the passed group.
    /**
     * Calls superposition to calculate the transformation matrix to
     * superimpose the current group onto the passed one, then applies the
     * transformation to the current group's coordinates.
     */
    GMatrix alignOnto(const AtomicGroup&);

  private:

    // *** Internal routines ***  See the .cpp file for details...
    void sorted(bool b) { _sorted = b; }

    int rangeCheck(int) const;

    void addAtom(pAtom pa) { atoms.push_back(pa); _sorted = false; }
    void deleteAtom(pAtom pa);

    boost::tuple<iterator, iterator> calcSubsetIterators(const int offset, const int len = 0);

    void copyCoordinatesById(AtomicGroup& g);

    // Some helper classes for using the STL
    struct CmpById {
      bool operator()(const pAtom& a, const pAtom& b) {
        return(a->id() < b->id());
      }
    };

    struct BindId {
      BindId(const int i) : id(i) { }
      bool operator()(const pAtom& a) { return(a->id() == id); }
      int id;
    };
    
    typedef std::tr1::unordered_set<int> HashInt;

    void walkBonds(AtomicGroup& mygroup, HashInt& seen, pAtom moi);


    double *coordsAsArray(void) const;
    double *transformedCoordsAsArray(const XForm&) const;

    bool _sorted;

  protected:
    std::vector<pAtom> atoms;
    SharedPeriodicBox box;

  };

  AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
  AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);


}

#endif
