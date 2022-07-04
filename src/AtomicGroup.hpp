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





#if !defined(LOOS_ATOMICGROUP_HPP)
#define LOOS_ATOMICGROUP_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <limits>

#include <boost/unordered_set.hpp>


#include <loos_defs.hpp>


#include <Atom.hpp>
#include <XForm.hpp>
#include <PeriodicBox.hpp>
#include <utils.hpp>
#include <Matrix.hpp>
#include <FormFactor.hpp>
#include <FormFactorSet.hpp>

#include <exceptions.hpp>


namespace loos {


  //! Virtual base-class for selecting atoms from a group

  struct AtomSelector {
    //! Predicate function for selecting atoms.  If true, then the passed
    //! Atom is selected for an operation (or addition to a new group).
    //! If false, then the passed Atom is skipped.
    virtual bool operator()(const pAtom& atom) const =0;
    virtual ~AtomSelector() { }
  };

  //! hash pairs in an unordered way, that is 
  //! hash(pair(a,b)) == hash(pair(b,a))
  template<typename T>
  struct unordered_pair_hash {
    std::size_t operator () (std::pair<T, T> const &pair) const {
      // order the elements in the pair
      auto min_max = std::minmax(pair.first, pair.second);
      // return either
      return static_cast<size_t>(min_max.second) << sizeof(T)*CHAR_BIT | min_max.first;
    }
  };

//! test equality of unordered pair. typename T must have == operator.
  template<typename T>
  struct unordered_pair_eq {
    bool operator() (std::pair<T, T> const &p1, 
                     std::pair<T, T> const&p2) const {
      return p1 == p2 || (p1.first == p2.second && p1.second == p2.first);
    }
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
    typedef pAtom                              value_type;

    // Threshold for catching effectively zero singular values in
    // the superposition code...
    static const double superposition_zero_singular_value;

  public:
    AtomicGroup() : _sorted(false) { }

    //! Creates a new AtomicGroup with \a n un-initialized atoms.
    /** The atoms will all have ascending atomid's beginning with 1, but
     *  otherwise no other properties will be set.
     */
    AtomicGroup(const int n) : _sorted(true) {
      assert(n >= 1 && "Invalid size in AtomicGroup(n)");
      for (int i=1; i<=n; i++) {
        pAtom pa(new Atom);
        pa->id(i);
        atoms.push_back(pa);
      }
    }

    //! Copy constructor (atoms and box shared)
    AtomicGroup(const AtomicGroup& g) :
      _sorted(g._sorted),
      atoms(g.atoms),
      box(g.box)
      { }


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


    uint length(void) const { return(atoms.size()); }
    uint size(void) const { return(atoms.size()); }
    bool empty(void) const { return(atoms.empty()); }

    //! Get the ith atom from this group.
    pAtom getAtom(const int i) const;

#if !defined(SWIG)
    //! Same as getAtom(i)
    pAtom& operator[](const int i);
    const pAtom& operator[](const int i) const;
#endif

    //! Append the atom onto the group
    AtomicGroup& append(pAtom pa) { atoms.push_back(pa); _sorted = false; return(*this); }
    //! Append a vector of atoms
    AtomicGroup& append(std::vector<pAtom> pas);
    //! Append an entire AtomicGroup onto this one (concatenation)
    AtomicGroup& append(const AtomicGroup& grp);

    //! Delete a single atom
    AtomicGroup& remove(pAtom pa) { deleteAtom(pa); return(*this); }
    //! Deletes a set of atoms
    AtomicGroup& remove(std::vector<pAtom> pas);
    //! Deletes all atoms in the passed grp that are also in the current group.
    AtomicGroup& remove(const AtomicGroup& grp);

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

    //! Inequality test for two groups
    bool operator!=(AtomicGroup& rhs) {
      return(!(operator==(rhs)));
    }


#if !defined(SWIG)

    //! Equality test for const groups
    /**Similar to the non-const version, but it will sort <I>copies</I>
     *of the atom lists if they are not already sorted...
     */
    bool operator==(const AtomicGroup& rhs) const;

    //! Inequality test for two groups
    bool operator!=(const AtomicGroup& rhs) const {
      return(!(operator==(rhs)));
    }

#endif // !defined(SWIG)



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

    template<class EqualsOp> bool contains(const pAtom& p, const EqualsOp& op) const {
      const_iterator ci = std::find_if(begin(), end(), bind2nd(op, p));
      return(ci != end());
    }

    //! Determines if a pAtom is contained in this group using the AtomEquals policy (ie the default comparison policy)
    bool contains(const pAtom& p) const { return(contains(p, AtomEquals())); }


    //! Determines if the passed group is a subset of the current group using the EqualsOp atom-equality policy
    template<class EqualsOp> bool contains(const AtomicGroup& g, const EqualsOp& op) const {
      for (const_iterator cj = g.begin(); cj != g.end(); ++cj)
        if (std::find_if(begin(), end(), bind2nd(op, *cj)) == end())
          return(false);
      return(true);
    }

    //! Determines if a group is a subset of the current group using the default AtomEquals policy
    bool contains(const AtomicGroup& g) const { return(contains(g, AtomEquals())); }


      //! Determines if a group contains any atom
      template<class EqualsOp> bool containsAny(const AtomicGroup& g, const EqualsOp& op) const
          {
              for (const_iterator cj = g.begin(); cj != g.end(); ++cj)
                  if (std::find_if(begin(), end(), bind2nd(op, *cj)) != end())
                      return(true);
              return(false);
          }


      //! Determines if a group contains any atom using the default AtomEquals policy
      bool containsAny(const AtomicGroup& g) const { return(containsAny(g, AtomEquals())); }

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
    /**
     * The groups that are returned will be in the same order that the segids appear
     * in the source AtomicGroup
     */
    std::vector<AtomicGroup> splitByUniqueSegid(void) const;

    //! Returns a vector of AtomicGroups split based on bond connectivity
    std::vector<AtomicGroup> splitByMolecule(void) const {
      AtomicGroup sortable = *this;
      return(sortable.sortingSplitByMolecule());
    }

    //! Takes selection string as argument to be applied to each group after splitting.
    //! Returns a vector of AtomicGroups split based on bond connectivity;
    std::vector<AtomicGroup> splitByMolecule(const std::string& selection) const {
      AtomicGroup sortable = *this;
      return(sortable.sortingSplitByMolecule(selection));
    }

    //! Returns a vector of AtomicGroups, each comprising a single residue
    std::vector<AtomicGroup> splitByResidue(void) const;

    //! Returns a vector of AtomicGroups, each containing atoms with the same name
    std::map<std::string, AtomicGroup> splitByName(void) const;


    //! Replace a group with the center of masses of contained molecules
    /**
     * The AtomicGroup is split into molecules.  A new group is constructed
     * where each atom is the center of mass of one molecule.  The atom
     * metadata is taken from the first atom of the associated molecule,
     * but with the atom name "CEN".
     */
    AtomicGroup centrifyByMolecule() const;

    //! Replace a group with the cente of masses of contained residues (see centrifyByMolecule())
    AtomicGroup centrifyByResidue() const;


    //! Find a contained atom by its atomid
    /**
     * The default behavior is to assume that the atoms in the
     * AtomicGroup are not in order of increasing atomid and to
     * therefore use a linear search.  If the atoms are sorted
     * (AtomicGroup::sort()), then the more efficient binary search
     * will be used.
     */
    pAtom findById(const int id) const;

    //! Create a new group from a vector of atomids
    AtomicGroup groupFromID(const std::vector<int> &id_list) const;

    //! Create a new group from a pair of atomids
    AtomicGroup groupFromID(const std::pair<int, int> &id_pair) const;

    //! Given an Atom, return a group of all the atoms contained by its
    //! containing residue
    AtomicGroup getResidue(pAtom res);

#if !defined(SWIG)
    //! Output the group in pseudo-XML format...
    friend std::ostream& operator<<(std::ostream& os, const AtomicGroup& grp);
#endif

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

    //! Do all the atoms in the group have coordinates?
    bool hasCoords(void) const;

    //! Remove any bonding information present in contained atoms
    void clearBonds(void);

    //! Attempt to prune connectivity (only retain bonds to atoms within this AtomicGroup)
    void pruneBonds();

    //! Reset the atom indices (used for interfacing with trajectories)
    void resetAtomIndices();

    //! Deduce atomic number from mass (if present), returning number of atoms assigned
    uint deduceAtomicNumberFromMass(const double tol = 0.1);

    //! Is the array of atoms already sorted???
    /**
     * While we make some effort to ensure that alterations to the AtomicGroup
     * will invalidate the sorted status, it's a good idea to
     * explicitly sort if you want to make sure that the group is in
     * fact sorted.
     */
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

    //! compute OCF for all atom-pairs in AG of distance offset from one another
    const greal ocf(uint offset);

    //! Provide access to the underlying shared periodic box...
    loos::SharedPeriodicBox sharedPeriodicBox() const { return(box); }

    //! Remove periodicity
    void removePeriodicBox() { box = SharedPeriodicBox(); }

    //! Translate the entire group so that the centroid is in the
    //! primary cell
    void reimage();

    //! Reimage atoms individually into the primary cell
    void reimageByAtom();

    //! Takes a group that's split across a periodic boundary and reimages it so it's all together.
    void mergeImage(pAtom &p);
    //! Takes a group that's split across a periodic boundary and reimages it so it's all together, using the first atom in the AtomicGroup as the reference
    void mergeImage();

    //! Find atoms in the current group that are within \a dist angstroms of any atom in \a grp
    AtomicGroup within(const double dist, AtomicGroup& grp) const {
      Distance2WithoutPeriodicity op;
      return(within_private(dist, grp, op));
    }

    //! Find atoms in \a grp that are within \a dist angstroms of atoms in the current group, considering periodicity
    AtomicGroup within(const double dist, AtomicGroup& grp, const GCoord& box) const {
      Distance2WithPeriodicity op(box);
      return(within_private(dist, grp, op));
    }


    //! Returns true if any atom of current group is within \a dist angstroms of \a grp
    /**
     * \a min is the minimum number of pair-wise contacts required to be considered
     * in contact
     */
    bool contactWith(const double dist, const AtomicGroup& grp, const uint min=1) const {
      Distance2WithoutPeriodicity op;
      return(contactwith_private(dist, grp, min, op));
    }

    //! Returns true if any atom of current group is within \a dist angstroms of \a grp
    /**
     * \a min is the minimum number of pair-wise contacts required to be considered
     * in contact
     */
    bool contactWith(const double dist, const AtomicGroup& grp, const GCoord& box, const uint min=1) const {
      Distance2WithPeriodicity op(box);
      return(contactwith_private(dist, grp, min, op));
    }
    
    //! return a list of atom ID pairs that correspond to all unique bonds.          
    std::vector<std::pair<int, int>> getBondsIDs() const;

    //! return a list of atom index pairs corresponding to all unique bonds.
    std::vector<AtomicGroup> getBondsAGs() const;

    //! Distance-based search for bonds
    /** Searches for bonds within an AtomicGroup based on distance.
     *  does NOT clear the existing bond list prior to building new
	 *  bonds.  The default distance cutoff is 1.65.  If a box (GCoord)
	 *  is passed, then periodicity is taken into consideration.
     */
	// Larger distances cause problems with hydrogens...
	void findBonds(const double dist, const GCoord& box) { findBondsImpl(dist, Distance2WithPeriodicity(box)); }
	void findBonds(const double dist) { findBondsImpl(dist, Distance2WithoutPeriodicity()); }
	void findBonds(const GCoord& box) { findBondsImpl(1.65, Distance2WithPeriodicity(box)); }
	void findBonds() { findBondsImpl(1.65, Distance2WithoutPeriodicity()); }



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
    iterator end(void) { return(atoms.end()); }

#if !defined(SWIG)
    const_iterator begin(void) const { return(atoms.begin()); }
    const_iterator end(void) const { return(atoms.end()); }
#endif

    // Statistical routines...
    //! Bounding box for the group...
    /**
     * Returns a vector containing 2 GCoords, one containing
     * (minx, miny, minz) and the other (maxx, maxy, maxz)
     */
    std::vector<GCoord> boundingBox(void) const;

    //! Translates the group so that the centroid is at the origin.
    /**
     * Returns the old centroid of the group
     */
    GCoord centerAtOrigin(void);

    //! Centroid of atoms (ignores mass, operates in group coordinates)
    GCoord centroid(void) const;

    //! Maximum radius from centroid of all atoms (not gyration)
    /**
     *  If optional argument is true, uses coordinates of atom 0 instead of centroid.
     *  Argument is false by default.
     */
    greal radius(const bool use_atom_as_reference=false) const;

    //! Center of mass of the group (in group coordinates)
    GCoord centerOfMass(void) const;

    //! Analogous to center of mass
    GCoord centerOfElectrons(void) const;

    //! Dipole moment, relative to group's centroid
    GCoord dipoleMoment(void) const;
    greal totalCharge(void) const;
    greal totalMass(void) const;
    greal radiusOfGyration(void) const;

    //! Spherical variance of group with respect to target atom
    greal sphericalVariance(const pAtom) const;
    greal sphericalVariance(const GCoord) const;

    //! Estimate stacking, as between two nucleobases
    /** Algorithm: n1=normal to self; n2=normal to other, dx = difference between
               centroids
     *         stacking = (n1*n2)^2 *[(n1 + n2)/2 * dx]/|dx| * 1/1 + (dx/threshold)^6
     */
     greal stacking(const AtomicGroup&, const GCoord& box, const double threshold) const;

    //! Compute the RMSD between two groups
    /** Assumes a 1:1 correspondence between ith atoms.
     *  Does NOT transform the coordinates in any way.
     */
    greal rmsd(const AtomicGroup&);

    //! Compute kinetic energy of group
    /**
      * Assumes mass and velocity have been set.
      * Output units are kcal/mol
      */
    greal kineticEnergy();

    // Geometric transformations...

    //! Returns a vector of coordinates transformed by the passed XForm
    /**
     * Does not alter the group's coordinates...
     */
    std::vector<GCoord> getTransformedCoords(const XForm&) const;

    //! Compute difference vectors between two AtomicGroups
    /**
        Does not align the coordinates first
     */
    std::vector<GCoord> differenceVectors(const AtomicGroup &other);

    //! Translate an atomic group by vector v
    void translate(const GCoord & v);

    //! Rotate group's coordinates (right-handed, about centroid)
    void rotate(const GCoord& axis, const greal angle_in_degrees);

    //! Apply the given transform to the group's coordinates...
    void applyTransform(const XForm&);

    //! Copy coordinates from a vector of GCoords using the atom index as an index into the vector.
    void copyCoordinatesWithIndex(const std::vector<GCoord>& coords);

    //! Copy velocities from a vector of GCoords using the atom index as an index into the vector.
    /**
     * This can be used to update a group's velocities if they come from a separate trajectory...
     * \code
     * pTraj trajcrds = createTrajectory('foo.dcd', model);
     * pTraj trajvels = createTrajectory('foo-velocities.dcd', model);
     *
     * while (trajcrds->readFrame()) {
     *    trajcrds->updateGroupCoords(model);
     *
     *    trajvels->readFrame();
     *    model.copyVelocitiesWithIndex(trajvels->coords());
     * }
     * \endcode
     */
    void copyVelocitiesWithIndex(const std::vector<GCoord>& velocities);


    //! Copy coordinates from g into current group
    /**
     * The offset is relative to the start of the current group
     * and length is the number of coordinates to copy.  If
     * length + offset is greater than the number of atoms
     * in the current group, then the excess coordinates will not be
     * copied.
     *
     * If length is 0, then all coordinates in g
     * will be copied.
     *
     * It is assumed that the atoms in g are in the appropriate
     * order relative to the current group for the copy to make
     * sense.
     */
    void copyCoordinatesFrom(const AtomicGroup& g, const uint offset = 0, const uint length = 0);

    //! Map the order of atoms in AtomicGroup g into the current group
    /**
     * Note that the order is only checked within a residue.  The residues
     * must appear in the same order between the two groups.  This
     * addresses edge issues such as when psfgen reorders the atoms within
     * a residue.  The map is an index into the AtomicGroup g that
     * puts g into the same order as the current group.
     */
    std::vector<uint> atomOrderMapFrom(const AtomicGroup& g);

    //! Given a mapping of atom order, copy the coordinates into the current group
    /**
     * See AtomicGroup::atomOrderMapFrom(const AtomicGroup& g) for
     * more information
     *
     * If you know that the atoms are in the same order in both
     * groups, then AtomicGroup::copyCoordinatesFrom() will be faster...
     */
    void copyMappedCoordinatesFrom(const AtomicGroup& g, const std::vector<uint>& order);

    //! Copy the coordinates from the group mapping the atom order
    /**
     * See AtomicGroup::atomOrderFrom(const AtomicGroup& g) for more information
     *
     * If you know that the atoms are in the same order in both
     * groups, then AtomicGroup::copyCoordinatesFrom() will be faster...
     */
    void copyMappedCoordinatesFrom(const AtomicGroup& g);

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
    /**
     * Calculates the principal moments and principal axes (from the
     * moment of inertia).  This is distinct from the principalAxes()
     * function which calculates the distribution of points about the
     * centroid.
     */
    std::vector<GCoord> momentsOfInertia(void) const;

    //! Calculates the transformation matrix for superposition of groups.
    /**
     * Uses the Kabsch alignment method (via SVD) to calculate the
     * transformation matrix that superimposes the current group onto
     * the passed group.  Returns the matrix.
     *
     * If too few atoms are given for aligning, the correlation matrix
     * may become singular and return fewer than three eigenpairs.  If
     * this is detected, superposition() will throw a NumericalError.
     * The threshold for a zero-eigenvalue (really, a zero singular value)
     * is set in AtomicGroup::superposition_zero_singular_value
     */
    GMatrix superposition(const AtomicGroup&);

    //! Superimposes the current group onto the passed group.
    /**
     * Calls superposition to calculate the transformation matrix to
     * superimpose the current group onto the passed one, then applies the
     * transformation to the current group's coordinates.
     */
    GMatrix alignOnto(const AtomicGroup&);

    //! Orient the principal axis of this group along the supplied vector
    /**
     * The supplied vector does not need to be normalized.
     */
    void orientAlong(const GCoord &);

    // Set coordinates to an array
    /**
     * This function is meant for Numpy/swig use in setting the model's
     * coordinates.  The passed array is row-major.
     */
    void setCoords(double* seq, int m, int n);

    // Return a newly allocated array containing the current group's coordinates
    /**
     * This function is meant for Numpy/swig use.  It will store the current
     * model's coordinates into a newly allocated array (using malloc).
     * The caller is expected to manage the memory.
     */
    void getCoords(double** outseq, int* m, int* n);

    std::vector<double> coordsAsVector() const;

    // Compute the packing score between 2 AtomicGroups
    /**
     * The packing score is the sum of 1/r^6 over all pairs of atoms,
     * respecting periodicity.
     * Quantity first defined in  Grossfield, A., et al,
     * Proc. Nat. Acad. Sci. USA, 2006, 103, 4888-4893
     */
    double packingScore(const AtomicGroup& other, const GCoord &box, bool norm) const;


    //* Logistic contact function between this group and another
    /**
        Compute the degree of contact between the centroid of this AG
        and the centroid of another AG, using a smooth logistic
        function
        S = 1/(1 + dist/radius)**sigma
     */
    double logisticContact(const AtomicGroup& group, double radius,
                           int sigma, const GCoord& box) const;

    //* Similar to logisticContact() but the distance between reference
    //  group centroid and another group centroid is 2D Euclidean distance
    //  instead of 3D, ignoring the z-component. Useful when calculating
    //  number of contacts in a plane
    /**
        Compute the number of contacts between the centroid of this AG
        and the centroid of another AG, using a smooth logistic
        function
        S = 1/(1 + dist/radius)**sigma
     */
    double logisticContact2D(const AtomicGroup& group, double radius,
                           int sigma, const GCoord& box) const;

    //* Hard contact function between this group and another
    /**
        Compute contact value of another AG with respect to
        the centroid of a given AG, using a hard step
        function
        S = 1; iff dist <= radius; else 0
     */
    double hardContact(const AtomicGroup& group, double radius,
                           const GCoord& box) const;

    //* Similar to hardContact() but the distance between reference
    //  group centroid and another group centroid is 2D Euclidean distance
    //  instead of 3D. Useful when calculating number of contacts or local
    //  neighbor density in a plane
    /**
        Compute contact value of another AG with respect to
        the centroid of a given AG, using a hard step
        function
        S = 1; iff dist <= radius; else 0
     */
    double hardContact2D(const AtomicGroup& group, double radius,
                           const GCoord& box) const;

    //* Compute x-ray scattering intensity from this group
    /**
        Computes X-ray scattering as a function of q, using
        I(q) = \sum_(atom pair) F_i(q) F_j(q) sin (q d_ij)/ (q d_ij)

        This approximates scattering off of individual atoms. If you use this
        with explicit solvent, you will get truncation artifacts from the periodic
        box (although the code computes all distances using periodicity).

        Form factors are from Szaloki, X-ray Spectrometry (1996), V25, 21-28

     */
    std::vector<double> scattering(const double qmin, const double max,
                                   const uint numValues,
                                   loos::FormFactorSet &formFactors);

  private:

	// These are functors for calculating distance between two coords
    // without and with periodicity.  These can be passed to functions
    // that need to support both ways of calculating distances, such
    // was within_private() below...
    struct Distance2WithoutPeriodicity {
      double operator()(const GCoord& a, const GCoord& b) const {
        return(a.distance2(b));
      }
    };

    struct Distance2WithPeriodicity {
      Distance2WithPeriodicity(const GCoord& box) : _box(box) { }

      double operator()(const GCoord& a, const GCoord& b) const {
        return(a.distance2(b, _box));
      }

      GCoord _box;
    };


    // This function is to to remove code duplication in
    // logisticContacts() and logisticContacts2D().
    // Handle even and odd powers separately -- even can
    // avoid the sqrt
    // Sigh, this doesnt' seem to make it much faster...
    double logisticFunc(const GCoord& cent, const GCoord& other, double radius, int sigma, const GCoord& box) const{
        double prod;
        if (sigma % 2 == 0) {
            double distance2 = cent.distance2(other, box);
            double ratio = distance2/(radius*radius);
            prod = ratio;
            for (int j=0; j<(sigma/2)-1; ++j) {
                prod *= ratio;
            }
        }
        else {
            double distance = cent.distance(other, box);
            double ratio = distance/radius;
            prod = ratio;
            for (int j=0; j < sigma-1; ++j) {
                prod *= ratio;
            }
        }
        double sum = 1./(1. + prod);
        return(sum);
    }


    // Find all atoms in the current group that are within dist
    // angstroms of any atom in the passed group.  The distance
    // calculation is determined by the passed functor so that the
    // same code can be used for both periodic and non-periodic
    // coordinates.

    template <typename DistanceCalc>
    AtomicGroup within_private(const double dist, AtomicGroup& grp, const DistanceCalc& distance_functor) const {

      AtomicGroup res;
      res.box = box;

      double dist2 = dist * dist;
      std::vector<uint> indices;

      for (uint j=0; j<size(); j++) {
        GCoord c = atoms[j]->coords();
        for (uint i=0; i<grp.size(); i++) {
          if (distance_functor(c, grp.atoms[i]->coords()) <= dist2) {
            indices.push_back(j);
            break;
          }
        }
      }

      if (indices.size() == 0)
        return(res);

      for (std::vector<uint>::const_iterator ci = indices.begin(); ci != indices.end(); ++ci)
        res.addAtom(atoms[*ci]);

      return(res);
    }


    template<typename DistanceCalc>
    bool contactwith_private(const double dist, const AtomicGroup& grp, const uint min_contacts, const DistanceCalc& distance_function) const {
      double dist2 = dist * dist;
      uint ncontacts = 0;

      for (uint j = 0; j<size(); ++j) {
        GCoord c = atoms[j]->coords();
	    for (uint i = 0; i<grp.size(); ++i)
            if (distance_function(c, grp.atoms[i]->coords()) <= dist2)
	          if (++ncontacts >= min_contacts)
	             return(true);
      }
      return(false);
    }


	  //! Internal implementation of find bonds.
	  /**
	   * Takes a functor for calculating distances.  This can be PBC aware or not
	   */
	  template<typename DistanceCalc>
	  void findBondsImpl(const double dist, const DistanceCalc& distance_function) {
		  iterator ij;
		  double dist2 = dist * dist;
		  double current_dist2;

		  for (ij = begin(); ij != end() - 1; ++ij) {
			  iterator ii;
			  GCoord u = (*ij)->coords();

			  for (ii = ij + 1; ii != end(); ++ii) {
				  current_dist2 = distance_function(u, (*ii)->coords());
				  if (current_dist2 < dist2) {
					  (*ij)->addBond(*ii);
					  (*ii)->addBond(*ij);
				  }
			  }
		  }
	  }




    std::vector<AtomicGroup> sortingSplitByMolecule();
    std::vector<AtomicGroup> sortingSplitByMolecule(const std::string& selection);

    // *** Internal routines ***  See the .cpp file for details...
    void sorted(bool b) { _sorted = b; }

    pAtom findById_linearSearch(const int id) const;
    pAtom findById_binarySearch(const int id) const;

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

    typedef boost::unordered_set<int> HashInt;

    void walkBonds(AtomicGroup& mygroup, HashInt& seen, AtomicGroup& working, pAtom& moi);


    double *coordsAsArray(void) const;
    double *transformedCoordsAsArray(const XForm&) const;

    bool _sorted;


  protected:

    void setGroupConnectivity();


    std::vector<pAtom> atoms;
    loos::SharedPeriodicBox box;

  };

  AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
  AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);


}

#endif
