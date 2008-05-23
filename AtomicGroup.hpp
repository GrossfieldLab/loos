/*
  AtomicGroup.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic class for groups of atoms...

  Notes:

    o Applying a transform to the coordinates of the group's atoms
      does NOT reset the transformation back to identity.  That is up
      to you to do if you want to do that...
*/




#if !defined(ATOMICGROUP_HPP)
#define ATOMICGROUP_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <tr1/memory>

#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace tr1;


#include <loos.hpp>
#include <Atom.hpp>
#include <XForm.hpp>



//! Virtual base-class for selecting atoms from a group

struct AtomSelector {
  //! Functor function for selecting atoms.  If true, then the passed
  //! Atom is selected for an operation (or addition to a new group).
  //! If false, then the passed Atom is skipped.
  virtual bool operator()(const pAtom& atom) const =0;
  virtual ~AtomSelector() { }
};


//! Class for handling groups of Atoms (pAtoms, actually)
//! This class contains a collection of shared pointers to Atoms
//! (i.e. pAtoms).  Copying an AtomicGroup is a light-copy.  You can,
//! however, perform a deep copy by using the AtomicGroup::clone()
//! method.  Note that atomid's are assumed to be unique for any given
//! AtomicGroup.
//!
//! The AtomicGroup also contains a XForm object for specifying
//! geometric transforms for the contained Atoms.  You can apply this
//! transform to the Atoms themselves, but be aware that doing so does
//! not reset the transform back to identity.
//!
//! Valid operators are '+' and '+=' and can combine either
//! AtomicGroup objects or pAtom objects.


class AtomicGroup {
protected:
  typedef vector<pAtom>::iterator AtomIterator;
  typedef vector<pAtom>::const_iterator ConstAtomIterator;

public:
  AtomicGroup() : _sorted(false) { }
  virtual ~AtomicGroup() { }

  //! Creates a deep copy of this group...
  virtual AtomicGroup* clone(void) const;


  int length(void) const { return(atoms.size()); }
  int size(void) const { return(atoms.size()); }

  //! Get the ith atom from this group.
  pAtom getAtom(const int i) const;

  //! Same as getAtom(i)
  pAtom& operator[](const int i);
  const pAtom& operator[](const int i) const;

  XForm& xform(void) { return(_xform); };

  //! Append the atom onto the group
  void append(pAtom pa) { atoms.push_back(pa); }
  //! Append a vector of atoms
  void append(vector<pAtom> pas);
  //! Append an entire AtomicGroup onto this one (concatenation)
  void append(const AtomicGroup& grp);

  //! Delete a single atom
  void remove(pAtom pa) { deleteAtom(pa); }
  //! Deletes a set of atoms
  void remove(vector<pAtom> pas);
  //! Deletes all atoms in the passed grp that are also in the current group.
  void remove(const AtomicGroup& grp);

  // Concatenation of groups and/or atoms
  AtomicGroup& operator+=(const AtomicGroup& rhs);
  AtomicGroup& operator+=(const pAtom& rhs);
  AtomicGroup operator+(const AtomicGroup& rhs);
  AtomicGroup operator+(const pAtom& rhs);

  //! subset() and excise() args are patterned after perl's substr...
  //! If offset is negative, then it's relative to the end of the
  //! group.  If length is 0, then everything from offset to the
  //! appropriate end is used...
  AtomicGroup subset(const int offset, const int len = 0);

  //! excise returns the excised atoms as a group...
  AtomicGroup excise(const int offset, const int len = 0);

  //! Intersection of two groups
  AtomicGroup intersect(const AtomicGroup& grp);

  //! Return a group consisting of atoms for which sel functor returns true...
  AtomicGroup select(AtomSelector& sel);

  //! Find a contained atom by its atomid
  pAtom findById(const int id);

  //! Given an Atom, return a group of all the atoms contained by its
  //! containing residue 
  AtomicGroup getResidue(pAtom res);

  //! Output the group in pseudo-XML format...
  friend ostream& operator<<(ostream& os, const AtomicGroup& grp);
  
  // Some misc support routines...

  //! Renumber the atomid's of the contained atoms...
  void renumber(const int start = 0, const int stride = 1);
  int minId(void) const;
  int maxId(void) const;
  int numberOfResidues(void) const;
  int numberOfChains(void) const;

  //! Is the array of atoms already sorted???
  bool sorted(void) const { return(_sorted); }

  //! Sort based on atomid
  void sort(void);


  // *** Helper classes...

  //! Our own simple iterator for stepping over all managed atoms.
  /*! Example:
      \verbatim
      AtomicGroup::Iterator iter(an_atomic_group);
      pAtom p;

      while (p = iter())
        do_something(p);
      \endverbatim

      Note that the shared atom returned is a copy of the shared-atom
      pointer stored, rather than a ref to the shared atom pointer...
      You should exercise GREAT care in modifying the atom while
      iterating, or performing any operations that modify the group
      you're iterating over.  In fact, don't do it, unless you are
      sure you know what you're doing.
   */
  class Iterator {
  public:
    Iterator(const AtomicGroup& grp) : iter(grp.atoms.begin()), final(grp.atoms.end()) { }
    pAtom operator()(void) {
      if (iter >= final)
	return(pAtom());
      return(*iter++);
    }
  private:
    vector<pAtom>::const_iterator iter, final;
  };


  //! Return value for boundingBox()
  struct BoundingBox {
    greal min[3], max[3];
  };

  // Statistical routines...
  //! Bounding box for the group...
  BoundingBox boundingBox(void) const;

  //! Centroid of atoms (ignores mass)
  GCoord centroid(void) const;

  //! Maximum radius from centroid of all atoms (not gyration)
  greal radius(void) const;

  GCoord centerOfMass(void) const;
  greal totalCharge(void) const;
  greal totalMass(void) const;
  greal radiusOfGyration(void) const;

  // Geometric transformations...
  
  //! Returns a vector of coordinates transformed by the current
  //! groups XForm
  vector<GCoord> transformedCoords(void) const;
  
  
  //! Applies the current XForm transformation to the contained Atom
  //! GCoord() coordinates.  This does NOT reset the XForm to the
  //! identity...
  void applyTransformation(void);


  //! Copy coordinates from one group into another...
  //! If the groups match in size, then a straight copy ensues.
  //! Otherwise, an attempt will be made to pick the correct
  //! coordinates...this could be a pretty costly operation...  Also
  //! note that this may change (sort) the group 'g'...
  
  void copyCoordinates(AtomicGroup& g);

private:

  // *** Internal routines ***  See the .cpp file for details...
  void sorted(bool b) { _sorted = b; }

  int rangeCheck(int) const;

  void addAtom(pAtom pa) { atoms.push_back(pa); }
  void deleteAtom(pAtom pa);

  AtomIterator findIteratorById(const int id);

  boost::tuple<AtomIterator, AtomIterator> calcSubsetIterators(const int offset, const int len = 0);

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


  bool _sorted;

protected:
  vector<pAtom> atoms;
  XForm _xform;

};

AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);

#endif
