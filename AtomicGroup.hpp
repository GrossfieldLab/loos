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



#if defined(__linux__)

extern "C" {

#include <atlas/cblas.h>
#include <atlas/clapack.h>

  void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
  void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);

}


// May need to fix this???
typedef int f77int;

#elif defined(__APPLE__)

#include <vecLib/vecLib.h>

typedef __CLPK_integer f77int;

#else

#warning Principal axes support will not be built in to the AtomicGroup class

#endif



using namespace std;
using namespace tr1;


#include <loos.hpp>
#include <Atom.hpp>
#include <XForm.hpp>
#include <UniqueStrings.hpp>


//! Virtual base-class for selecting atoms from a group

struct AtomSelector {
  //! Predicate function for selecting atoms.  If true, then the passed
  //! Atom is selected for an operation (or addition to a new group).
  //! If false, then the passed Atom is skipped.
  virtual bool operator()(const pAtom& atom) const =0;
  virtual ~AtomSelector() { }
};


//! Class for handling groups of Atoms (pAtoms, actually)
/** This class contains a collection of shared pointers to Atoms
 * (i.e. pAtoms).  Copying an AtomicGroup is a light-copy.  You can,
 * however, perform a deep copy by using the AtomicGroup::clone()
 * method.  Note that atomid's are assumed to be unique for any given
 * AtomicGroup.
 *
 * The AtomicGroup also contains a XForm object for specifying
 * geometric transforms for the contained Atoms.  You can apply this
 * transform to the Atoms themselves, but be aware that doing so does
 * not reset the transform back to identity.  Also, some geometric
 * operations work in group coordinate space (i.e. untransformed)
 * whereas some will operate in world coordinates (i.e. transformed).
 *
 * If you need to fetch the coordinates for an atom in a group, it is
 * <I>strongly</I> recommended that you use the
 * getAtomsTransformedCoord() method.  This will return a GCoord that
 * has been transformed by the current transformation matrix.  If you
 * access the atom directly, it does not know about the transformation
 * and you will get incorrect coordinates.
 * 
 * Valid operators are '+' and '+=' and can combine either
 * AtomicGroup objects or pAtom objects.
 *
 * AtomicGroups also support periodic boundary conditions via the
 * periodicBox() method.  If a box has been set, then isPeriodic()
 * will return true.
*/


class AtomicGroup {
protected:
  typedef vector<pAtom>::iterator AtomIterator;
  typedef vector<pAtom>::const_iterator ConstAtomIterator;

public:
  AtomicGroup() : _sorted(false), _periodic(false) { }
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

  //! subset() and excise() args are patterned after perl's substr...
  /** If offset is negative, then it's relative to the end of the
   * group.  If length is 0, then everything from offset to the
   * appropriate end is used...
   */
  AtomicGroup subset(const int offset, const int len = 0);

  //! excise returns the excised atoms as a group...
  AtomicGroup excise(const int offset, const int len = 0);

  //! Intersection of two groups
  AtomicGroup intersect(const AtomicGroup& grp);

  //! Return a group consisting of atoms for which sel predicate returns true...
  AtomicGroup select(const AtomSelector& sel);

  //! Returns a vector of AtomicGroups split from the current group based on segid
  vector<AtomicGroup> splitByUniqueSegid(void) const;

  //! Find a contained atom by its atomid
  pAtom findById(const int id);

  //! Create a new group from a vector of atomids
  AtomicGroup groupFromID(const vector<int> &id_list);

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
  int minResid(void) const;
  int maxResid(void) const;
  int numberOfResidues(void) const;
  int numberOfChains(void) const;

  //! Is the array of atoms already sorted???
  bool sorted(void) const { return(_sorted); }

  //! Sort based on atomid
  void sort(void);

  //! Test whether or not periodic boundary conditions are set
  bool isPeriodic(void) const { return(_periodic); }
  
  //! Fetch the periodic boundary conditions
  GCoord periodicBox(void) const { return(box); }
  
  //! Set the periodic boundary conditions
  void periodicBox(const GCoord& c) { _periodic = true; box = c; }
  


  // *** Helper classes...

  //! Our own simple iterator for stepping over all managed atoms.
  /** Example:
   *  \verbatim
   *  AtomicGroup::Iterator iter(an_atomic_group);
   *  pAtom p;
   *
   *   while (p = iter())
   *    do_something(p);
   *  \endverbatim
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

  //! Translates (via XForm) the group so that the centroid is at the origin.
  void centerAtOrigin(void);

  //! Centroid of atoms (ignores mass, operates in group coordinates)
  GCoord centroid(void) const;

  //! Maximum radius from centroid of all atoms (not gyration)
  greal radius(void) const;

  //! Center of mass of the group (in group coordinates)
  GCoord centerOfMass(void) const;
  GCoord centerOfCharge(void) const;
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
  
  //! Returns a vector of coordinates transformed by the current
  //! groups XForm
  vector<GCoord> transformedCoords(void) const;

  //! Return the ith atom's coordinates transformed by the current XForm
  GCoord getAtomsTransformedCoord(int) const;
  
  
  //! Applies the current XForm transformation to the contained Atom
  //! GCoord() coordinates.  This does NOT reset the XForm to the
  //! identity...
  void applyTransform(void);


  //! Copy coordinates from one group into another...
  /** Requires that the groups be the same size and that the ith atom
   * in group g matches the ith atom in the current group.
   */
  
  void copyCoordinates(AtomicGroup& g) {
    AtomIterator i, j;

    for (i = atoms.begin(), j = g.atoms.begin(); i != atoms.end(); i++, j++)
      (*i)->coords((*j)->coords());
  }

  //! Each atom is moved in a random direction by a vector of the passed size
  void perturbCoords(const greal);


#if defined(__linux__) || defined(__APPLE__)
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
   *
   *  - Operates in world coordinates (i.e. transformed group coordinates)
   * 
   */
  vector<GCoord> principalAxes(void) const;

  //! Superimposes the current group atop another group.
  /**Computes the Kabsch superposition using an SVD of the current
   *group with the passed group.  The current XForm is altered to affect
   *the superposition.  The transformation matrix (not including
   *centering the current group) is returned as a GMatrix.
   */
  GMatrix alignOnto(AtomicGroup&);

#endif

private:

  // *** Internal routines ***  See the .cpp file for details...
  void sorted(bool b) { _sorted = b; }

  int rangeCheck(int) const;

  void addAtom(pAtom pa) { atoms.push_back(pa); }
  void deleteAtom(pAtom pa);

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


  void dumpMatrix(const string, double*, int, int) const;
  double *coordsAsArray(void) const;
  double *transformedCoordsAsArray(void) const;

  bool _sorted;
  bool _periodic;

protected:
  vector<pAtom> atoms;
  XForm _xform;
  GCoord box;

};

AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);

#endif
