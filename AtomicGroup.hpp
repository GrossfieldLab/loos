/*
  AtomicGroup.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic class for groups of atoms...
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



// Subclass this for programmatic selection of atoms...

struct AtomSelector {
  virtual bool operator()(const pAtom& atom) const = 0;
  virtual ~AtomSelector() { }
};




class AtomicGroup {
protected:
  typedef vector<pAtom>::iterator AtomIterator;
  typedef vector<pAtom>::const_iterator ConstAtomIterator;

public:
  AtomicGroup() : _sorted(false) { }
  virtual ~AtomicGroup() { }

  // clone() creates a deep copy of this group...
  virtual AtomicGroup* clone(void) const;

  int length(void) const { return(atoms.size()); }
  int size(void) const { return(atoms.size()); }

  pAtom getAtom(const int i) const;

  pAtom& operator[](const int i);
  const pAtom& operator[](const int i) const;

  void append(pAtom pa) { atoms.push_back(pa); }   // Add an atom...
  void append(vector<pAtom> pas);
  void append(const AtomicGroup& grp);

  void remove(pAtom pa) { deleteAtom(pa); }        // Delete an atom...
  void remove(vector<pAtom> pas);
  void remove(const AtomicGroup& grp);

  // Concatenation of groups and/or atoms
  AtomicGroup& operator+=(const AtomicGroup& rhs);
  AtomicGroup& operator+=(const pAtom& rhs);
  AtomicGroup operator+(const AtomicGroup& rhs);
  AtomicGroup operator+(const pAtom& rhs);

  // subset() and excise() args are patterned after perl's substr...
  // If offset is negative, then it's relative to the end of the
  // group.  If length is 0, then everything from offset to the
  // appropriate end is used...
  AtomicGroup subset(const int offset, const int len = 0);

  // excise returns the excised atoms as a group...
  AtomicGroup excise(const int offset, const int len = 0);

  AtomicGroup intersect(const AtomicGroup& grp);

  // Return a group consisting of atoms for which sel() returns true...
  AtomicGroup select(const AtomSelector& sel);

  // Find a contained atom by its atomid
  pAtom findById(const int id);

  // Given an Atom, return a group of all the atoms contained by its
  // containing residue 
  AtomicGroup getResidue(pAtom res);

  // Output...  This is in an XML'ish format...
  friend ostream& operator<<(ostream& os, const AtomicGroup& grp);
  
  // Some misc support routines...

  void renumber(const int start = 0, const int stride = 1);
  int minId(void) const;
  int maxId(void) const;
  int numberOfResidues(void) const;
  int numberOfChains(void) const;

  // Is the array of atoms already sorted???
  bool sorted(void) const { return(_sorted); }

  // Sort based on atomid
  void sort(void);


  // *** Helper classes...

  // Iterator lets you iterate over all contained atoms...
  // AtomicGroup::Iterator(grp)  -- constructs an iterator for the
  // group "grp" 
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


  // Return value for boundingBox()
  struct BoundingBox {
    greal min[3], max[3];
  };

  // Statistical routines...
  BoundingBox boundingBox(void) const;

  // Centroid of atoms (ignores mass)
  GCoord centroid(void) const;

  // Maximum radius from centroid of all atoms (not gyration)
  greal radius(void) const;

  GCoord centerOfMass(void) const;
  greal totalCharge(void) const;
  greal totalMass(void) const;
  greal radiusOfGyration(void) const;

private:

  // *** Internal routines ***  See the .cpp file for details...
  void sorted(bool b) { _sorted = b; }

  int rangeCheck(int) const;

  void addAtom(pAtom pa) { atoms.push_back(pa); }
  void deleteAtom(pAtom pa);

  AtomIterator findIteratorById(const int id);

  boost::tuple<AtomIterator, AtomIterator> calcSubsetIterators(const int offset, const int len = 0);

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

};

AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);

#endif
