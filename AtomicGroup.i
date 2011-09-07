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






%header %{
#include <AtomicGroup.hpp>
#include <sfactories.hpp>
#include <Trajectory.hpp>

#include <Selectors.hpp>
#include <Parser.hpp>
#include <utils.hpp>

#include <Kernel.hpp>
#include <KernelValue.hpp>
#include <KernelActions.hpp>
#include <KernelStack.hpp>
#include <Selectors.hpp>

#include <pdb_remarks.hpp>

#include <sstream>
%}

%wrapper %{
typedef double    greal;
typedef loos::Matrix44<double>   GMatrix;
typedef unsigned int   uint;
typedef unsigned long  ulong;
 %}

%template(AtomicGroupVector) std::vector<loos::AtomicGroup>;


namespace loos {


  struct AtomSelector {
    virtual bool operator()(const pAtom& atom) const =0;
    virtual ~AtomSelector() { }
  };



  class AtomicGroup {
  public:
    typedef std::vector<pAtom>::iterator       iterator;
    typedef std::vector<pAtom>::const_iterator const_iterator;

  public:
    AtomicGroup();
    AtomicGroup(const int n);
    virtual ~AtomicGroup();
    AtomicGroup copy(void) const;
    virtual AtomicGroup* clone(void) const;

    uint length(void) const;
    uint size(void) const;
    bool empty(void) const;
    pAtom getAtom(const int i) const;
    AtomicGroup& append(pAtom pa);
    AtomicGroup& append(std::vector<pAtom> pas);
    AtomicGroup& append(const AtomicGroup& grp);
    AtomicGroup& remove(pAtom pa);
    AtomicGroup& remove(std::vector<pAtom> pas);
    AtomicGroup& remove(const AtomicGroup& grp);
    AtomicGroup& operator+=(const AtomicGroup& rhs);
    AtomicGroup& operator+=(const pAtom& rhs);
    AtomicGroup operator+(const AtomicGroup& rhs);
    AtomicGroup operator+(const pAtom& rhs);
    bool operator==(AtomicGroup& rhs);
    bool operator!=(AtomicGroup& rhs);
    AtomicGroup subset(const int offset, const int len = 0);
    AtomicGroup excise(const int offset, const int len = 0);
    bool contains(const pAtom& p);
    bool contains(const AtomicGroup& g);
    AtomicGroup intersect(const AtomicGroup& g);
    AtomicGroup merge(const AtomicGroup& g);
    AtomicGroup select(const AtomSelector& sel) const;
    std::vector<AtomicGroup> splitByUniqueSegid(void) const;
    std::vector<AtomicGroup> splitByMolecule(void);
    std::vector<AtomicGroup> splitByResidue(void) const;
    std::map<std::string, AtomicGroup> splitByName(void) const;
    pAtom findById(const int id);
    AtomicGroup groupFromID(const std::vector<int> &id_list);
    AtomicGroup getResidue(pAtom res);
    void renumber(const int start = 1, const int stride = 1);
    int minId(void) const;
    int maxId(void) const;
    int minResid(void) const;
    int maxResid(void) const;
    int numberOfResidues(void) const;
    int numberOfSegids(void) const;


    bool allHaveProperty(const Atom::bits& property) const;
    bool anyHaveProperty(const Atom::bits& property) const;
    bool hasBonds(void) const;
    bool hasCoords(void) const;
    void clearBonds(void);
    void pruneBonds();

    uint deduceAtomicNumberFromMass(const double tol = 0.1);
    bool sorted(void) const;
    void sort(void);

    bool isPeriodic(void) const { return(box.isPeriodic()); }
    GCoord periodicBox(void) const { return(box.box()); }
    void periodicBox(const GCoord& c) { box.box(c); }
    void periodicBox(const greal x, const greal y, const greal z);

    //! Provide access to the underlying shared periodic box...
    //SharedPeriodicBox sharedPeriodicBox() const { return(box); }
    void reimage();
    void reimageByAtom();
    void mergeImage(pAtom &p);
    void mergeImage();
    void findBonds(const double dist = 1.65);

    std::vector<GCoord> boundingBox(void) const;
    GCoord centerAtOrigin(void);
    GCoord centroid(void) const;
    greal radius(void) const;
    GCoord centerOfMass(void) const;
    GCoord centerOfElectrons(void) const;
    GCoord dipoleMoment(void) const;
    greal totalCharge(void) const;
    greal totalMass(void) const;
    greal radiusOfGyration(void) const;
    greal rmsd(const AtomicGroup&);

    std::vector<GCoord> getTransformedCoords(const XForm&) const;
    void translate(const GCoord & v);
    void rotate(const GCoord& axis, const greal angle_in_degrees);
    void applyTransform(const XForm&);
    void copyCoordinates(AtomicGroup& g);
    void perturbCoords(const greal);
    std::vector<loos::Coord<double> > principalAxes(void) const;
    std::vector<GCoord> momentsOfInertia(void) const;
    GMatrix superposition(const AtomicGroup&);
    GMatrix alignOnto(const AtomicGroup&);
  };

  AtomicGroup operator+(const pAtom& lhs, const pAtom& rhs);
  AtomicGroup operator+(const pAtom& lhs, const AtomicGroup& rhs);

  %extend AtomicGroup {
    pAtom __getitem__(const int i) {
      return((*$self)[i]);
    }
    
    void __setitem__(const int i, const pAtom& d) {
      (*$self)[i] = d;
    }

    // Will this leak?
    char* __str__() {
      std::ostringstream oss;
      oss << *$self;
      size_t n = oss.str().size();
      char* buf = new char[n+1];
      strncpy(buf, oss.str().c_str(), n+1);
      return(buf);
    }
  };

  %rename(__add__)  loos::AtomicGroup::operator+;

};

