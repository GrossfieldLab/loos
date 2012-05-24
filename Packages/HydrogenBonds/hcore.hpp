/*
  Core code for hbonds utilties
*/



/*

  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2010, Tod D. Romo
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

#if !defined(LOOS_HCORE_HPP)
#define LOOS_HCORE_HPP

#include <loos.hpp>



typedef loos::Math::Matrix<int, loos::Math::RowMajor>   BondMatrix;


// Our own exception so we can provide a little more helpful
// information when we throw-up...

struct ErrorWithAtom : public std::exception {
  std::string _msg;
  
  ErrorWithAtom(const loos::pAtom& a, const std::string& msg) {
    std::stringstream ss;
    ss << msg << std::endl << *a << std::endl;
    _msg = ss.str();
  }

  const char* what() const throw() { return(_msg.c_str()); }

  ~ErrorWithAtom() throw() { }
};



// Track the atoms that may participate in a hydrogen bond, and the
// atoms they're attached to (if necessary) to compute the bond angle.
// Also encapsulates the operations for determining if an h-bond
// exists...
//
// Note that we hook into the parent group's SharedPeriodicBox so we
// always have current periodic boundary info...


class SimpleAtom {
public:
  SimpleAtom(const loos::pAtom& a) : atom(a), isHydrogen(false), usePeriodicity(false) { }
  SimpleAtom(const loos::pAtom& a, const loos::SharedPeriodicBox& b, const bool c = true) : atom(a), isHydrogen(false), usePeriodicity(c), sbox(b) { }

  void attach(const loos::pAtom&a) { attached_to = a; }
  loos::pAtom attachedTo() const { return(attached_to); }

  double distance2(const SimpleAtom& s) const;
  double angle(const SimpleAtom& s) const;

  static bool debuggingMode() { return(debugging); }
  static void debuggingMode(const bool b) { debugging = b; }

  static double innerRadius()  { return(sqrt(inner)); }
  static void innerRadius(const double r)  { inner = r*r; }

  static double outerRadius()  { return(sqrt(outer)); }
  static void outerRadius(const double r)  { outer = r*r; }

  static double maxDeviation()  { return(deviation); }
  static void maxDeviation(const double d)  { deviation = d; }


  // Tests whether two SimpleAtoms have a potential hydrogen-bond
  // between them.
  bool hydrogenBond(const SimpleAtom& other) const;

  // Returns an AtomicGroup of all atoms that may have hydrogen-bonds
  // to the current SimpleAtom
  loos::AtomicGroup findHydrogenBonds(const std::vector<SimpleAtom>& group, const bool findFirstOnly = true);

  std::vector<uint> findHydrogenBondsVector(const std::vector<SimpleAtom>& group);
  
  // Returns a matrix where the rows represent time (frames in the
  // trajectory) and columns represent acceptors (i.e. the passed
  // group).  Wherever there is a hydrogen-bond, U_ij is 1, and 0
  // otherwise.
  //
  // maxt determines the maximum time (frame #) that is considered.
  BondMatrix findHydrogenBondsMatrix(const std::vector<SimpleAtom>& group, loos::pTraj& traj, loos::AtomicGroup& model, const uint maxt) const;
  BondMatrix findHydrogenBondsMatrix(const std::vector<SimpleAtom>& group, loos::pTraj& traj, loos::AtomicGroup& model) const {
    return(findHydrogenBondsMatrix(group, traj, model, traj->nframes()));
  }


  // Converts an AtomicGroup into a vector of SimpleAtom's based on
  // the passed selection.  The use_periodicity is applied to all
  // created SimpleAtoms...they also shared the PeriodicBox with the
  // passed AtomicGroup.

  static std::vector<SimpleAtom> processSelection(const std::string& selection, const loos::AtomicGroup& system, const bool use_periodicity = false);

  friend std::ostream& operator<<(std::ostream& os, const SimpleAtom& s) {
    os << "<SimpleAtom>\n";
    os << *(s.atom) << std::endl;
    os << "<isHydrogen " << s.isHydrogen << "/>\n";
    os << "<usePeriodicity " << s.usePeriodicity << "/>\n";
    if (s.usePeriodicity)
      os << "<PeriodicBox>" << s.sbox.box() << "</PeriodicBox>\n";
    if (s.attached_to != 0)
      os << *(s.attached_to) << std::endl;
    os << "</SimpleAtom>";
    return(os);
  }


private:

  loos::pAtom atom;
  bool isHydrogen;
  bool usePeriodicity;

  static double inner, outer, deviation;
  static bool debugging;

  loos::SharedPeriodicBox sbox;
  loos::pAtom attached_to;

};



// Typedefs to make life easier in the tools...

typedef SimpleAtom    SAtom;
typedef std::vector<SAtom> SAGroup;


#endif
