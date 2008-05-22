/*
  Fmt.hpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Formatter routines, adapted from Stroustrup's book...
*/


#if !defined(FMT_HPP)
#define FMT_HPP

#include <ios>
#include <sstream>
#include <iostream>

using namespace std;


class BoundFmt;

//! Output formatter class, adapted from Stroustrup's book
class Fmt {
  friend ostream& operator<<(ostream&, const BoundFmt&);

public:
  //! Alignment of the text
  enum FmtAlignment { LEFT, RIGHT, INTERNAL };

  //! Default is for precision width 6, no zeros, padding with spaces,
  //! and left aligned (and general formatting)
  explicit Fmt(int p = 6) : prc(p), wdth(0), fil(' '), trl(false), pos(false), ali(LEFT) {
    fmt = ~(ios_base::fixed|ios_base::scientific);
  }

  //! Returns the bound formatter
  BoundFmt operator()(double d) const;
  
  //! Output in scientific format
  Fmt& scientific() { fmt = ios_base::scientific; return(*this); }
  //! Output in fixed-point
  Fmt& fixed() { fmt = ios_base::fixed; return(*this); }
  //! Output normally (i.e. general)
  Fmt& general() { fmt = ~(ios_base::fixed|ios_base::scientific); return(*this); }

  //! Set the precision
  Fmt& precision(const int p) { prc = p; return(*this); }
  //! Set the output field width
  Fmt& width(const int w) { wdth = w; return(*this); }
  //! Set the fill character
  Fmt& fill(const char c) { fil = c; return(*this); }

  //! Determines whether or not trailing zeros are shown
  Fmt& trailingZeros(bool b=true) { trl = b; return(*this); }
  //! Prepend plus sign?
  Fmt& plus(bool b=true) { pos = b; return(*this); }

  //! Align left
  Fmt& left(void) { ali = LEFT; return(*this); }
  //! Align right
  Fmt& right(void) { ali = RIGHT; return(*this); }
  //! Align "internal" (see C++ ref)
  Fmt& internal(void) { ali = INTERNAL; return(*this); }

private:

private:
  ios_base::fmtflags fmt;
  int prc, wdth;
  char fil;
  bool trl, pos;
  FmtAlignment ali;

};


//! Internal helper class to bind formatting state
struct BoundFmt {
  BoundFmt(const Fmt& ff, double v) : f(ff), val(v) { }
  const Fmt& f;
  double val;
};

//! Return the bound formatter
BoundFmt Fmt::operator()(double d) const { return(BoundFmt(*this, d)); }

//! Create the output with the specified formatter
ostream& operator<<(ostream& os, const BoundFmt& bf) {
  ostringstream s;

  s.precision(bf.f.prc);
  s.setf(bf.f.fmt, ios_base::floatfield);
  s.width(bf.f.wdth);
  s.fill(bf.f.fil);

  switch(bf.f.ali) {
  case Fmt::LEFT:
    s.setf(ios_base::left, ios_base::adjustfield); break;
  case Fmt::RIGHT:
    s.setf(ios_base::right, ios_base::adjustfield); break;
  default:
    s.setf(ios_base::internal, ios_base::adjustfield); break;
  }

  if (bf.f.trl)
    s.setf(ios_base::showpoint);
  else
    s.unsetf(ios_base::showpoint);

  if (bf.f.pos)
    s.setf(ios_base::showpos);
  else
    s.unsetf(ios_base::showpos);

  s << bf.val;
  return(os << s.str());
}


#endif
