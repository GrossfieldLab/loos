/*
  Fmt.cpp
  (c) 2008 Tod D. Romo

  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Formatter routines, adapted from Stroustrup's book...
*/


#include <Fmt.hpp>


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

