/*
  Fmt.hpp
  Tod D. Romo

  Formatter based heavily on example from Stroustrup's C++ book...
*/




#if !defined(LOOSFMT_HPP)
#define LOOSFMT_HPP

#include <ios>
#include <sstream>
#include <iostream>

namespace loos {

  class BoundFmt;

  //! Output formatter class, adapted from Stroustrup's book
  class Fmt {
    friend std::ostream& operator<<(std::ostream&, const BoundFmt&);

  public:
    //! Alignment of the text
    enum FmtAlignment { LEFT, RIGHT, INTERNAL };

    //! Default is for precision width 6, no zeros, padding with spaces,
    //! and left aligned (and general formatting)
    explicit Fmt(int p = 6) : prc(p), wdth(0), fil(' '), trl(false), pos(false), ali(LEFT) {
      fmt = ~(std::ios_base::fixed|std::ios_base::scientific);
    }

    //! Returns the bound formatter
    BoundFmt operator()(double d) const;
  
    //! Output in scientific format
    Fmt& scientific() { fmt = std::ios_base::scientific; return(*this); }
    //! Output in fixed-point
    Fmt& fixed() { fmt = std::ios_base::fixed; return(*this); }
    //! Output normally (i.e. general)
    Fmt& general() { fmt = ~(std::ios_base::fixed|std::ios_base::scientific); return(*this); }

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
    std::ios_base::fmtflags fmt;
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

  std::ostream& operator<<(std::ostream&, const BoundFmt&);


}

#endif
