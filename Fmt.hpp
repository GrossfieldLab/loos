/*
  Fmt.hpp
  Tod D. Romo

  Formatter based heavily on example from Stroustrup's C++ book...
*/




#if !defined(LOOS_FMT_HPP)
#define LOOS_FMT_HPP

#include <ios>
#include <sstream>
#include <iostream>

#include <loos_defs.hpp>

namespace loos {

  struct BoundFmt;

  //! Output formatter class, adapted from Stroustrup's book
  class Fmt {
    friend std::ostream& operator<<(std::ostream&, const BoundFmt&);

  public:
    //! Alignment of the text
    enum FmtAlignment { LEFT, RIGHT, INTERNAL };

    //! Default is for precision width 6, no zeros, padding with spaces,
    //! and left aligned (and general formatting)
    explicit Fmt(uint p = 6) : prc(p), wdth(0), fil(' '), trl(false), pos(false), ali(LEFT) {
      fmt = ~(std::ios_base::fixed|std::ios_base::scientific);
    }

    //! Returns the bound formatter
    BoundFmt operator()(double d) const;
  
    //! Output in scientific format
    Fmt& scientific(void);
    //! Output in fixed-point
    Fmt& fixed(void);
    //! Output normally (i.e. general)
    Fmt& general(void);

    //! Set the precision
    Fmt& precision(const uint);
    //! Set the output field width
    Fmt& width(const uint);
    //! Set the fill character
    Fmt& fill(const char);

    //! Determines whether or not trailing zeros are shown
    Fmt& trailingZeros(bool b = true);
    //! Prepend plus sign?
    Fmt& plus(bool b = true);

    //! Align left
    Fmt& left(void);
    //! Align right
    Fmt& right(void);
    //! Align "internal" (see C++ ref)
    Fmt& internal(void);

  private:

  private:
    std::ios_base::fmtflags fmt;
    uint prc, wdth;
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
