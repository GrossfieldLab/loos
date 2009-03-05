/*
  Fmt.cpp

  Formatter based heavily on example from Stroustrup's C++ book...
*/



#include <Fmt.hpp>


namespace loos {
  Fmt& Fmt::scientific() { fmt = std::ios_base::scientific; return(*this); }
    Fmt& Fmt::fixed() { fmt = std::ios_base::fixed; return(*this); }
    Fmt& Fmt::general() { fmt = ~(std::ios_base::fixed|std::ios_base::scientific); return(*this); }

    Fmt& Fmt::precision(const int p) { prc = p; return(*this); }
    Fmt& Fmt::width(const int w) { wdth = w; return(*this); }
    Fmt& Fmt::fill(const char c) { fil = c; return(*this); }

    Fmt& Fmt::trailingZeros(bool b) { trl = b; return(*this); }
    Fmt& Fmt::plus(bool b) { pos = b; return(*this); }

    Fmt& Fmt::left(void) { ali = LEFT; return(*this); }
    Fmt& Fmt::right(void) { ali = RIGHT; return(*this); }
    Fmt& Fmt::internal(void) { ali = INTERNAL; return(*this); }


  //! Return the bound formatter
  BoundFmt Fmt::operator()(double d) const { return(BoundFmt(*this, d)); }

  //! Create the output with the specified formatter
  std::ostream& operator<<(std::ostream& os, const BoundFmt& bf) {
    std::ostringstream s;

    s.precision(bf.f.prc);
    s.setf(bf.f.fmt, std::ios_base::floatfield);
    s.width(bf.f.wdth);
    s.fill(bf.f.fil);

    switch(bf.f.ali) {
    case Fmt::LEFT:
      s.setf(std::ios_base::left, std::ios_base::adjustfield); break;
    case Fmt::RIGHT:
      s.setf(std::ios_base::right, std::ios_base::adjustfield); break;
    default:
      s.setf(std::ios_base::internal, std::ios_base::adjustfield); break;
    }

    if (bf.f.trl)
      s.setf(std::ios_base::showpoint);
    else
      s.unsetf(std::ios_base::showpoint);

    if (bf.f.pos)
      s.setf(std::ios_base::showpos);
    else
      s.unsetf(std::ios_base::showpos);

    s << bf.val;
    return(os << s.str());
  }


}
