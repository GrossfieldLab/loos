%module Coord

%{
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <Coord.hpp>
%}

namespace loos {

  template<class T>
  class Coord {
    enum { X=0, Y=1, Z=2, MAXCOORD } CoordIndex;

    //! The threshold for vector equality
    static const double epsilon = 1e-16;

  public:



    // Constructors

    Coord() { zero(); }
    Coord(const T ax, const T ay, const T az) { set(ax, ay, az); }
    Coord(const Coord<T>& o) { copy(o); }
    Coord(const T x) { for (int i=0; i<MAXCOORD; ++i) v[i] = x; v[MAXCOORD] = 1.0; }

    // ---------------------------------------
    // Accessors

    // T& x(void);
    // const T& x(void) const;
    void x(const T ax);
    T x() const { return(v[0]); }

    // T& y(void);
    // const T& y(void) const;
    void y(const T ay);
    T y() const { return(v[1]); }

    // T& z(void);
    // const T& z(void) const;
    void z(const T az);
    T z() const { return(v[2];) }

    //! Retrieve an element from the Coord with range-checking
    //  T& operator[](const unsigned int i);

    //! Retrieve an element from a const Coord with range-checking 
    //const T& operator[](const unsigned int i) const;

    //! Short-cut to set the cartesian coordinates...
    void set(const T x, const T y, const T z);


    // ---------------------------------------
    // I/O

    //! Output the coordinate in pseudo-XML
    // friend std::ostream& operator<<(std::ostream& os, const Coord<T>&o) { 
    //   os << "(";
    //   int i;
    //   for (i=0; i<MAXCOORD; i++)
    //     os << o.v[i] << (i < MAXCOORD-1 ? "," : "");
    //   os << ")";
    //   return(os);
    // }

    // friend std::istream& operator>>(std::istream& is, Coord<T>& i) {
    //   char c;
    //   is >> c;
    //   if (c != '(')
    //     throw(std::runtime_error("Invalid Coord conversion"));
    
    //   is >> i.x();
    //   is >> c;
    //   if (c != ',')
    //     throw(std::runtime_error("Invalid Coord conversion"));

    
    //   is >> i.y();
    //   is >> c;
    //   if (c != ',')
    //     throw(std::runtime_error("Invalid Coord conversion"));

    
    //   is >> i.z();
    //   is >> c;
    //   if (c != ')')
    //     throw(std::runtime_error("Invalid Coord conversion"));

    //   return(is);
    // }

    // ---------------------------------------
    // Operators




    //const Coord<T>& operator=(const Coord<T>& c);


    //! Handle addition
    Coord<T>& operator+=(const Coord<T>& rhs);

    Coord<T> operator+(const Coord<T>& rhs) const;

    //! Handle the case of T + Coord<T>
    //  friend Coord<T> operator+(const T lhs, const Coord<T>& rhs);

    //! Subtraction
    Coord<T>& operator-=(const Coord<T>& rhs);

    Coord<T> operator-(const Coord<T>& rhs) const;

    //! Unary negation
    Coord<T> operator-();
  
    //! Handle the case of T - Coord<T>
    //friend Coord<T> operator-(const T lhs, const Coord<T>& rhs);

    //! For matrix-vector multiply
    //friend Coord<T> operator*<>(const Matrix44<T>&, const Coord<T>&);


    //! Multiplication by a constant
    Coord<T>& operator*=(const T rhs);

    Coord<T> operator*(const T rhs) const;

    //! Handle T * Coord<T>
    //  friend Coord<T> operator*(const T lhs, const Coord<T>& rhs);


    //! Division by a constant
    Coord<T>& operator/=(const T rhs);

    Coord<T> operator/(const T rhs) const;

    //!  T / Coord<T> case... This may not actually be a good idea? 
    //friend Coord<T> operator/(const T lhs, const Coord<T>& rhs);

  
    //! Dot product
    T dot(const Coord<T>& rhs) const;

    T operator*(const Coord<T>rhs) const;

    //! Cross-product.  Returns a new Coord<T>
    Coord<T> cross(const Coord<T>& rhs) const;

    //! Mutating cross-product (note precedence issues)
    Coord<T>& operator^=(const Coord<T>& rhs);

    //! Cross-product (note precedence issues)
    Coord<T> operator^(const Coord<T>& rhs) const;


    //! Modulo of each component of the Coord with a constant
    Coord<T>& operator%=(const Coord<T>& rhs);

    Coord<T> operator%(const Coord<T>& rhs) const;

    //-----------------------------
    // Misc

    //! Handle coordinates with periodic boundary conditions.
    void reimage(const Coord<T>& box);

    //! Length of the Coord (as a vector) squared
    T length2(void) const;

    //! Length of the coordinate (as a vector)
    T length(void) const;
  

    //! Distance squared between two coordinates
    T distance2(const Coord<T>& o) const;

    //! Distance squared between two coordinates considering periodic
    //! boundary conditions
    T distance2(const Coord<T>& o, const Coord<T>& box) const;

    //! Distance between two coordinates.
    T distance(const Coord<T>& o) const;

    //! Distance between two coordinates considering periodic boundary
    //! conditions
    T distance(const Coord<T>& o, const Coord<T>& box) const;

    //! Zero out the coordinates (while keeping it homogenous)
    void zero(void);

    //! Compute equality based on norm(u-v) < epsilon
    bool operator==(const Coord<T>& rhs) const;

    //! Compute inequality based on ! ==
    bool operator!=(const Coord<T>& rhs) const;




  private:
    void copy(const Coord<T>& c);

    T v[MAXCOORD+1];
  };



  %extend Coord<T> {
    T __getitem__(const int i) {
      if (i < 0 || i >= 3)
        return(0);
      return((*$self)[i]);
    }
    
    void __setitem__(const int i, const T d) {
      if (3 >= i && i >= 0)
        (*$self)[i] = d;
    }

  };

};


%extend loos::Coord<double> {
  char* __str__() {
    static char buf[1024];
    std::ostringstream oss;
    oss << *$self;
    strncpy(buf, oss.str().c_str(), sizeof(buf));
    return(buf);
  }

 };

%template(GCoord)  loos::Coord<double>;
%rename(__add__)  loos::Coord<double>::operator+;
%rename(__sub__) loos::Coord<double>::operator-;
%rename(__mul__) loos::Coord<double>::operator*;
%rename(__div__) loos::Coord<double>::operator/;
%rename(__mod__) loos::Coord<double>::operator%;
%rename(__pow__) loos::Coord<double>::operator^;
