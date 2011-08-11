%module Coord

%{
#include <iostream>
#include <string>
#include <stdexcept>
#include "Coord.hpp"
%}


%rename(GCoord)  Coord<double>;
%rename(__add__)  Coord<double>::operator+;
%rename(__sub__) Coord<double>::operator-;
%rename(__mul__) Coord<double>::operator*;
%rename(__div__) Coord<double>::operator/;
%rename(__mod__) Coord<double>::operator%;
%rename(__pow__) Coord<double>::operator^;

// Forward declarations for matrix-vector multiply...

template<class T> class Coord;
template<class T> class Matrix44;
template<class T> Coord<T> operator*(const Matrix44<T>&, const Coord<T>&);


//! Basic 3-D coordinates class.
/**
 * Coordinates are stored internally as homogenous coordinates in an
 * array of T.  There is some awkward support for making coordinates in
 * higher-dimensions, but caveat programmer...
 *
 *
 * Warnings:
 *
 *  - The modulus operator assumes that T can be converted to double
 *    and back since it does this internally to use fmod().
 *
 *  - The distance() and length() methods assume that T can be
 *    converted into a double (and return a double as the result).
 *
 *
 *
 *Notes:
 *
 *  - The size of the stored coords is determined by the CoordIndex
 *    enum.  Places where The size is effectively hard-coded are in
 *    the constructor, the individual accessors, and in the
 *    cross-product.  Internally, the coords are homogenous, being of
 *    size n+1 where the n+1th element is always 1...
 *
 *  - Performance will probably suffer until you set the optimization
 *    level high enough to do some loop-unrolling.
 *
 *  - The +, -, and * operation are symmetric with respect to T,
 *    i.e. 1 + v and v + 1 are the same
 *
 *  - The precedence for the cross-product operator ^ is low, so be
 *    careful
 *
 *  - The * operator has dual use...  If either side is a T, then that
 *    value will be multiplied across all elements of the coord.  If
 *    it is a Coord<T>, then the dot- product will be computed.
 *    What's a little confusion amongst friends???
 */

class Coord<double> {
  enum { X=0, Y=1, Z=2, MAXCOORD } CoordIndex;

  //! The threshold for vector equality
  static const double epsilon = 1e-16;

public:



  // Constructors

  Coord<double>() { zero(); }
  Coord<double>(const double ax, const double ay, const double az) { set(ax, ay, az); }
  Coord<double>(const Coord<double>& o) { copy(o); }
  Coord<double>(const double x) { for (int i=0; i<MAXCOORD; ++i) v[i] = x; v[MAXCOORD] = 1.0; }

  // ---------------------------------------
  // Accessors

  // double& x(void);
  // const double& x(void) const;
  void x(const double ax);
  double x() const { return(v[0]); }

  // double& y(void);
  // const double& y(void) const;
  void y(const double ay);
  double y() const { return(v[1]); }

  // double& z(void);
  // const double& z(void) const;
  void z(const double az);
  double z() const { return(v[2];) }

  //! Retrieve an element from the Coord with range-checking
  //  double& operator[](const unsigned int i);

  //! Retrieve an element from a const Coord with range-checking 
  //const double& operator[](const unsigned int i) const;

  //! Short-cut to set the cartesian coordinates...
  void set(const double x, const double y, const double z);


  // ---------------------------------------
  // I/O

  //! Output the coordinate in pseudo-XML
  // friend std::ostream& operator<<(std::ostream& os, const Coord<double>&o) { 
  //   os << "(";
  //   int i;
  //   for (i=0; i<MAXCOORD; i++)
  //     os << o.v[i] << (i < MAXCOORD-1 ? "," : "");
  //   os << ")";
  //   return(os);
  // }

  // friend std::istream& operator>>(std::istream& is, Coord<double>& i) {
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




  //const Coord<double>& operator=(const Coord<double>& c);


  //! Handle addition
  Coord<double>& operator+=(const Coord<double>& rhs);

  Coord<double> operator+(const Coord<double>& rhs) const;

  //! Handle the case of double + Coord<double>
  //  friend Coord<double> operator+(const double lhs, const Coord<double>& rhs);

  //! Subtraction
  Coord<double>& operator-=(const Coord<double>& rhs);

  Coord<double> operator-(const Coord<double>& rhs) const;

  //! Unary negation
  Coord<double> operator-();
  
  //! Handle the case of T - Coord<double>
  //friend Coord<double> operator-(const double lhs, const Coord<double>& rhs);

  //! For matrix-vector multiply
  //friend Coord<double> operator*<>(const Matrix44<double>&, const Coord<double>&);


  //! Multiplication by a constant
  Coord<double>& operator*=(const double rhs);

  Coord<double> operator*(const double rhs) const;

  //! Handle double * Coord<double>
  //  friend Coord<double> operator*(const double lhs, const Coord<double>& rhs);


  //! Division by a constant
  Coord<double>& operator/=(const double rhs);

  Coord<double> operator/(const double rhs) const;

  //!  double / Coord<double> case... This may not actually be a good idea? 
  //friend Coord<double> operator/(const double lhs, const Coord<double>& rhs);

  
  //! Dot product
  double dot(const Coord<double>& rhs) const;

  double operator*(const Coord<double>rhs) const;

  //! Cross-product.  Returns a new Coord<double>
  Coord<double> cross(const Coord<double>& rhs) const;

  //! Mutating cross-product (note precedence issues)
  Coord<double>& operator^=(const Coord<double>& rhs);

  //! Cross-product (note precedence issues)
  Coord<double> operator^(const Coord<double>& rhs) const;


  //! Modulo of each component of the Coord with a constant
  Coord<double>& operator%=(const Coord<double>& rhs);

  Coord<double> operator%(const Coord<double>& rhs) const;

  //-----------------------------
  // Misc

  //! Handle coordinates with periodic boundary conditions.
  void reimage(const Coord<double>& box);

  //! Length of the Coord (as a vector) squared
  double length2(void) const;

  //! Length of the coordinate (as a vector)
  double length(void) const;
  

  //! Distance squared between two coordinates
  double distance2(const Coord<double>& o) const;

  //! Distance squared between two coordinates considering periodic
  //! boundary conditions
  double distance2(const Coord<double>& o, const Coord<double>& box) const;

  //! Distance between two coordinates.
  double distance(const Coord<double>& o) const;

  //! Distance between two coordinates considering periodic boundary
  //! conditions
  double distance(const Coord<double>& o, const Coord<double>& box) const;

  //! Zero out the coordinates (while keeping it homogenous)
  void zero(void);

  //! Compute equality based on norm(u-v) < epsilon
  bool operator==(const Coord<double>& rhs) const;

  //! Compute inequality based on ! ==
  bool operator!=(const Coord<double>& rhs) const;




private:
  void copy(const Coord<double>& c);

  double v[MAXCOORD+1];
};



%extend Coord<double> {
  double __getitem__(const int i) {
    if (i < 0 || i >= 3)
      return(0);
    return((*self)[i]);
  }

  void __setitem__(const int i, const double d) {
    if (3 >= i && i >= 0)
      (*self)[i] = d;
  }


};
