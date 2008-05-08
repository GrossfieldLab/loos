/*
  Coords.h
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Basic 3-D coordinates class.



  Warnings:

    o The modulus operator assumes that T can be converted to double
      and back since it does this internally to use fmod().

    o The distance() and length() methods assume that T can be
      converted into a double (and return a double as the result).



  Notes:

    o The size of the stored coords is determined by the CoordIndex
      enum.  Places where The size is effectively hard-coded are in
      the constructor, the individual accessors, and in the
      cross-product.  Internally, the coords are homogenous, being of
      size n+1 where the n+1th element is always 1...

    o Performance will probably suffer until you set the optimization
      level high enough to do some loop-unrolling.

    o The +, -, and * operation are symmetric with respect to T,
      i.e. 1 + v and v + 1 are the same

    o The precedence for the cross-product operator ^ is low, so be
      careful

    o The * operator has dual use...  If either side is a T, then that
      value will be multiplied across all elements of the coord.  If
      it is a Coord<T>, then the dot- product will be computed.
      What's a little confusion amongst friends???

*/


#if !defined(COORDS_HPP)
#define COORDS_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include <math.h>

using namespace std;

// Forward declarations for matrix-vector multiply...

template<class T> class Coord;
template<class T> class Matrix44;
template<class T> Coord<T> operator*(const Matrix44<T>&, const Coord<T>&);

template <class T>
class Coord {
  enum { X=0, Y=1, Z=2, MAXCOORD } CoordIndex;
public:



  // Constructors

  Coord() { zero(); }


  Coord(const T ax, const T ay, const T az) { set(ax, ay, az); }

  Coord(const Coord<T>& o) { copy(o); }

  Coord(const T x) {
    int i;
    for (i=0; i<MAXCOORD; i++)
      v[i] = x;
    v[i] = 1.0;
  }


  // ---------------------------------------
  // Accessors

  T x(void) const { return(v[X]); }
  void x(const T ax) { v[X] = ax; }

  T y(void) const { return(v[Y]); }
  void y(const T ay) { v[Y] = ay; }

  T z(void) const { return(v[Z]); }
  void z(const T az) { v[Z] = az; }


  T& operator[](const unsigned int i) {
    if (i>=MAXCOORD)
      throw out_of_range("Index into Coord<T> is out of range.");
    return(v[i]);
  }

  const T& operator[](const unsigned int i) const {
    if (i>=MAXCOORD)
      throw out_of_range("Index into Coord<T> is out of range.");
    return(v[i]);
  }

  void set(const T x, const T y, const T z) {
    v[X] = x;
    v[Y] = y;
    v[Z] = z;
    v[MAXCOORD] = 1;
  }


  // ---------------------------------------
  // I/O

  friend ostream& operator<<(ostream& os, const Coord<T>&o) { 
    os << "(";
    int i;
    for (i=0; i<MAXCOORD; i++)
      os << o.v[i] << (i < MAXCOORD-1 ? "," : "");
    os << ")";
    return(os);
  }

  // ---------------------------------------
  // Operators




  const Coord<T>& operator=(const Coord<T>& c) { copy(c); return(*this); }


  // Handle addition
  Coord<T>& operator+=(const Coord<T>& rhs) {
    int i;

    for (i=0; i<MAXCOORD; i++)
      v[i] += rhs.v[i];

    return(*this);
  }

  Coord<T> operator+(const Coord<T>& rhs) const {
    Coord<T> res(*this);
    res += rhs;
    return(res);
  }

  // Handle the case of T + Coord<T>
  friend Coord<T> operator+(const T lhs, const Coord<T>& rhs) {
    Coord<T> res(rhs);
    res += lhs;
    return(res);
  }


  // Subtraction
  Coord<T>& operator-=(const Coord<T>& rhs) {
    int i;

    for (i=0; i<MAXCOORD; i++)
      v[i] -= rhs.v[i];

    return(*this);
  }

  Coord<T> operator-(const Coord<T>& rhs) const {
    Coord<T> res(*this);
    res -= rhs;
    return(res);
  }

  // Unary negation
  Coord<T> operator-() {
    Coord<T> res(*this);
    int i;

    for (i=0; i<MAXCOORD; i++)
      res.v[i] = -res.v[i];
    return(res);
  }
  
  // Handle the case of T - Coord<T>
  friend Coord<T> operator-(const T lhs, const Coord<T>& rhs) {
    int i;
    Coord<T> res;

    for (i=0; i<MAXCOORD; i++)
      res.v[i] = lhs - rhs.v[i];

    return(res);
  }


  // For matrix-vector multiply
  friend Coord<T> operator*<>(const Matrix44<T>&, const Coord<T>&);


  // Multiplication by a constant
  Coord<T>& operator*=(const T rhs) {
    int i;

    for (i=0; i<MAXCOORD; i++)
      v[i] *= rhs;

    return(*this);
  }

  Coord<T> operator*(const T rhs) const {
    Coord<T> res(*this);
    res *= rhs;
    return(res);
  }

  // And handle T * Coord<T>
  friend Coord<T> operator*(const T lhs, const Coord<T>& rhs) {
    Coord<T> res(rhs);
    res *= lhs;
    return(res);
  }



  // Division by a constant
  Coord<T>& operator/=(const T rhs) {
    int i;

    for (i=0; i<MAXCOORD; i++)
      v[i] /= rhs;

    return(*this);
  }

  Coord<T> operator/(const T rhs) const {
    Coord<T> res(*this);
    res /= rhs;
    return(res);
  }

  // And the T / Coord<T> case...
  // This may not actually be a good idea...
  friend Coord<T> operator/(const T lhs, const Coord<T>& rhs) {
    Coord<T> res;
    int i;

    for (i=0; i<MAXCOORD; i++)
      res.v[i] = lhs / rhs.v[i];

    return(res);
  }

  
  // Dot product
  T dot(const Coord<T>& rhs) const {
    int i;

    T s = v[0] * rhs.v[0];
    for (i=1; i<MAXCOORD; i++)
      s += v[i] * rhs.v[i];

    return(s);
  }


  T operator*(const Coord<T>rhs) const {
    return(dot(rhs));
  }

  // Cross-product
  // Note precedence issues...
  void cross(const Coord<T>& rhs) {
    Coord<T> lhs(*this);
    v[X] = lhs.v[Y] * rhs.v[Z] - lhs.v[Z] * rhs.v[Y];
    v[Y] = lhs.v[Z] * rhs.v[X] - lhs.v[X] * rhs.v[Z];
    v[Z] = lhs.v[X] * rhs.v[Y] - lhs.v[Y] * rhs.v[X];
  }
  
  Coord<T>& operator^=(const Coord<T>& rhs) {
    cross(rhs);
    return(*this);
  }

  Coord<T> operator^(const Coord<T>& rhs) const {
    Coord<T> res(*this);
    res ^= rhs;
    return(res);
  }


  // Modulo
  Coord<T>& operator%=(const Coord<T>& rhs) {
    int i;

    for (i=0; i<MAXCOORD; i++)
      v[i] = fmod(v[i], rhs.v[i]);

    return(*this);
  }

  Coord<T> operator%(const Coord<T>& rhs) const {
    Coord<T> res(*this);

    res %= rhs;
    return(res);
  }



  //-----------------------------
  // Misc


  void reimage(const Coord<T>& box) {
    int i;
    int n;
    
    for (i=0; i<MAXCOORD; i++) {
      n = (int)(fabs(v[i]) / box.v[i] + 0.5);
      v[i] = (v[i] >= 0) ? v[i] - n*(box.v[i]) : v[i] + n*(box.v[i]);
    }
  }


  void canonical(const Coord<T>& box) {
    int i;
    for (i=0; i<MAXCOORD; i++)
      v[i] = fmod(v[i] + 1.5 * box.v[i], box.v[i]) - (box.v[i] / 2.0);
  }


  double length2(void) const {
    double d = 0.0;
    int i;

    for (i=0; i<MAXCOORD; i++)
      d += v[i]*v[i];

    return(d);
  }


  double length(void) const {
    return(sqrt(length2()));
  }
  

  double distance2(const Coord<T>& o) {
    Coord<T> d = o - *this;
    return(d.length2());
  }

  double distance2(const Coord<T>& o, const Coord<T>& box) {
    Coord<T> d = o - *this;
    d.reimage(box);
    return(d.length2());
  }

  double distance(const Coord<T>& o) {
    return(sqrt(distance2(o)));
  }

  double distance(const Coord<T>& o, const Coord<T>& box) {
    return(sqrt(distance2(o, box)));
  }



  void zero(void) {
    int i;
    for (i=0; i<MAXCOORD; i++)
      v[i] = 0;
    v[i] = 1.0;
  }


private:
  void copy(const Coord<T>& c) {
    if (&c == this)
      return;

    int i;
    for (i=0; i<MAXCOORD; i++)
      v[i] = c.v[i];

    v[i] = 1.0;
  };
  

  T v[MAXCOORD+1];
};



#endif
