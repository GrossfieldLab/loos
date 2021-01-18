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



#if !defined(LOOS_COORDS_HPP)
#define LOOS_COORDS_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <utils_random.hpp>

namespace loos {

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

  template <class T>
  class Coord {
    enum { X=0, Y=1, Z=2, MAXCOORD } CoordIndex;

    //! The threshold for vector equality
    static const double epsilon;

  public:

    typedef T      element_type;

    // Constructors

    Coord() { zero(); }


    Coord(const T ax, const T ay, const T az) { set(ax, ay, az); }

    Coord(const Coord<T>& o) { copy(o); }

    Coord(const T x) {
      uint i;
      for (i=0; i<MAXCOORD; ++i)
        v[i] = x;
      v[i] = 1.0;
    }


    // ---------------------------------------
    // Accessors

    void x(const T ax) { v[X] = ax; }
    void y(const T ay) { v[Y] = ay; }
    void z(const T az) { v[Z] = az; }


    const T& x(void) const { return(v[X]); }
    const T& y(void) const { return(v[Y]); }
    const T& z(void) const { return(v[Z]); }


#if !defined(SWIG)
    T& x(void) { return(v[X]); }
    T& y(void) { return(v[Y]); }
    T& z(void) { return(v[Z]); }


    //! Retrieve an element from the Coord with range-checking
    T& operator[](const unsigned int i) throw(std::out_of_range) {
      if (i>=MAXCOORD)
        throw std::out_of_range("Index into Coord<T> is out of range.");
      return(v[i]);
    }

    //! Retrieve an element from a const Coord with range-checking 
    const T& operator[](const unsigned int i) const throw(std::out_of_range) {
      if (i>=MAXCOORD)
        throw std::out_of_range("Index into Coord<T> is out of range.");
      return(v[i]);
    }
#endif // !defined(SWIG)

    //! Short-cut to set the cartesian coordinates...
    void set(const T x, const T y, const T z) {
      v[X] = x;
      v[Y] = y;
      v[Z] = z;
      v[MAXCOORD] = 1;
    }


    // ---------------------------------------
    // I/O

#if !defined(SWIG)
    //! Output the coordinate in pseudo-XML
    friend std::ostream& operator<<(std::ostream& os, const Coord<T>&o) { 
      os << "(";
      for (uint i=0; i<MAXCOORD; ++i)
        os << o.v[i] << (i < MAXCOORD-1 ? "," : "");
      os << ")";
      return(os);
    }

    friend std::istream& operator>>(std::istream& is, Coord<T>& i) throw(std::runtime_error){
      char c;
      is >> c;
      if (c != '(')
        throw(std::runtime_error("Invalid Coord conversion"));
    
      is >> i.x();
      is >> c;
      if (c != ',')
        throw(std::runtime_error("Invalid Coord conversion"));

    
      is >> i.y();
      is >> c;
      if (c != ',')
        throw(std::runtime_error("Invalid Coord conversion"));

    
      is >> i.z();
      is >> c;
      if (c != ')')
        throw(std::runtime_error("Invalid Coord conversion"));

      return(is);
    }

    // ---------------------------------------
    // Operators



    const Coord<T>& operator=(const Coord<T>& c) { copy(c); return(*this); }


    //! Handle the case of T + Coord<T>
    friend Coord<T> operator+(const T lhs, const Coord<T>& rhs) {
      Coord<T> res(rhs);
      res += lhs;
      return(res);
    }


  
    //! Handle the case of T - Coord<T>
    friend Coord<T> operator-(const T lhs, const Coord<T>& rhs) {
      Coord<T> res;

      for (uint i=0; i<MAXCOORD; ++i)
        res.v[i] = lhs - rhs.v[i];

      return(res);
    }


    //! For matrix-vector multiply
    friend Coord<T> operator*<>(const Matrix44<T>&, const Coord<T>&);


    //! Handle T * Coord<T>
    friend Coord<T> operator*(const T lhs, const Coord<T>& rhs) {
      Coord<T> res(rhs);
      res *= lhs;
      return(res);
    }



    //! Division by a constant
    Coord<T>& operator/=(const T rhs) {
      for (uint i=0; i<MAXCOORD; ++i)
        v[i] /= rhs;

      return(*this);
    }

    Coord<T> operator/(const T rhs) const {
      Coord<T> res(*this);
      res /= rhs;
      return(res);
    }

    //!  T / Coord<T> case... This may not actually be a good idea? 
    friend Coord<T> operator/(const T lhs, const Coord<T>& rhs) {
      Coord<T> res;

      for (uint i=0; i<MAXCOORD; ++i)
        res.v[i] = lhs / rhs.v[i];

      return(res);
    }



#endif // !defined(SWIG)


    Coord<T>& operator+=(const T rhs) {
      for (uint i=0; i<MAXCOORD; ++i)
        v[i] += rhs;

      return(*this);
    }

    Coord<T> operator+(const T rhs) const {
      Coord<T> res(*this);
      res += rhs;
      return(res);
    }

    //! Handle addition
    Coord<T>& operator+=(const Coord<T>& rhs) {
      for (uint i=0; i<MAXCOORD; i++)
        v[i] += rhs.v[i];

      return(*this);
    }

    Coord<T> operator+(const Coord<T>& rhs) const {
      Coord<T> res(*this);
      res += rhs;
      return(res);
    }

    //! Subtraction
    Coord<T>& operator-=(const Coord<T>& rhs) {
      for (uint i=0; i<MAXCOORD; ++i)
        v[i] -= rhs.v[i];

      return(*this);
    }

    Coord<T> operator-(const Coord<T>& rhs) const {
      Coord<T> res(*this);
      res -= rhs;
      return(res);
    }

    Coord<T>& operator-=(const T rhs) {
      for (uint i=0; i<MAXCOORD; ++i)
        v[i] -= rhs;
      
      return(*this);
    }

    Coord<T> operator-(const T rhs) const {
      Coord<T> res(*this);

      res -= rhs;
      return(res);
    }

    //! Unary negation
    Coord<T> operator-() {
      Coord<T> res(*this);

      for (uint i=0; i<MAXCOORD; ++i)
        res.v[i] = -res.v[i];
      return(res);
    }



    //! Multiplication by a constant
    Coord<T>& operator*=(const T rhs) {

      for (uint i=0; i<MAXCOORD; ++i)
        v[i] *= rhs;

      return(*this);
    }

    Coord<T> operator*(const T rhs) const {
      Coord<T> res(*this);
      res *= rhs;
      return(res);
    }

  
    //! Dot product
    T dot(const Coord<T>& rhs) const {
      T s = v[0] * rhs.v[0];
      for (uint i=1; i<MAXCOORD; ++i)
        s += v[i] * rhs.v[i];

      return(s);
    }

    //! Unit vector dot product (cosine of projective angle)
    T uvdot(const Coord<T>& rhs) const {
      // define quantities to accumulate with first step of loop.
      T s = v[0] * rhs.v[0];
      T d = v[0] * v[0];
      T drhs = rhs.v[0] * rhs.v[0];
      
      for (uint i = 1; i < MAXCOORD; ++i){
        // dot product
        s += v[i] * rhs.v[i];
        // length2 of this
        d += v[i] * v[i];
        // length2 of rhs
        drhs += rhs.v[i] * rhs.v[i];
      }
      return(s / sqrt(d * drhs));
    }   



    T operator*(const Coord<T>rhs) const {
      return(dot(rhs));
    }


    //! Cross-product.  Returns a new Coord<T>
    Coord<T> cross(const Coord<T>& rhs) const {
      Coord<T> res;
      res.v[X] = v[Y] * rhs.v[Z] - v[Z] * rhs.v[Y];
      res.v[Y] = v[Z] * rhs.v[X] - v[X] * rhs.v[Z];
      res.v[Z] = v[X] * rhs.v[Y] - v[Y] * rhs.v[X];

      return(res);
    }

    //! Mutating cross-product (note precedence issues)
    Coord<T>& operator^=(const Coord<T>& rhs) {
      Coord<T> tmp = cross(rhs);
      *this = tmp;
      return(*this);
    }

    //! Cross-product (note precedence issues)
    Coord<T> operator^(const Coord<T>& rhs) const {
      return(cross(rhs));
    }


    //! Modulo of each component of the Coord with a constant
    Coord<T>& operator%=(const Coord<T>& rhs) {

      for (uint i=0; i<MAXCOORD; ++i)
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

    //! Handle coordinates with periodic boundary conditions.
    void reimage(const Coord<T>& box) {
      for (uint i=0; i<MAXCOORD; ++i) {
        int n = (int)(fabs(v[i]) / box.v[i] + 0.5);
        v[i] = (v[i] >= 0) ? v[i] - n*(box.v[i]) : v[i] + n*(box.v[i]);
      }
    }


    //! Length of the Coord (as a vector) squared
    double length2(void) const {
      double d = 0.0;

      for (uint i=0; i<MAXCOORD; ++i)
        d += v[i]*v[i];

      return(d);
    }

    //! Length of the coordinate (as a vector)
    double length(void) const {
      return(sqrt(length2()));
    }
  

    //! Distance squared between two coordinates
    double distance2(const Coord<T>& o) const {
      Coord<T> d = o - *this;
      return(d.length2());
    }

    //! Distance squared between two coordinates considering periodic
    //! boundary conditions
    double distance2(const Coord<T>& o, const Coord<T>& box) const {
      Coord<T> d = o - *this;
      d.reimage(box);
      return(d.length2());
    }

    //! Distance between two coordinates.
    double distance(const Coord<T>& o) const {
      return(sqrt(distance2(o)));
    }

    //! Distance between two coordinates considering periodic boundary
    //! conditions
    double distance(const Coord<T>& o, const Coord<T>& box) const {
      return(sqrt(distance2(o, box)));
    }

    //! Generate a random vector on a unit sphere
    /**
     *  Note: this method uses a singleton random number generator set up
     *        by LOOS, but doesn't set the initial seed.  Instead, the calling
     *        program needs to do that, either by calling randomSeedRNG()
     *        (see utils_random.cpp), which seeds off the time, or by generating 
     *        the seed another way.
     */
    void random(void) {
        base_generator_type& rng = rng_singleton();
        boost::uniform_on_sphere<> uni(3);
        boost::variate_generator<base_generator_type&, 
                                 boost::uniform_on_sphere<> > func(rng, uni);
        std::vector<double> vec = func();
        for (uint i=0; i<MAXCOORD; i++) {
            v[i] = static_cast<T>(vec[i]);
        }

    }
    


    //! Zero out the coordinates (while keeping it homogenous)
    void zero(void) {
      uint i;
      for (i=0; i<MAXCOORD; ++i)
        v[i] = 0;
      v[i] = 1;
    }

    //! Compute equality based on norm(u-v) < epsilon
    bool operator==(const Coord<T>& rhs) const {
      double d = distance2(rhs);
      if (d < epsilon * epsilon)
        return(true);
      return(false);
    }

    //! Compute inequality based on ! ==
    bool operator!=(const Coord<T>& rhs) const {
      return(!(operator==(rhs)));
    }




  private:
    void copy(const Coord<T>& c) {
      if (&c == this)
        return;

      uint i;
      for (i=0; i<MAXCOORD; ++i)
        v[i] = c.v[i];

      v[i] = 1;
    };
  

    T v[MAXCOORD+1];
  };

  template<typename T>
  const double Coord<T>::epsilon = 1e-15;

}

#endif
