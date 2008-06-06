/*
  Geometry.hpp
  (c) 2008 Alan Grossfield


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/


#if !defined(GEOM_HPP)
#define GEOM_HPP

#include <iostream>
#include <string>
#include <stdexcept>

#include <math.h>

#include "Coord.hpp"
#include "Atom.hpp"

using namespace std;

const double DEGREES = 180 / M_PI ;

// TODO: do 2 versions, taking Coord and const pAtom types


//! Compute the angle in degrees assuming the middle is the vertex
template<class T> 
T angle(const Coord<T> a, const Coord<T> b, const Coord<T> c) {
    // TODO: check to make sure the sign is right
    Coord<T> ba = b - a;
    Coord<T> bc = b - c;
    T cosine = (ba * bc) / (ba.length() * bc.length());
    return (acos(cosine) * DEGREES);
}

//! Compute the angle in degrees assuming the middle is the vertex
double atom_angle(const pAtom & a, const pAtom & b, const pAtom & c) {
    return(angle(a->coords(), b->coords(), c->coords()));
}

//! Compute the torsion in degrees 
template<class T> 
T torsion(const Coord<T> a, const Coord<T> b, const Coord<T> c, 
          const Coord<T> d) {
    // TODO: check to make sure the sign is right
    Coord<T> ba = b - a;
    Coord<T> cb = c - b;
    Coord<T> dc = d - c;

    Coord<T> norm1 = ba^cb;
    Coord<T> norm2 = cb^dc;
    T cosine = norm1*norm2 / (norm1.length() * norm2.length());
    T angle = acos(cosine) * DEGREES;
    T sign = ba * norm2;
    if (sign < 0) angle = -angle;
    return(angle);
}

double atom_torsion(const pAtom & a, const pAtom & b, const pAtom & c, 
               const pAtom & d) {
    return(torsion(a->coords(), b->coords(), c->coords(), d->coords()));
}

#endif
