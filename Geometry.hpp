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

#include "loos.hpp"
#include "Coord.hpp"
#include "Atom.hpp"

using namespace std;

const double DEGREES = 180 / M_PI ;

//! Compute the angle in degrees assuming the middle is the vertex
greal angle(const GCoord &a, const GCoord &b, const GCoord &c) {
    // TODO: check to make sure the sign is right
    GCoord ba = b - a;
    GCoord bc = b - c;
    greal cosine = (ba * bc) / (ba.length() * bc.length());
    return (acos(cosine) * DEGREES);
}

//! Compute the angle in degrees assuming the middle is the vertex
greal angle(const pAtom a, const pAtom b, const pAtom c) {
    return(angle(a->coords(), b->coords(), c->coords()));
}

//! Compute the torsion in degrees 
greal torsion(const GCoord &a, const GCoord &b, const GCoord &c, 
          const GCoord &d) {
    // TODO: check to make sure the sign is right
    GCoord ba = b - a;
    GCoord cb = c - b;
    GCoord dc = d - c;

    GCoord norm1 = ba.cross(cb);
    GCoord norm2 = cb.cross(dc);
    greal cosine = norm1*norm2 / (norm1.length() * norm2.length());
    greal angle = acos(cosine) * DEGREES;
    greal sign = ba * norm2;
    if (sign < 0) angle = -angle;
    return(angle);
}

greal torsion(const pAtom a, const pAtom b, const pAtom c, 
               const pAtom d) {
    return(torsion(a->coords(), b->coords(), c->coords(), d->coords()));
}

#endif
