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





#if !defined(PERIODICBOX_HPP)
#define PERIODICBOX_HPP



#include <boost/shared_ptr.hpp>

#include <loos.hpp>

//! Class for managing periodic box information.
/** This is the fundamental object that gets shared amongst related
 *  groups.  It contains the GCoord representing the box size and a
 *  flag that indicates whether or not the box has actually been set.
 *  The client will not interact with this class/object directly, but
 *  will use the SharedPeriodicBox instead.
 */

class PeriodicBox {
public:
  PeriodicBox() : thebox(99999,99999,99999), box_set(false) { }
  PeriodicBox(const GCoord& c) : thebox(c), box_set(true) { }

  GCoord box(void) const { return(thebox); }
  void box(const GCoord& c) { thebox = c; box_set = true; }

  bool isPeriodic(void) const { return(box_set); }
  void setPeriodic(const bool b) { box_set = b; }

private:
  GCoord thebox;
  bool box_set;
};


//! This class manages a shared Periodicbox
/** This is what most clients should use.  It maintains a shared
 * resource (the shared PeriodicBox) and forwards member function
 * calls to it.  The copy() member function creates a new
 * (i.e. dissociated PeriodicBox) and returns the associated
 * SharedPeriodicBox
 */

class SharedPeriodicBox {
public:
  SharedPeriodicBox() : pbox(new PeriodicBox) { }
  GCoord box(void) const { return(pbox->box()); }
  void box(const GCoord& c) { pbox->box(c); }
  bool isPeriodic(void) const { return(pbox->isPeriodic()); }


  SharedPeriodicBox copy(void) const {
    SharedPeriodicBox thecopy;

    if (isPeriodic())
      thecopy.box(box());
    
    return(thecopy);
  }

private:
  boost::shared_ptr<PeriodicBox> pbox;
};


#endif
