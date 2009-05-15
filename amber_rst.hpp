/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2009, Tod D. Romo, Alan Grossfield
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

#if !defined(AMBER_RST_HPP)
#define AMBER_RST_HPP


#include <string>

#include <loos_defs.hpp>
#include <Trajectory.hpp>


namespace loos {

  //! Class for reading amber restart files as a single-frame trajectory

  class AmberRst : public Trajectory {
  public:
    explicit AmberRst(const std::string& s, const int na) : Trajectory(s), _natoms(na),
                                                            frame_offset(0), frame_size(0),
                                                            periodic(false), seek_flag(false) {
      if (!parseFrame())
        throw(std::runtime_error("Could not open the Amber RST file"));
    }



    explicit AmberRst(const char* p, const int na) : Trajectory(p), _natoms(na),
                                                     frame_offset(0), frame_size(0),
                                                     periodic(false), seek_flag(false) {
      if (!parseFrame())
        throw(std::runtime_error("Could not open the Amber RST file"));
    }


    virtual void rewind(void) { seek_flag = false; cached_first = true; }
    virtual uint nframes(void) const { return(1); }
    virtual uint natoms(void) const { return(_natoms); }
    virtual std::vector<GCoord> coords(void) { return(frame); }
    virtual void updateGroupCoords(AtomicGroup&);

    virtual void seekNextFrame(void);
    virtual void seekFrame(const uint);
    virtual bool parseFrame(void);

    virtual bool hasPeriodicBox(void) const { return(periodic); }
    virtual GCoord periodicBox(void) const { return(box); }

    virtual greal timestep(void) const { return(0.0); }  // Dummy routine...
    greal currentTime(void) const { return(current_time); }
    

  private:
    uint _natoms, _nframes;
    unsigned long frame_offset, frame_size;
    greal current_time;
    bool periodic;
    GCoord box;
    std::vector<GCoord> frame;

    bool seek_flag;
  };


}
#endif
