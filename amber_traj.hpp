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

#if !defined(LOOS_AMBER_TRAJ_HPP)
#define LOOS_AMBER_TRAJ_HPP


#include <string>

#include <loos_defs.hpp>
#include <Coord.hpp>
#include <Trajectory.hpp>


namespace loos {

  //! Class for reading amber coordinate trajectories
  /*!
   * This class will read in the first frame of the trajectory upon
   * instantiation.  It will also scan the file to determine how many
   * frames there are.
   *
   * Since the Amber trajectory format does not store the # of atoms
   * present, this must be passed to the AmberTraj constructor.
   *
   * Note that the Amber timestep is (presumably) defined in the parmtop
   * file, not in the trajectory file.  So we return a null-value here...
   */

  class AmberTraj : public Trajectory {
  public:
    explicit AmberTraj(const std::string& s, const int na) : Trajectory(s),
                                                             _natoms(na), frame_offset(0),
                                                             frame_size(0), periodic(false) { init(); }

    explicit AmberTraj(const char* p, const int na) : Trajectory(p), _natoms(na),
                                                      frame_offset(0), frame_size(0),
                                                      periodic(false) { init(); }

    explicit AmberTraj(std::istream& is, const int na) : Trajectory(is), _natoms(na),
                                                     frame_offset(0), frame_size(0),
                                                     periodic(false) { init(); }

    std::string description() const { return("Amber trajectory"); }
    static pTraj create(const std::string& fname, const AtomicGroup& model) {
      return(pTraj(new AmberTraj(fname, model.size())));
    }


    virtual uint nframes(void) const { return(_nframes); }
    virtual uint natoms(void) const { return(_natoms); }
    virtual std::vector<GCoord> coords(void) { return(frame); }

    virtual bool hasPeriodicBox(void) const { return(periodic); }
    virtual GCoord periodicBox(void) const { return(box); }

    /*!
     * As stated above, Amber does not store the timestep in the
     * trajectory, but in the parmtop instead.  So we return a
     * null-value here...
     */
    virtual float timestep(void) const { return(0.0); }  // Dummy routine...

    virtual bool parseFrame(void);


  private:
    void init(void);
    virtual void rewindImpl(void) { ifs()->clear(); ifs()->seekg(frame_offset); }
    virtual void seekNextFrameImpl(void) { }
    virtual void seekFrameImpl(const uint);
    virtual void updateGroupCoordsImpl(AtomicGroup&);


  private:
    uint _natoms, _nframes;
    unsigned long frame_offset, frame_size;
    bool periodic;
    GCoord box;
    std::vector<GCoord> frame;

  };


}
#endif
