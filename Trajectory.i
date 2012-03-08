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

%shared_ptr(loos::Trajectory)

%header %{
#include <Trajectory.hpp>
#include <StreamWrapper.hpp>
%}



namespace loos {

  // Note, this is actually derived from boost::noncopyable
  class Trajectory;
  typedef boost::shared_ptr<Trajectory> pTraj;
  class Trajectory {
  public:
    Trajectory();
    Trajectory(const std::string& s);
    Trajectory(const char* s);
    virtual ~Trajectory();

    virtual uint natoms(void) const =0;
    virtual float timestep(void) const =0;
    virtual uint nframes(void) const =0;
    void rewind(void);
    virtual bool hasPeriodicBox(void) const =0;
    virtual GCoord periodicBox(void) const =0;
    virtual std::vector<GCoord> coords(void) =0;
    virtual void updateGroupCoords(AtomicGroup& g) =0;
    void seekNextFrame(void);
    void seekFrame(const uint i);
    virtual bool parseFrame(void) =0;
    bool readFrame(void);
    bool readFrame(const int i);
  };

}
