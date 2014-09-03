/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2014, Tod D. Romo, Alan Grossfield
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

%shared_ptr(loos::XTCWriter)


%header %{
#include <xtcwriter.hpp>
%}


namespace loos {
  class XTCWriter : public TrajectoryWriter {
  public:


    XTCWriter(const std::string fname, const bool append = false);
    XTCWriter(const std::string fname, const double dt, const uint steps_per_frame, const bool append = false);
    XTCWriter(const std::string fname, const double dt, const uint steps_per_frame, const float precision, const bool append = false);
    ~XTCWriter();

    double timePerStep() const;
    void timePerStep(const double dt);
    uint stepsPerFrame() const;
    void stepsPerFrame(const uint s);
    uint currentStep() const;
    void currentStep(const uint s);
    void writeFrame(const AtomicGroup& model);
    void writeFrame(const AtomicGroup& model, const uint step, const double time);
    uint framesWritten() const;
  };
}
