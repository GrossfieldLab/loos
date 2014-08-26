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

%shared_ptr(loos::TrajectoryWriter)

%header %{
#include <trajwriter.hpp>
#include <exception>
%}


namespace loos {



  class TrajectoryWriter;
  typedef boost::shared_ptr<TrajectoryWriter> pTrajectoryWriter;

  class TrajectoryWriter {

  public:

    struct WriteError : public std::exception {};

    TrajectoryWriter(const std::string& fname, const bool append = false);
    virtual ~TrajectoryWriter();
    virtual void setComments(const std::vector<std::string>& comments);
    virtual void setComments(const std::string& s);
    virtual void writeFrame(const AtomicGroup& model) =0;
    virtual void writeFrame(const AtomicGroup& model, const uint step, const double time);
    virtual bool hasFrameStep() const;
    virtual bool hasFrameTime() const;
    virtual bool hasComments() const;
    virtual uint framesWritten() const =0;
  };
}



