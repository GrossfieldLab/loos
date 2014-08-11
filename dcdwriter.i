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


%header %{
#include <dcdwriter.hpp>
%}

namespace loos {

  // Note, this is actually derived from boost::noncopyable
  class DCDWriter {
  public:
    explicit DCDWriter(const std::string& s, const bool append = false);
    DCDWriter(const std::string& s, const std::vector<AtomicGroup>& grps, const bool append = false);

    DCDWriter(const std::string& s, const std::vector<AtomicGroup>& grps, const std::string& comment, const bool append = false);
    DCDWriter(const std::string& s, const std::vector<AtomicGroup>& grps, const std::vector<std::string>& comments, const bool append = false);
    ~DCDWriter();
    void setHeader(const int na, const int ns, const greal ts, const bool bf);
    void setTitles(const std::vector<std::string>& titles);
    void setTitle(const std::string& s);
    void addTitle(const std::string& s);
    void writeFrame(const AtomicGroup& grp);
    void writeFrames(const std::vector<AtomicGroup>& grps);
    void writeHeader(void);
    int framesWritten(void) const { return(_current); }
  };

}

