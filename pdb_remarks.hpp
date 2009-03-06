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



#if !defined(PDB_REMARKS_HPP)
#define PDB_REMARKS_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

namespace loos {

  //! Class for handling PDB Remarks
  /**
   *This class just manages a vector of strings, but it will
   *truncate/pad the input strings to the appropriate length for a PDB
   *file and then output them with record numbers.
   */
  class Remarks {
  public:
    int numberOf(void) const;
    int size(void) const;
    //! Access the ith remark
    std::string get(const int i) const;
    //! Add a remark
    void add(const std::string s);
    //! Erase the ith remark
    void erase(const int i);

    //! Access the ith remark
    std::string& operator[](const int i);
    //! Access the ith remark
    const std::string& operator[](const int i) const;

    //! Returns a copy of the remarks vector
    std::vector<std::string> allRemarks(void) const { return(remarks); }

    //! Output the Remark(s) in PDB format
    friend std::ostream& operator<<(std::ostream& os, const Remarks& r);

  private:
    void rangeCheck(const unsigned int i) const;
    std::string sanitize(const std::string s) const;

  private:
    std::vector<std::string> remarks;
  };


}

#endif
