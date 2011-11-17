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



#if !defined(LOOS_UNIQUESTRINGS_HPP)
#define LOOS_UNIQUESTRINGS_HPP


#include <string>
#include <ext/slist>


namespace loos {

  //! Class for uniquifying strings...
  /**  This class uses an slist and just does a linear search to see if
   *   we've already encountered the passed string before.  This probably
   *   should be updated to use a faster method in the future...
   */
  class UniqueStrings {
  public:

    //! Adds a string to the unique string list
    void add(const std::string& s) {
      if (find(s) < 0)
        uniques.push_front(s);
    }

    //! Number of unique strings found...
    int size(void) const { return(uniques.size()); }

    //! Returns the raw slist of strings...
    __gnu_cxx::slist<std::string> strings(void) const { return(uniques); }

    //! Checks to see if we've encountered this string before...
    /** Returns an index (a unique int) representing this string.
     *  If the string is not found, returns -1.
     */
    int find(const std::string& s) {
      int j = 0;
      __gnu_cxx::slist<std::string>::const_iterator i;
      for (i = uniques.begin(); i != uniques.end(); i++, j++)
        if (*i == s)
          return(j);
      return(-1); 
    }

  private:
    __gnu_cxx::slist<std::string> uniques;
  };

}
#endif
