/*
  UniqueStrings.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Find unique strings...
*/

#if !defined(UNIQUESTRINGS_HPP)
#define UNIQUESTRINGS_HPP


#include <string>
#include <ext/slist>

using namespace std;
using namespace __gnu_cxx;

//! Class for uniquifying strings...
/**  This class uses an slist and just does a linear search to see if
 *   we've already encountered the passed string before.  This probably
 *   should be updated to use a faster method in the future...
 */
class UniqueStrings {
public:

  //! Adds a string to the unique string list
  void add(const string& s) {
    if (find(s) < 0)
      uniques.push_front(s);
  }

  //! Number of unique strings found...
  int size(void) const { return(uniques.size()); }

  //! Returns the raw slist of strings...
  slist<string> strings(void) const { return(uniques); }

  //! Checks to see if we've encountered this string before...
  /** Returns an index (a unique int) representing this string.
   *  If the string is not found, returns -1.
   */
  int find(const string& s) {
    int j = 0;
    slist<string>::const_iterator i;
    for (i = uniques.begin(); i != uniques.end(); i++, j++)
      if (*i == s)
	return(j);
    return(-1); 
  }

private:
  slist<string> uniques;
};

#endif
