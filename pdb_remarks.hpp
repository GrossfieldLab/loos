/*
  pdb_remarks.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/

#if !defined(PDB_REMARKS_HPP)
#define PDB_REMARKS_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>


using namespace std;


#include <loos.hpp>


class Remarks {
public:
  int numberOf(void) const { return(remarks.size()); }   // Compat with PERL
  int size(void) const { return(remarks.size()); }
  string get(const int i) const;
  void add(const string s);
  void erase(const int i);

  string& operator[](const int i);
  const string& operator[](const int i) const;

  friend ostream& operator<<(ostream& os, const Remarks& r);

private:
  void rangeCheck(const unsigned int i) const;
  string sanitize(const string s) const;

private:
  vector<string> remarks;
};


#endif
