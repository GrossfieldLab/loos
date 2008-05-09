/*
  pdb_remarks.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/


#include <ios>
#include <sstream>
#include <iostream>
#include <iomanip>


#include <pdb_remarks.hpp>



void Remarks::rangeCheck(const unsigned int i) const {
  if (i >= remarks.size())
    throw(range_error("Bad indices into remarks"));
}


string Remarks::get(const int i) const {
  rangeCheck(i);
  return(remarks[i]);
}


void Remarks::add(const string s) {
  string t = sanitize(s);
  remarks.push_back(t);
}


void Remarks::erase(const int i) {
  rangeCheck(i);
  remarks.erase(remarks.begin() + i);
}


string& Remarks::operator[](const int i) {
  rangeCheck(i);
  return(remarks[i]);
}


const string& Remarks::operator[](const int i) const {
  rangeCheck(i);
  return(remarks[i]);
}



string Remarks::sanitize(const string s) const {
  string t;
  int n = s.size();

  if (n > 58)
    t = s.substr(0, 58);
  else if (n < 58)
    t = s + string(58 - n, ' ');

  return(t);
}



ostream& operator<<(ostream& os, const Remarks& r) {
  vector<string>::const_iterator i;
  int n = 1;

  for (i = r.remarks.begin(); i != r.remarks.end(); i++) {
    string s("REMARK                                                                \n");
    ostringstream ss;

    ss << setw(3) << setfill('0') << n++;
    s.replace(7,3, ss.str());
    s.replace(11,58, *i);
    os << s;
  }

  return(os);
}
