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
  string t(s);
  while (t.size() > 58) {
    remarks.push_back(t.substr(0, 58));
    t = t.substr(58);
  }

  string u = sanitize(t);
  remarks.push_back(sanitize(t));
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
  string t(s);
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
