/*
  KernelStack.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/



#if !defined(KERNELSTACK_HPP)
#define KERNELSTACK_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <string.h>

#include "KernelValue.hpp"


using namespace std;

namespace loos {


  class ValueStack {
    vector<Value> values;

    void notEmpty(void) const {
      if (values.size() == 0)
	throw(logic_error("Operation requested on an empty stack."));
    }

  public:
    void push(const Value& val) { values.push_back(val); };

    Value pop(void) {
      notEmpty();
      Value val = values.back();
      values.pop_back();
      return(val); 
    }

    void dup(void) {
      notEmpty();
      Value val = values.back();
      push(val);
    }

    void drop(void) {
      notEmpty();
      values.pop_back();
    }

    Value peek(int i) {
      if (i < 0)
	i += values.size();
      if ((unsigned int)i >= values.size())
	throw(logic_error("Peeking beyond the stack!"));

      return(values[i]);
    }

    unsigned int size(void) const { return(values.size()); }

    void clear(void) { values.clear(); }

    friend ostream& operator<<(ostream& os, const ValueStack& s) {
      os << "<STACK>\n";
      vector<Value>::const_iterator i;

      for (i=s.values.begin(); i != s.values.end(); i++)
	os << "  " << *i << endl;

      return(os);
    }
  };

};


#endif
