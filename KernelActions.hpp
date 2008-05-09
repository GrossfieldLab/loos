/*
  KernelActions.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/



#if !defined(KERNELACTIONS_HPP)
#define KERNELACTIONS_HPP

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <string.h>

#include <boost/regex.hpp>

#include "Atom.hpp"

#include "KernelValue.hpp"
#include "KernelStack.hpp"


using namespace std;

namespace loos {
  class Action {
  protected:
    ValueStack *stack;
    pAtom atom;

    int binComp(void) {
      Value v1 = stack->pop();
      Value v2 = stack->pop();
      return(compare(v1, v2));
    }

    void hasAtom(void) {
      if (atom == 0)
	throw(runtime_error("No atom set"));
    }

  public:
    Action() : stack(0) { }
    Action(ValueStack *stk) : stack(stk) { }

    void setStack(ValueStack* ptr) { stack=ptr; }
    void setAtom(pAtom pa) { atom = pa; }
    

    virtual void execute(void) = 0;
    virtual ~Action() { }
  };


  class pushString : public Action {
    Value val;
  public:
    pushString(const string str) : val(str) { }
    void execute(void) { stack->push(val); }
  };

  class pushInt : public Action {
    Value val;
  public:
    pushInt(const int i) : val(i) { }
    void execute(void) { stack->push(val); }
  };

  class pushFloat : public Action {
    Value val;
  public:
    pushFloat(const float f) : val(f) { }
    void execute(void) { stack->push(val); }
  };

  class drop : public Action {
  public:
    void execute(void) { stack->drop(); }
  };

  class dup : public Action {
  public:
    void execute(void) { stack->dup(); }
  };

  class equals : public Action {
  public:
    void execute(void) {
      Value v(binComp() == 0);
      stack->push(v);
    }
  };

  class lessThan : public Action {
  public:
    void execute(void) {
      Value v(binComp() < 0);
      stack->push(v);
    }
  };

  class lessThanEquals : public Action {
  public:
    void execute(void) {
      Value v(binComp() <= 0);
      stack->push(v);
    }
  };

  class greaterThan : public Action {
  public:
    void execute(void) {
      Value v(binComp() > 0);
      stack->push(v);
    }
  };

  class greaterThanEquals : public Action {
  public:
    void execute(void) {
      Value v(binComp() >= 0);
      stack->push(v);
    }
  };

  class matchRegex : public Action {
    boost::regex regexp;
  public:
    matchRegex(const string s) : regexp(s, boost::regex::perl|boost::regex::icase) { }
    void execute(void) { 
      Value v = stack->pop();
      Value r(0);
      if (boost::regex_search(v.getString(), regexp))
	r.setInt(1);

      stack->push(r);
    }
  };

  class matchStringAsRegex : public Action {
  public:
    void execute() {
      Value v = stack->pop();
      boost::regex re(v.getString(), boost::regex::perl|boost::regex::icase);
      Value u = stack->pop();
      Value r(0);
      
      if (boost::regex_search(u.getString(), re))
	r.setInt(1);

      stack->push(r);
    }
  };


  class pushAtomName : public Action {
  public:
    void execute() {
      hasAtom();
      Value v(atom->name());
      stack->push(v);
    }
  };


  class pushAtomId : public Action {
  public:
    void execute() {
      hasAtom();
      Value v(atom->id());
      stack->push(v);
    }
  };

  class pushAtomResname : public Action {
  public:
    void execute() {
      hasAtom();
      Value v(atom->resname());
      stack->push(v);
    }
  };

  class pushAtomResid : public Action {
  public:
    void execute() {
      hasAtom();
      Value v(atom->resid());
      stack->push(v);
    }
  };

  class pushAtomSegid : public Action {
  public:
    void execute() {
      hasAtom();
      Value v(atom->segid());
      stack->push(v);
    }
  };

  class logicalAnd : public Action {
  public:
    void execute() {
      Value v1 = stack->pop();
      Value v2 = stack->pop();

      if (!(v1.type == Value::INT && v2.type == Value::INT))
	throw(runtime_error("Invalid operands to logicalAnd"));

      Value u(v1.itg && v2.itg);
      stack->push(u);
    }
  };


  class logicalOr : public Action {
  public:
    void execute() {
      Value v1 = stack->pop();
      Value v2 = stack->pop();

      if (!(v1.type == Value::INT && v2.type == Value::INT))
	throw(runtime_error("Invalid operands to logicalOr"));

      Value u(v1.itg || v2.itg);
      stack->push(u);
    }
  };


  class logicalNot : public Action {
  public:
    void execute() {
      Value v1 = stack->pop();

      if (v1.type != Value::INT)
	throw(runtime_error("Invalid operand to logicalNot"));

      Value u(!v1.itg);
      stack->push(u);
    }
  };



};



#endif

