/*
  KernelActions.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

*/



#if !defined(KERNELACTIONS_HPP)
#define KERNELACTIONS_HPP

#include <sstream>
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
    string my_name;

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
    Action(const string s) : stack(0), atom(pAtom()), my_name(s) { }

    void setStack(ValueStack* ptr) { stack=ptr; }
    void setAtom(pAtom pa) { atom = pa; }

    virtual string name(void) const { return(my_name); }

    virtual void execute(void) = 0;
    virtual ~Action() { }

  };


  class pushString : public Action {
    Value val;
  public:
    pushString(const string str) : Action("pushString"), val(str) { }
    void execute(void) { stack->push(val); }
    string name(void) const {
      stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }

  };

  class pushInt : public Action {
    Value val;
  public:
    pushInt(const int i) : Action("pushInt"), val(i) { }
    void execute(void) { stack->push(val); }
    string name(void) const {
      stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }
  };

  class pushFloat : public Action {
    Value val;
  public:
    pushFloat(const float f) : Action("pushFloat"), val(f) { }
    void execute(void) { stack->push(val); }
    string name(void) const {
      stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }
  };

  class drop : public Action {
  public:
    drop() : Action("drop") { }
    void execute(void) { stack->drop(); }
  };

  class dup : public Action {
  public:
    dup() : Action("dup") { }
    void execute(void) { stack->dup(); }
  };

  class equals : public Action {
  public:
    equals() : Action("==") { }
    void execute(void) {
      Value v(binComp() == 0);
      stack->push(v);
    }
  };

  class lessThan : public Action {
  public:
    lessThan() : Action("<") { }
    void execute(void) {
      Value v(binComp() < 0);
      stack->push(v);
    }
  };

  class lessThanEquals : public Action {
  public:
    lessThanEquals() : Action("<=") { }
    void execute(void) {
      Value v(binComp() <= 0);
      stack->push(v);
    }
  };

  class greaterThan : public Action {
  public:
    greaterThan() : Action(">") { }
    void execute(void) {
      Value v(binComp() > 0);
      stack->push(v);
    }
  };

  class greaterThanEquals : public Action {
  public:
    greaterThanEquals() : Action(">=") { }
    void execute(void) {
      Value v(binComp() >= 0);
      stack->push(v);
    }
  };

  class matchRegex : public Action {
    boost::regex regexp;
  public:
    matchRegex(const string s) : Action("matchRegex"), regexp(s, boost::regex::perl|boost::regex::icase) { }
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
    matchStringAsRegex() : Action("matchStringAsRegex") { }
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
    pushAtomName() : Action("pushAtomName") { }
    void execute() {
      hasAtom();
      Value v(atom->name());
      stack->push(v);
    }
  };


  class pushAtomId : public Action {
  public:
    pushAtomId() : Action("pushAtomId") { }
    void execute() {
      hasAtom();
      Value v(atom->id());
      stack->push(v);
    }
  };

  class pushAtomResname : public Action {
  public:
    pushAtomResname() : Action("pushAtomResname") { }
    void execute() {
      hasAtom();
      Value v(atom->resname());
      stack->push(v);
    }
  };

  class pushAtomResid : public Action {
  public:
    pushAtomResid() : Action("pushAtomResid") { }
    void execute() {
      hasAtom();
      Value v(atom->resid());
      stack->push(v);
    }
  };

  class pushAtomSegid : public Action {
  public:
    pushAtomSegid() : Action("pushAtomSegid") { }
    void execute() {
      hasAtom();
      Value v(atom->segid());
      stack->push(v);
    }
  };

  class logicalAnd : public Action {
  public:
    logicalAnd() : Action("&&") { }
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
    logicalOr() : Action("||") { }
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
    logicalNot() : Action("!") { }
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

