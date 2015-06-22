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





#if !defined(LOOS_KERNELACTIONS_HPP)
#define LOOS_KERNELACTIONS_HPP


#include <iostream>
#include <string>
#include <stdexcept>
#include <string.h>


#include <loos_defs.hpp>
#include <boost/regex.hpp>

#include <exceptions.hpp>

#include "KernelValue.hpp"
#include "KernelStack.hpp"



namespace loos {

  class BackboneSelector;

  //! Loos esoterica.
  /** You probably don't want to look in here unless you want to
   *  program with the virtual machine for atom selections, but I'd
   *  advise against that...
   */
  namespace internal {

    //! Base class for all commands...
    /** All subclasses must implement the execute() method, which will
     *  operate on the data stack pointer.
     *
     *  Subclasses may also override the name() method if they want to
     *  augment the command-name string (i.e. to show additional internal
     *  data)
     */

  
    class Action {
    protected:
      //! Pointer to the data stack
      ValueStack *stack;
      //! Pointer to the atom we'll be working on...
      pAtom atom;
      //! Record of command-name (for printing)
      std::string my_name;

      // Some utility functions...

      //! Compare the top two items on the stack...
      int binComp(void);

      bool negativeOperand(void);

      void binaryFalseResult(void);

      //! Check to make sure an atom has been set...
      void requireAtom(void);

    public:
      Action(const std::string s) : stack(0), atom(pAtom()), my_name(s) { }

      void setStack(ValueStack*);
      void setAtom(pAtom&);

      virtual std::string name(void) const;

      virtual void execute(void) =0;
      virtual ~Action() { }

    };



    //! Push a std::string onto the data stack
    class pushString : public Action {
      Value val;
    public:
      explicit pushString(const std::string str) : Action("pushString"), val(str) { }
      void execute(void);
      std::string name(void) const;
    };

    //! Push an integer onto the data stack
    class pushInt : public Action {
      Value val;
    public:
      explicit pushInt(const long i) : Action("pushInt"), val(i) { }
      void execute(void);
      std::string name(void) const;
    };

    //! Push a float onto the data stack
    class pushFloat : public Action {
      Value val;
    public:
      explicit pushFloat(const float f) : Action("pushFloat"), val(f) { }
      void execute(void);
      std::string name(void) const;
    };


    //! Drop the top item from the stack
    class drop : public Action {
    public:
      drop() : Action("drop") { }
      void execute(void);
    };

    //! Duplicate the top item on the stack
    class dup : public Action {
    public:
      dup() : Action("dup") { }
      void execute(void);
    };


    //! Test for equality: ARG1 ARG ==
    class equals : public Action {
    public:
      equals() : Action("==") { }
      void execute(void);
    };

    //! Relation operators...:  ARG1 ARG2 <
    class lessThan : public Action {
    public:
      lessThan() : Action("<") { }
      void execute(void);
    };

    //! ARG1 ARG2 <=
    class lessThanEquals : public Action {
    public:
      lessThanEquals() : Action("<=") { }
      void execute(void);
    };

    //! ARG1 ARG2 >
    class greaterThan : public Action {
    public:
      greaterThan() : Action(">") { }
      void execute(void);
    };

    //! ARG1 ARG2 >=
    class greaterThanEquals : public Action {
    public:
      greaterThanEquals() : Action(">=") { }
      void execute(void);
    };

    //! Regular expression matching: ARG1 regexp(S)
    /** Compiles the passed string into a regex pattern at instantiation,
     * then at execution matches the top stack entry against the
     * pattern...
     */

    class matchRegex : public Action {
      boost::regex regexp;
    public:
      explicit matchRegex(const std::string s) : Action("matchRegex"), regexp(s, boost::regex::perl|boost::regex::icase), pattern(s) { }
      void execute(void);
      std::string name(void) const;
    
    private:
      std::string pattern;
    };
  

    //! Regular expression matching: ARG1 regexp(ARG2)
    /** Takes the top item on the stack and compiles this into a regular
     * expression, then matches it against the next item on the stack.
     * This is likely to be pretty inefficient, so it's better to use
     * matchRegex instead if you can.
     */
    class matchStringAsRegex : public Action {
    public:
      matchStringAsRegex() : Action("matchStringAsRegex") { }
      void execute(void);
    };
  
  
    //! Extracts a number for a string on the stack using a regular expression: ARG1 regexp(S)
    /** Compiles the passed string into a regex pattern at
     * instantiation.  At execution, examines each matched capture for
     * the first one that converts to an integer and pushes that value
     * onto the data stack.  If no match is found (or no numeric
     * conversion works), then "-1" is pushed onto the stack.
     */
    class extractNumber : public Action {
    public:
      explicit extractNumber(const std::string s) : Action("extractNumber"),
                                                    regexp(s, boost::regex::perl|boost::regex::icase),
                                                    pattern(s) { }

      void execute(void);
      std::string name(void) const;

    private:
      boost::regex regexp;
      std::string pattern;
    };




    //! Push atom name onto the stack
    class pushAtomName : public Action {
    public:
      pushAtomName() : Action("pushAtomName") { }
      void execute(void);
    };

    //! Push atom id onto the stack
    class pushAtomId : public Action {
    public:
      pushAtomId() : Action("pushAtomId") { }
      void execute(void);
    };

    //! Push atom'ss residue name onto the stack
    class pushAtomResname : public Action {
    public:
      pushAtomResname() : Action("pushAtomResname") { }
      void execute(void);
    };

    //! Push atom's residue id onto the stack
    class pushAtomResid : public Action {
    public:
      pushAtomResid() : Action("pushAtomResid") { }
      void execute(void);
    };

    //! Push atom's segid onto the stack
    class pushAtomSegid : public Action {
    public:
      pushAtomSegid() : Action("pushAtomSegid") { }
      void execute(void);
    };

    //! Push atom's chain ID onto the stack
    class pushAtomChainId : public Action {
    public:
      pushAtomChainId() : Action("pushAtomChainId") { }
      void execute(void);
    };



    // Logical operations...  Assumes stack args are ints...

    //! ARG1 ARG2 &&
    class logicalAnd : public Action{
    public:
      logicalAnd() : Action("&&") { }
      void execute(void);
    };

    //! ARG1 ARG2 ||
    class logicalOr : public Action {
    public:
      logicalOr() : Action("||") { }
      void execute(void);
    };


    //! ARG1 !
    class logicalNot : public Action {
    public:
      logicalNot() : Action("!") { }
      void execute(void);
    };

    //! Always returns true...
    class logicalTrue : public Action {
    public:
      logicalTrue() : Action("TRUE") { }
      void execute(void);
    };


    //! Shortcut for checking for hydrogens...
    class Hydrogen : public Action {
    public:
      Hydrogen() : Action("Hydrogen") { }
      void execute(void);
    };

    //! Shortcut for checking for backbone atoms...
    class Backbone : public Action {
      static BackboneSelector bbsel;    // Make this class-level since we only need one to wrap
    public:
      Backbone() : Action("Backbone") { }
      void execute(void);
    };
  

  }

}



#endif

