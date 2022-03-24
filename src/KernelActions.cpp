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

#include <KernelActions.hpp>
#include <Atom.hpp>
#include <sstream>
#include <Selectors.hpp>


namespace loos {

  namespace internal {


    int Action::binComp(void) {
      Value v1 = stack->pop();
      Value v2 = stack->pop();
      return(compare(v2, v1));
    }

    bool Action::negativeOperand(void) {
      Value v1 = stack->peek(-1);
      Value v2 = stack->peek(-2);

      if ( (v1.type == Value::INT && v1.itg < 0) ||
           (v2.type == Value::INT && v2.itg < 0) )
        return(true);

      return(false);
    }


    void Action::binaryFalseResult(void) {
      Value v(0);

      stack->drop();
      stack->drop();
      stack->push(v);
    }


    void Action::requireAtom(void) {
      if (atom == 0)
        throw(LOOSError("No atom set"));
    }


    void Action::setStack(ValueStack* ptr) { stack=ptr; }
    void Action::setAtom(pAtom& pa) { atom = pa; }

    std::string Action::name(void) const { return(my_name); }

    //-------------------------------------------------------------


    void pushString::execute(void) { stack->push(val); }
    std::string pushString::name(void) const {
      std::stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }

    void pushInt::execute(void) { stack->push(val); }
    std::string pushInt::name(void) const {
      std::stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }

    void pushFloat::execute(void) { stack->push(val); }
    std::string pushFloat::name(void) const {
      std::stringstream s;
      s << my_name << "(" << val << ")";
      return(s.str());
    }


    void drop::execute(void) { stack->drop(); }

    void dup::execute(void) { stack->dup(); }

    void equals::execute(void) {
      Value v(binComp() == 0);
      stack->push(v);
    }

    void lessThan::execute(void) {
      if (negativeOperand())
        binaryFalseResult();
      else {
        Value v(binComp() < 0);
        stack->push(v);
      }
    }

    void lessThanEquals::execute(void) {
      if (negativeOperand())
        binaryFalseResult();
      else {
        Value v(binComp() <= 0);
        stack->push(v);
      }
    }


    void greaterThan::execute(void) {
      Value v(binComp() > 0);
      stack->push(v);
    }

    void greaterThanEquals::execute(void) {
      Value v(binComp() >= 0);
      stack->push(v);
    }

    void matchRegex::execute(void) {
      Value v = stack->pop();
      Value r(0);
      if (boost::regex_search(v.getString(), regexp))
        r.setInt(1);

      stack->push(r);
    }

    std::string matchRegex::name(void) const {
      return(my_name + "(" + pattern + ")");
    }



    void matchStringAsRegex::execute(void) {
      Value v = stack->pop();
      boost::regex re(v.getString(), boost::regex::perl|boost::regex::icase);
      Value u = stack->pop();
      Value r(0);

      if (boost::regex_search(u.getString(), re))
        r.setInt(1);

      stack->push(r);
    }

    void extractNumber::execute(void) {
      Value v = stack->pop();
      Value r(-1);
      boost::smatch what;

      if (boost::regex_search(v.getString(), what, regexp)) {
        unsigned i;
        int val;
        for (i=0; i<what.size(); i++) {
          if ((std::stringstream(what[i]) >> val)) {
            r.setInt(val);
            break;
          }
        }
      }

      stack->push(r);
    }

    std::string extractNumber::name(void) const {
      return(my_name + "(" + pattern + ")");
    }


    void pushAtomName::execute(void) {
      requireAtom();
      Value v(atom->name());
      stack->push(v);
    }

    void pushAtomId::execute(void) {
      requireAtom();
      Value v(atom->id());
      stack->push(v);
    }

    // Beware of overflows here!!!
    void pushAtomIndex::execute(void) {
      requireAtom();
      Value v(static_cast<long>(atom->index()));
      stack->push(v);
    }

    void pushAtomResname::execute(void) {
      requireAtom();
      Value v(atom->resname());
      stack->push(v);
    }

    void pushAtomResid::execute(void) {
      requireAtom();
      Value v(atom->resid());
      stack->push(v);
    }


    void pushAtomSegid::execute(void) {
      requireAtom();
      Value v(atom->segid());
      stack->push(v);
    }

    void pushAtomChainId::execute(void) {
      requireAtom();
      Value v(atom->chainId());
      stack->push(v);
    }


    void logicalAnd::execute(void) {
      Value v2 = stack->pop();
      Value v1 = stack->pop();

      if (!(v1.type == Value::INT && v2.type == Value::INT))
        throw(LOOSError("Invalid operands to logicalAnd"));

      Value u(v1.itg && v2.itg);
      stack->push(u);
    }


    void logicalOr::execute(void) {
      Value v1 = stack->pop();
      Value v2 = stack->pop();

      if (!(v1.type == Value::INT && v2.type == Value::INT))
        throw(LOOSError("Invalid operands to logicalOr"));

      Value u(v1.itg || v2.itg);
      stack->push(u);
    }


    void logicalNot::execute(void) {
      Value v1 = stack->pop();

      if (v1.type != Value::INT)
        throw(LOOSError("Invalid operand to logicalNot"));

      Value u(!v1.itg);
      stack->push(u);
    }

    void logicalTrue::execute(void) {
      Value v((int)1);
      stack->push(v);
    }


    void Hydrogen::execute(void) {
      requireAtom();

      bool masscheck = true;
      if (atom->checkProperty(Atom::massbit))
        // Note: Checking for mass < 4.1 allows use of hydrogen selector
        //       even when the system has hydrogen mass repartitioning.
        //       False positive if someone has He in the system.
        masscheck = (atom->mass() < 4.1);

      std::string n = atom->name();
      Value v;
      v.setInt( (n[0] == 'H' && masscheck) );
      stack->push(v);
    }



    // Provide storage for class-level selector
    BackboneSelector Backbone::bbsel;


    void Backbone::execute(void) {
      requireAtom();
      Value v(bbsel(atom));

      stack->push(v);
    }


  }

}
