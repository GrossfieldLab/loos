/*
  Kernel.hpp
  (c) 2008 Tod D. Romo


  Grossfield Lab
  Department of Biochemistry and Biophysics
  University of Rochester Medical School

  Kernel for handling atom selections...
*/



#if !defined(KERNEL_HPP)
#define KERNEL_HPP

#include <deque>

#include "KernelValue.hpp"
#include "KernelStack.hpp"
#include "KernelActions.hpp"

using namespace std;

namespace loos {

  class Kernel {
    deque<Action*> actions;
    ValueStack val_stack;

  public:
    
    ~Kernel() {
      deque<Action*>::iterator i;
      for (i = actions.begin(); i != actions.end(); i++)
	delete (*i);
    }

    void add(Action *act) {
      act->setStack(&val_stack);
      actions.push_back(act); 
    }

    void execute(void) {
      while (actions.size() > 0) {
	//	cout << val_stack << "-------\n";
	Action* act = actions.front();
	actions.pop_front();
	act->execute();
      }
    }

    ValueStack& stack(void) { return(val_stack); }

  };

};


#endif
