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

#include <vector>

#include "KernelValue.hpp"
#include "KernelStack.hpp"
#include "KernelActions.hpp"

using namespace std;

namespace loos {

  class Kernel {
    vector<Action*> actions;
    ValueStack val_stack;

  public:
    
    ~Kernel() {
      vector<Action*>::iterator i;
      for (i = actions.begin(); i != actions.end(); i++)
	delete (*i);
    }

    void push(Action *act) {
      act->setStack(&val_stack);
      actions.push_back(act); 
    }

    void pop(void) { actions.pop_back(); }

    void execute(void) {
      vector<Action*>::iterator i;
      for (i=actions.begin(); i != actions.end(); i++)
	(*i)->execute();
    }

    
    void clearActions(void) { actions.clear(); }

    ValueStack& stack(void) { return(val_stack); }
    

  };

};


#endif
