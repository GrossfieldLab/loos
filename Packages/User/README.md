		     * LOOS USER TOOL TEMPLATES *

This directory contains some tool-templates designed to make writing
new tools with LOOS easier.  Pick the tool-template that is closest to
what you need, then copy it into a new file and look for the regions
marked "***EDIT***".  Finally, add your tool into the list of "apps"
in the SConscript file.  You can now build your tool alongside LOOS.
In the top-level LOOS directory, type "scons user".  This will build
the LOOS library and everything in the user packages directory.

Templates:
   
* model_calc.cpp *

Calculates a statistic from a single structure.  Supports command-line
options using the LOOS OptionsFramework infrastructure.


* simple_model_calc.cpp *

Calculates a statistic from a single structure using a very simple
command-line format.


* simple_model_transform.cpp *

Performs a transformation of a single structure.  The new structure is
written to standard out as a PDB.  A very simple command-line format
is used.


* traj_cal.cpp *

Calculates a statistic from a single structure.  Command-line options
using the LOOS OptionsFramework infrastructure is supported.


* traj_transform.cpp *

Performs a transformation of a single structure.  The new structure is
written to standard out as a PDB.  Command-line options using the LOOS
OptionsFramework infrastructure is supported.
