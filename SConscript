#!/usr/bin/env python3
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2016, Tod D. Romo
#  Department of Biochemistry and Biophysics
#  School of Medicine & Dentistry, University of Rochester
#
#  This package (LOOS) is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation under version 3 of the License.
#
#  This package is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


import sys
import os

Import('env')

# Handle installation...
PREFIX = env['PREFIX']

env.Install(os.path.join(PREFIX, 'lib'), 'loos')


# Setup environment script(s)

script_sh = env.Scripts('setup.sh', 'setup.sh-pre')
script_csh = env.Scripts('setup.csh', 'setup.csh-pre')
scripts = [ script_sh, script_csh ]

script_sh_inst = env.Scripts(os.path.join(PREFIX, 'setup.sh'), 'setup.sh-pre')
script_csh_inst = env.Scripts(os.path.join(PREFIX, 'setup.csh'), 'setup.csh-pre')
scripts_inst = [ script_sh_inst, script_csh_inst ]


Return('scripts')
