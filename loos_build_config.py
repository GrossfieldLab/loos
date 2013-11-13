#!/usr/bin/env python
#  This file is part of LOOS.
#
#  LOOS (Lightweight Object-Oriented Structure library)
#  Copyright (c) 2013, Tod D. Romo
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


loos_version = '2.1.0'
default_lib_path = '/usr/lib64'
min_boost_version = 36    # 1.xx

host_type = 'unknown'
linux_type = 'unknown'
suffix = 'so'

# List of packages to build and install
# The keywords become build targets for SCons
package_list = { 'ENM': 'ElasticNetworks',
                 'HBonds' : 'HydrogenBonds',
                 'Conv' : 'Convergence',
                 'Density' : 'DensityTools',
                 'User': 'User',
                 'Python': 'PyLOOS' }



