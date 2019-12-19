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


loos_version = '3.1'
min_boost_version = '1_36'
required_boost_libraries = ('regex', 'system', 'program_options', 'thread', 'filesystem')

min_swig_version = '2.0'

host_type = 'unknown'
linux_type = 'unknown'
suffix = 'so'

versions = {}

# List of packages to build and install
# The keywords become build targets for SCons
package_list = { 'ENM': 'ElasticNetworks',
                 'HBonds' : 'HydrogenBonds',
                 'Conv' : 'Convergence',
                 'Clustering' : 'Clustering',
                 'Density' : 'DensityTools',
                 'OMG' : 'OptimalMembraneGenerator',
                 'Voronoi' : 'Voronoi',
                 'User': 'User',
                 'Python': 'PyLOOS' }


user_libdirs = {}
user_boost_flag = 0     # 1 = user has overridden boost location
