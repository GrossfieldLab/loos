
import sys
import os
import glob
import platform
import re
from subprocess import *
from time import strftime
import shutil
import distutils.sysconfig
import distutils.spawn
from string import Template

import SCons

import loos_globals


# Attempt to canonicalize system name, type and other related info...
def canonicalizeSystem():
    linux_type = 'nonlinux'
    host_type = platform.system()
    suffix = 'so'
# Detect CYGWIN & canonicalize linux type, setting defaults...
    if (re.search("(?i)cygwin", host_type)):
        host_type = 'Cygwin'
    elif (host_type == 'Linux'):
        # Determine linux variant...
        linux_type = platform.platform()
        
        if (re.search("(?i)ubuntu", linux_type)):
            linux_type = 'ubuntu'
            loos_globals.default_lib_path = '/usr/lib'
        elif (re.search("(?i)suse", linux_type)):
            linux_type = 'suse'
        elif (re.search("(?i)debian", linux_type)):
            linux_type = 'debian'
    # MacOS is special (of course...)
    elif (host_type == 'Darwin'):
        loos_globals.default_lib_path = '/usr/bin'
        suffix = 'dylib'
    return(host_type, linux_type, suffix)




### Create a revision file for linking against.
def setupRevision(env):

    # Divine the current revision...
    revision = loos_globals.loos_version + " " + strftime("%y%m%d")

    # Now, write this out to a cpp file that can be linked in...this avoids having
    # to recompile everything when building on a new date.  We also rely on SCons
    # using the MD5 checksum to detect changes in the file (even though it's always
    # rewritten)
    revfile = open('revision.cpp', 'w')
    revfile.write('#include <string>\n')
    revfile.write('std::string revision_label = "')
    revfile.write(revision)
    revfile.write('";\n')
    revfile.close()


### Let environment variables override or modify some build paramaters...
def environOverride(conf):
    # Allow overrides from environment...
    if 'CXX' in os.environ:
        conf.env.Replace(CXX = os.environ['CXX'])
        print '*** Using compiler ' + os.environ['CXX']
    
    if 'CCFLAGS' in os.environ:
        conf.env.Append(CCFLAGS = os.environ['CCFLAGS'])
        print '*** Appending custom build flags: ' + os.environ['CCFLAGS']
        
    if 'LDFLAGS' in os.environ:
        conf.env.Append(LINKFLAGS = os.environ['LDFLAGS'])
        print '*** Appending custom link flag: ' + os.environ['LDFLAGS']



### Builder for setup scripts

# This copies the environment setup script while changing the directory
# that's used for setting up PATH and [DY]LD_LIBRARY_PATH.  If LOOS
# is being built in a directory, the env script will be setup to use
# the built-in-place distribution.  If LOOS is being installed, then
# it will use the installation directory instead.

def script_builder_python(target, source, env):
   if 'LOOS_PATH' in env:
       dir_path = env['LOOS_PATH']
   else:
       dir_path = env.Dir('.').abspath

   libpaths = env['LIBPATH']
   libpaths.pop(0)
   libpaths.insert(0, env.Dir('.').abspath)

   cpppaths = env['CPPPATH']
   cpppaths.pop(0)
   cpppaths.insert(0, env.Dir('.').abspath)

   if not 'install' in SCons.Script.COMMAND_LINE_TARGETS:
       toolpath = '$LOOS/Tools:' + ':'.join(['$LOOS/Packages/' + s for s in [loos_globals.package_list[i] for i in loos_globals.package_list]])

   else:
       toolpath = dir_path + '/' + 'bin'
       

   file = open(str(source[0]), 'r')
   script = file.read()
   script_template = Template(script)
   script = script_template.substitute(loos_path = dir_path,
                                       tool_path = toolpath,
                                       libpath = ':'.join(libpaths),
                                       cpppath = ':'.join(cpppaths),
                                       linkflags = env['LINKFLAGS'],
                                       libs = ':'.join(env['LIBS']),
                                       ccflags = env['CCFLAGS'])

   outfile = open(str(target[0]), 'w')
   outfile.write(script)

   return None



# Verify that we have swig and it's v2.0+
# Returns the path to swig
def CheckForSwig(conf):
    conf.Message('Checking for Swig v2.0+ ...')
    swig_location = distutils.spawn.find_executable('swig', conf.env['ENV']['PATH'])
    if swig_location == None:
        conf.Result('no')
    else:
        swig_check = Popen([swig_location, "-version"], stdout=PIPE).communicate()[0]
        swig_version = swig_check.split('\n')[1].split(' ')[2]
        swig_major = swig_version.split('.')[0]
        if int(swig_major) < 2:
            conf.Result('no')
            swig_location = None
        else:
            conf.Result('yes')
    return(swig_location)



# See if we need gfortran in order to build code with atlas/lapack
# Returns the full list of libraries used...
def CheckAtlasBuild(conf, libs):
   test_code = """
extern "C"{void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);}
int main(int argc, char *argv[]) { char C[1]; double D[1];int I[1];dgesvd_(C, C, I, I, D, I, D, D, I, D, I, D, I, I); }
"""

   conf.Message('Checking if ATLAS/LAPACK needs gfortran...')
   
   lastLIBS = conf.env['LIBS']
   conf.env.Append(LIBS = libs)
   result = conf.TryLink(test_code, '.cpp')
   if not result:
      conf.env.Append(LIBS = option_libs)
      result = conf.TryLink(test_code)
      if not result:
         conf.Result('error')
         return([])
      conf.Result('yes')
      conf.env.Replace(LIBS = lastLIBS)
      return(libs + 'gfortran')
   conf.Result('no')
   conf.env.Replace(LIBS = lastLIBS)
   return(libs)


# Check for existince of boost library with various naming variants
# Will return a tuple containing the correct name and a flag indicating
# whether this is the threaded or non-threaded version.
def CheckForBoostLibrary(conf, name, path, suffix):
   conf.Message('Checking for Boost library %s...' % name)

   if (os.path.isfile(os.path.join(path, 'lib%s-mt.%s' % (name , suffix)))):
      conf.Result('yes')
      return(name + '-mt', 1)

   if (os.path.isfile(os.path.join(path, 'lib%s.%s' % (name , suffix)))):
      conf.Result('yes')
      return(name, 0)

   def sortByLength(w1,w2):
      return len(w1)-len(w2)

    # Now check for names lib libboost_regex-gcc43-mt.so ...
   files = glob.glob(os.path.join(path, 'lib%s-*-mt.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result('yes')
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 1)

   files = glob.glob(os.path.join(path, 'lib%s-*.%s' % (name, suffix)))
   files.sort(cmp=sortByLength)
   if files:
      conf.Result('yes')
      name = os.path.basename(files[0])[3:-(len(suffix)+1)]
      return(name, 0)


   conf.Result('no')
   return('', -1)

            
