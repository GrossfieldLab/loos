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


# This file contains support code for SCons for building LOOS

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
import string
from distutils.version import LooseVersion
from sysconfig import get_paths, get_config_var

import SCons

import loos_build_config


default_lib_path = None
conda_path = None

# Attempt to canonicalize system name, type and other related info...
# Note: this exports to globals rather than being contained within the check
# framework.
def CheckSystemType(conf):
    conf.Message("Determining platform ...")
    loos_build_config.host_type = platform.system()

    # Detect CYGWIN & canonicalize linux type, setting defaults...
    if re.search("(?i)cygwin", loos_build_config.host_type):
        loos_build_config.host_type = "Cygwin"
        loos_build_config.suffix = "dll.a"
    elif loos_build_config.host_type == "Linux":
        # Determine linux variant...
        linux_type = platform.platform()

        if re.search("(?i)ubuntu", linux_type):
            linux_type = "ubuntu"
        elif re.search("(?i)suse", linux_type):
            linux_type = "suse"
        elif re.search("(?i)debian", linux_type):
            linux_type = "debian"
        elif re.search("(?i)centos", linux_type):
            linux_type = "centos"
        elif re.search("(?i)fedora", linux_type):
            linux_type = "fedora"

        loos_build_config.linux_type = linux_type

    # MacOS is special (of course...)
    elif loos_build_config.host_type == "Darwin":
        loos_build_config.suffix = "dylib"

    typemsg = loos_build_config.host_type
    if typemsg == "Linux":
        typemsg = typemsg + " [" + loos_build_config.linux_type + "]"

    conf.Result(typemsg)


# Create a revision file for linking against.
def setupRevision(env):

    # Divine the current revision...
    revision = loos_build_config.loos_version + " " + strftime("%y%m%d")

    # Now, write this out to a cpp file that can be linked in. This avoids
    # having to recompile everything when building on a new date.  We also
    # rely on SCons using the MD5 checksum to detect changes in the file
    # (even though it's always rewritten)
    revfile = open("src/revision.cpp", "w")
    revfile.write("#include <string>\n")
    revfile.write('std::string revision_label = "')
    revfile.write(revision)
    revfile.write('";\n')
    revfile.close()


# Let environment variables override or modify some build paramaters...
def environOverride(conf):
    # Allow overrides from environment...
    if "CXX" in os.environ:
        conf.env.Replace(CXX=os.environ["CXX"])
        print(("*** Using compiler " + os.environ["CXX"]))

    if "CCFLAGS" in os.environ:
        conf.env.Append(CCFLAGS=os.environ["CCFLAGS"])
        print(("*** Appending custom build flags: " + os.environ["CCFLAGS"]))

    if "LDFLAGS" in os.environ:
        conf.env.Append(LINKFLAGS=os.environ["LDFLAGS"])
        print(("*** Appending custom link flag: " + os.environ["LDFLAGS"]))


# Builder for setup scripts


def expand_scons_paths(path, topdir):
    newpath = []
    for item in path:
        item = item.replace("#", topdir + "/")
        newpath.append(item)
    return newpath


# This copies the environment setup script while changing the directory
# that's used for setting up PATH and [DY]LD_LIBRARY_PATH.  If LOOS
# is being built in a directory, the env script will be setup to use
# the built-in-place distribution.  If LOOS is being installed, then
# it will use the installation directory instead.


def script_builder_python(target, source, env):

    libpaths = env["LIBPATH"]
    libpaths.pop(0)

    cpppaths = env["CPPPATH"]
    cpppaths.pop(0)

    ldlibrary = list(loos_build_config.user_libdirs.values())

    if "install" not in SCons.Script.COMMAND_LINE_TARGETS:
        toolpath = "$LOOS/Tools:" + ":".join(
            [
                "$LOOS/Packages/" + s
                for s in [
                    loos_build_config.package_list[i]
                    for i in loos_build_config.package_list
                ]
            ]
        )
        loos_dir = env.Dir(".").abspath
        libpaths.insert(0, loos_dir)
        cpppaths.insert(0, loos_dir)
        ldlibrary.insert(0, loos_dir)

        libpaths = expand_scons_paths(libpaths, loos_dir)
        cpppaths = expand_scons_paths(cpppaths, loos_dir)
        ldlibrary = expand_scons_paths(ldlibrary, loos_dir)

        loos_pythonpath = loos_dir

    else:
        loos_dir = env["PREFIX"]
        toolpath = loos_dir + "/bin"
        libpaths.insert(0, loos_dir + "/lib")
        cpppaths.insert(0, loos_dir + "/include")
        ldlibrary.insert(0, loos_dir + "/lib")
        loos_pythonpath = loos_dir + "/lib"

    file = open(str(source[0]), "r")
    script = file.read()
    script_template = string.Template(script)
    script = script_template.substitute(
        loos_path=loos_dir,
        tool_path=toolpath,
        libpath=":".join(libpaths),
        cpppath=":".join(cpppaths),
        linkflags=env["LINKFLAGS"],
        libs=":".join(env["LIBS"]),
        ccflags=env["CCFLAGS"],
        loos_cxx=env["CXX"],
        loos_pythonpath=loos_pythonpath,
        ldlibrary=":".join(ldlibrary),
    )

    outfile = open(str(target[0]), "w")
    outfile.write(script)

    return None


# Verify that we have swig and it's v2.0+
# Returns the path to swig
def CheckForSwig(conf, min_version):
    conf.Message("Checking for Swig...")
    # Need to use has_key() for older distros...
    if "SWIGVERSION" in conf.env:
        if LooseVersion(conf.env["SWIGVERSION"]) >= LooseVersion(min_version):
            conf.Result("yes [%s]" % (conf.env["SWIGVERSION"]))
            return 1
        else:
            conf.Result(
                "too old [%s, requires at least %s; pyloos disabled]"
                % (conf.env["SWIGVERSION"], min_version)
            )
            return 0

    conf.Result("no [pyloos disabled]")
    return 0


# See if a library requires another to link...
def CheckAtlasRequires(conf, name, lib, required):

    conf.Message("Checking if %s requires %s ... " % (name, required))
    lastLIBS = list(conf.env["LIBS"])
    conf.env.Append(LIBS=lib)

    #test_code = """
    # extern "C"{void dgesvd_(char*, char*, int*, int*, double*, int*, double*,
    # double*, int*, double*, int*#, double*, int*, int*);}
    # int main(int argc, char *argv[]) { char C[1]; double D[1];int
    # I[1];dgesvd_(C, C, I, I, D, I, D, D, I, #D, I, D, I, I); }
    # """
    test_code = "int main(){return(0);}"

    result = conf.TryLink(test_code, ".cpp")
    if not result:
        conf.env.Append(LIBS=required)
        result = conf.TryLink(test_code, ".cpp")
        conf.env.Replace(LIBS=lastLIBS)
        if not result:
            conf.Result("fail")
            return ()
        conf.Result("yes")
        return (lib, required)

    conf.env.Replace(LIBS=lastLIBS)
    conf.Result("no")
    return lib


# Check for IEC-559 compliance
def CheckForIEC559(conf):

    conf.Message("Checking for IEC-559/IEE-754 support... ")

    test_code = """
#include <iostream>
#include <limits>

int main(int argc, char *argv[]) {
  if (std::numeric_limits<double>::is_iec559 && std::numeric_limits<float>::is_iec559)
     std::cout << "yes";
  else
     std::cout << "no";
}
"""

    result = conf.TryRun(test_code, ".cpp")
    if not result[0]:
        conf.Result("unable to check")
        return 0

    if result[1] != "yes":
        conf.Result("no")
        return 0

    conf.Result("yes")

    return 1


# Check for existince of boost library with various naming variants
# Will return a tuple containing the correct name and a flag indicating
# whether this is the threaded or non-threaded version.
# This will only search specified paths, not the built-in paths for g++,
# so some libraries may be missed...
def CheckForBoostLibrary(conf, name, path, suffix):
    conf.Message("Checking for Boost library %s..." % name)
    name = "boost_" + name

    def sortByLength(w):
        return len(w)

    # Now check for names lib libboost_regex-gcc43-mt.so ...
    files = glob.glob(os.path.join(path, "lib%s*-mt.%s" % (name, suffix)))
    files.sort(key=sortByLength)
    if files:
        conf.Result(name + "-mt")
        name = os.path.basename(files[0])[3 : -(len(suffix) + 1)]
        return (name, 1)

    files = glob.glob(os.path.join(path, "lib%s*.%s" % (name, suffix)))
    files.sort(key=sortByLength)
    if files:
        conf.Result(name)
        name = os.path.basename(files[0])[3 : -(len(suffix) + 1)]
        return (name, 0)

    conf.Result("missing")
    return ("", -1)


# Check for Boost include files...
def CheckBoostHeaders(conf):
    test_code = """
#include <boost/version.hpp>
int main(int argc, char *argv[]) { return(0); }
"""

    conf.Message("Checking for Boost... ")
    result = conf.TryLink(test_code, ".cpp")
    if not result:
        conf.Result("no")
        return 0

    conf.Result("yes")
    return 1


# Check for version of Boost includes
def CheckBoostHeaderVersion(conf, min_boost_version):
    source_code = """
#include <iostream>
#include <boost/version.hpp>
int main(int argc, char *argv[]) { std::cout << BOOST_LIB_VERSION; return(0); }
"""

    conf.Message("Checking Boost version... ")
    result = conf.TryRun(source_code, ".cpp")
    if not result[0]:
        conf.Result("boost missing or incomplete?")
        return 0
    boost_version = result[1]
    loos_build_config.versions["boost"] = boost_version

    if LooseVersion(boost_version) < LooseVersion(min_boost_version):
        conf.Result(
            "%s [too old, LOOS requires at least %s]"
            % (boost_version, min_boost_version)
        )
        return 0

    conf.Result("%s [ok]" % boost_version)
    return 1


# Check for presence of a directory
def CheckDirectory(conf, dirname):

    conf.Message("Checking for directory %s... " % dirname)
    if os.path.isdir(dirname):
        conf.Result("yes")
        return 1
    conf.Result("no")
    return 0


def CheckNumpy(conf, pythonpath):
    global default_lib_path
    conf.Message("Checking for numpy... ")

    if conf.env.USING_CONDA:
        python_lib_dir = get_config_var("MACHDESTLIB")
        # some older compilers require this directory explicitly in the
        # include path
        numpy_inc_path = python_lib_dir + "/site-packages/numpy/core/include"
        if CheckDirectory(conf, numpy_inc_path):
            conf.env.Append(CPPPATH=numpy_inc_path)
            conf.Result("yes")
            return 1
        else:
            print("Numpy not found inside conda, looking elsewhere...")

    env = conf.env["ENV"]

    ok = checkForPythonHeader(conf, "numpy/arrayobject.h")
    if ok:
        conf.Result("yes")
        return 1
    newpaths = []

    if "PYTHON_PATH" in conf.env:
        envpath = conf.env["PYTHON_PATH"]
        # Catch cases where PYTHON_PATH is present but null...
        if len(envpath) > 1:
            newpaths.extend(envpath.split(":"))

    newpaths.append(default_lib_path)
    for dir in newpaths:
        for p, d, f in os.walk(dir):
            for file in f:
                if file == "arrayobject.h":
                    (prefix, numpydir) = os.path.split(p)
                    ok = checkForPythonHeaderInPath(
                        conf, "numpy/arrayobject.h", [prefix]
                    )
                    if ok:
                        conf.Result("yes")
                        return 1

    # Special handling for MacOS
    if loos_build_config.host_type == "Darwin":
        ok = checkForPythonHeaderInPath(
            conf,
            "numpy/arrayobject.h",
            [
"/System/Library/Frameworks/Python.framework/Versions/Current/Extras/lib/python/numpy/core/include"
            ],
        )
        if ok:
            conf.Result("yes")
            return 1

    conf.Result("no")
    return 0


def SetupOpenBabelPaths(env):
    """
    If you supplied OPENBABEL_INCLUDE or OPENBABEL_LIBPATH, use them.
    If you didn't but did supply OPENBABEL, derive include and lib paths
    from it.
    If we're in conda and you didn't supply anything, build what we think
    conda wants.
    Otherwise, do nothing.

    """

    #OPENBABEL = env["OPENBABEL"]
    #OPENBABEL_INCLUDE = env["OPENBABEL_INCLUDE"]
    #OPENBABEL_LIBPATH = env["OPENBABEL_LIBPATH"]

    openbabel_libpath = ""
    openbabel_include = ""

    if "OPENBABEL_INCLUDE" in env:
        openbabel_include = env["OPENBABEL_INCLUDE"]
    if "OPENBABEL_LIBPATH" in env:
        openbabel_libpath = env["OPENBABEL_LIBPATH"]

    if "OPENBABEL" in env:
        if not openbabel_include:
            openbabel_include = os.path.join(env["OPENBABEL"], "include", "openbabel3")
        if not openbabel_libpath:
            openbabel_libpath = os.path.join(OPENBABEL, "lib")

    if "CONDA_PREFIX" in env:
        if not openbabel_include:
            openbabel_include = os.path.join(env["CONDA_PREFIX"],
                                        "include", "openbabel3")
        if not openbabel_libpath:
            openbabel_libpath = os.path.join(env["CONDA_PREFIX"], "lib")

    if openbabel_libpath:
        env.Prepend(LIBPATH=[openbabel_libpath])
        env["OPENBABEL_LIBPATH"] = openbabel_libpath

    if openbabel_include:
        env.Prepend(CPPPATH=[openbabel_include])
        env["OPENBABEL_INCLUDE"] = openbabel_include

def SetupBoostPaths(env):

    BOOST = env["BOOST"]
    BOOST_INCLUDE = env["BOOST_INCLUDE"]
    BOOST_LIBPATH = env["BOOST_LIBPATH"]
    BOOST_LIBS = env["BOOST_LIBS"]

    # If boost is not set but we're inside a conda environment,
    # automatically redirect boost into here...
    if not BOOST and env.USING_CONDA:
        BOOST = env["CONDA_PREFIX"]

    boost_libpath = ""
    boost_include = ""

    if BOOST:
        boost = BOOST
        boost_include = boost + "/include"
        boost_libpath = boost + "/lib"
        loos_build_config.user_libdirs["BOOST"] = boost_libpath
        loos_build_config.user_boost_flag = 1

    if BOOST_INCLUDE:
        boost_include = BOOST_INCLUDE
    if BOOST_LIBPATH:
        boost_libpath = BOOST_LIBPATH
        loos_build_config.user_libdirs["BOOST"] = boost_libpath
        loos_build_config.user_boost_flag = 1

    if boost_libpath:
        env.Prepend(LIBPATH=[boost_libpath])
        env["BOOST_LIBPATH"] = boost_libpath
    if boost_include:
        env.Prepend(CPPPATH=[boost_include])
        env["BOOST_INCLUDE"] = boost_include


def SetupNetCDFPaths(env):
    NETCDF = env["NETCDF"]
    NETCDF_INCLUDE = env["NETCDF_INCLUDE"]
    NETCDF_LIBPATH = env["NETCDF_LIBPATH"]
    NETCDF_LIBS = env["NETCDF_LIBS"]

    # If netcdf is not set but we're inside a conda environment,
    # automatically redirect netcdf into here...
    if not NETCDF and env.USING_CONDA:
        NETCDF=env['CONDA_PREFIX']

    netcdf_libpath = ""
    netcdf_include = ""

    if NETCDF:
        netcdf = NETCDF
        netcdf_include = netcdf + "/include"
        netcdf_libpath = netcdf + "/lib"
        loos_build_config.user_libdirs["NETCDF"] = netcdf_libpath

    if NETCDF_INCLUDE:
        netcdf_include = NETCDF_INCLUDE
    if NETCDF_LIBPATH:
        netcdf_libpath = NETCDF_LIBPATH
        loos_build_config.user_libdirs["NETCDF"] = netcdf_libpath

    if netcdf_libpath:
        env.Prepend(LIBPATH=[netcdf_libpath])
    if netcdf_include:
        env.Prepend(CPPPATH=[netcdf_include])


def AutoConfigSystemBoost(conf):
    boost_libs = []
    first = 1
    thread_suffix = 0

    for libname in loos_build_config.required_boost_libraries:
        if first:
            first = 0
            full_libname = "boost_" + libname + "-mt"
            result = conf.CheckLib(full_libname, autoadd=0)
            if result:
                boost_libs.append(full_libname)
                thread_suffix = 1
            else:
                full_libname = "boost_" + libname
                result = conf.CheckLib(full_libname, autoadd=0)
                if result:
                    boost_libs.append(full_libname)
        else:
            full_libname = "boost_" + libname
            if thread_suffix:
                full_libname += "-mt"
            result = conf.CheckLib(full_libname, autoadd=0)
            if result:
                boost_libs.append(full_libname)
            else:
                print(("Error- missing Boost library %s" % libname))
                conf.env.Exit(1)

    return boost_libs


def AutoConfigUserBoost(conf):
    boost_libs = []
    first = 1
    thread_suffix = 0

    for libname in loos_build_config.required_boost_libraries:
        result = conf.CheckForBoostLibrary(
            libname, conf.env["BOOST_LIBPATH"], loos_build_config.suffix
        )
        if not result[0]:
            print(("Error- missing Boost library %s" % libname))
            conf.env.Exit(1)
        if first:
            thread_suffix = result[1]
        else:
            if thread_suffix and not result[1]:
                print(("Error- expected %s-mt but found %s" % (libname, libname)))
                conf.env.Exit(1)
            elif not thread_suffix and result[1]:
                print(("Error- expected %s but found %s-mt" % (libname, libname)))
                conf.env.Exit(1)
        boost_libs.append(result[0])

    return boost_libs


# Check to see if a function exists with current libs
# If gfortran is present and function is not found, then
# try appending gfortran to lib list...


def checkForFunction(context, funcname, libs, has_gfortran):
    old_libs = list(context.env["LIBS"])
    context.env.Append(LIBS=libs)
    requires_gf = 0

    ok = context.CheckFunc(funcname)
    if not ok and has_gfortran:
        print("Trying again with gfortran...")
        context.env.Append(LIBS=["gfortran"])
        ok = context.CheckFunc(funcname)
        if ok:
            requires_gf = 1

    context.env["LIBS"] = old_libs
    return (ok, requires_gf)


# Tries adding each library from list to the build lib list
# and sees if funcname is present.  Any lib in the excludelist
# is ignored.  No special handling for gfortran...


def checkLibsForFunction(context, funcname, liblist, excludelist):
    for lib in liblist:
        if lib in excludelist:
            continue
        old_libs = list(context.env["LIBS"])
        context.env.Append(LIBS=lib)
        print(("> Checking in %s ..." % lib))
        ok = context.CheckFunc(funcname)
        context.env["LIBS"] = old_libs
        if ok:
            return lib
    return ""


def checkForPythonHeader(context, header):
    test_code = (
        """
#include <Python.h>
#include <%s>
"""
        % header
    )

    oldcpp = None
    if "CPPFLAGS" in context.env:
        oldcpp = context.env["CPPFLAGS"]
    if "CPPPATH" in context.env and not ("CONDA_PREFIX" in context.env):
        for dir in context.env["CPPPATH"]:
            context.env.Append(CPPFLAGS="-I%s " % dir)
            print("CPPFLAGS", context.env["CPPFLAGS"])

    ok = context.TryCompile(test_code, ".cpp")

    if oldcpp:
        context.env["CPPFLAGS"] = oldcpp

    return ok


def checkForPythonHeaderInPath(context, header, pathlist):

    for path in pathlist:
        oldcpp = None
        if "CPPPATH" in context.env:
            oldcpp = context.env["CPPPATH"]
        context.env.Append(CPPPATH=[path])
        ok = checkForPythonHeader(context, header)
        if ok:
            return True
        if oldcpp:
            context.env["CPPPATH"] = oldcpp
    return False


def AutoConfiguration(env):
    global default_lib_path
    global conda_path

    conf = env.Configure(
        custom_tests={
            "CheckForSwig": CheckForSwig,
            "CheckBoostHeaders": CheckBoostHeaders,
            "CheckForBoostLibrary": CheckForBoostLibrary,
            "CheckBoostHeaderVersion": CheckBoostHeaderVersion,
            "CheckDirectory": CheckDirectory,
            "CheckAtlasRequires": CheckAtlasRequires,
            "CheckForIEC559": CheckForIEC559,
            "CheckSystemType": CheckSystemType,
            "CheckNumpy": CheckNumpy,
        }
    )

    use_threads = int(env["threads"])

    # Get system information
    conf.CheckSystemType()

    conf.env["host_type"] = loos_build_config.host_type
    conf.env["linux_type"] = loos_build_config.linux_type

    if env.GetOption("clean") or env.GetOption("help"):
        env["HAS_NETCDF"] = 1
    else:
        has_netcdf = 0

        if env.USING_CONDA:
            conda_path = env["CONDA_PREFIX"]
            default_lib_path = conda_path + "/lib"
            if loos_build_config.host_type != "Darwin":
                conf.env.Append(RPATH=default_lib_path)
        else:
            default_lib_path = "/usr/lib"

            # if we're not in conda, add system library directory
            if not conf.CheckDirectory("/usr/lib64"):
                if not conf.CheckDirectory("/usr/lib"):
                    print("Fatal error- cannot find your system library directory")
                    conf.env.Exit(1)
            else:
                # /usr/lib64 is found, so make sure we link against this
                # (and not against any 32-bit libs)
                default_lib_path = "/usr/lib64"
        conf.env.Append(LIBPATH=default_lib_path)

        # Only setup ATLAS if we're not on a Mac and we're not using conda
        if loos_build_config.host_type != "Darwin" and not env.USING_CONDA:
            atlas_libpath = ""
            ATLAS_LIBPATH = env["ATLAS_LIBPATH"]
            ATLAS_LIBS = env["ATLAS_LIBS"]
            if not ATLAS_LIBPATH:
                # Some distros have atlas in /atlas-base, so must check that...
                if conf.CheckDirectory(default_lib_path + "/atlas-base"):
                    atlas_libpath = default_lib_path + "/atlas-base"
                elif conf.CheckDirectory(default_lib_path + "/atlas"):
                    atlas_libpath = default_lib_path + "/atlas"
                else:
                    print("Warning: Could not find an atlas directory! ")
            else:
                atlas_libpath = ATLAS_LIBPATH
                loos_build_config.user_libdirs["ATLAS"] = atlas_libpath

            if atlas_libpath:
                conf.env.Prepend(LIBPATH=[atlas_libpath])

        if not conf.CheckLib("pthread"):
            print("Error- LOOS requires a pthread library installed")

        # Now that we know the default library path, setup Boost, NetCDF, and
        # ATLAS based on the environment or custom.py file
        SetupBoostPaths(conf.env)
        SetupNetCDFPaths(conf.env)
        SetupOpenBabelPaths(conf.env)

        # Check for standard typedefs...
        if not conf.CheckType("ulong", "#include <sys/types.h>\n"):
            conf.env.Append(CCFLAGS="-DREQUIRES_ULONG")
        if not conf.CheckType("uint", "#include <sys/types.h>\n"):
            conf.env.Append(CCFLAGS="-DREQUIRES_UINT")

        # Check for floating point format...
        if not conf.CheckForIEC559():
            print("Error- your system must use the IEC559/IEEE754 floating point")
            print("       format for Gromacs support in LOOS.  Check your compiler")
            print("       options or contact the LOOS developers at")
            print("       loos.maintainer@gmail.com")
            conf.env.Exit(1)

        # --- NetCDF Autoconf
        has_netcdf = 0
        if conf.env["NETCDF_LIBS"]:
            netcdf_libs = env["NETCDF_LIBS"]
            conf.env.Append(CCFLAGS=["-DHAS_NETCDF"])
            has_netcdf = 1
        else:
            if conf.CheckLibWithHeader(
                "netcdf", "netcdf.h", "c"
            ):  # Should we check C or C++?
                netcdf_libs = "netcdf"
                conf.env.Append(CCFLAGS=["-DHAS_NETCDF"])
                has_netcdf = 1

        conf.env["HAS_NETCDF"] = has_netcdf

        # --- Swig Autoconf (unless user requested NO PyLOOS)
        if int(env["pyloos"]):
            if conf.CheckForSwig(loos_build_config.min_swig_version):
                conf.env["pyloos"] = 1
                pythonpath = get_paths()['include']
                if "PYTHON_INC" in conf.env:
                    if conf.env["PYTHON_INC"] != "":
                        pythonpath = conf.env["PYTHON_INC"]

                conf.env.Append(CPPPATH=[pythonpath])
                if not conf.CheckNumpy(pythonpath):
                    print("ERROR- PyLOOS build requires NumPy")
                    conf.env.Exit(1)
            else:
                conf.env["pyloos"] = 0

        # --- Boost Autoconf
        if not conf.CheckBoostHeaders():
            conf.env.Exit(1)

        if not conf.CheckBoostHeaderVersion(loos_build_config.min_boost_version):
            conf.env.Exit(1)

        if conf.env["BOOST_LIBS"]:
            boost_libs = env.Split(env["BOOST_LIBS"])
        if env.USING_CONDA:
            boost_libs = AutoConfigUserBoost(conf)
        elif not loos_build_config.user_boost_flag:
            boost_libs = AutoConfigSystemBoost(conf)
        else:
            boost_libs = AutoConfigUserBoost(conf)

        env.Append(LIBS=boost_libs)

        # --- Check for ATLAS/LAPACK and how to build

        if loos_build_config.host_type != "Darwin" and not env.USING_CONDA:
            atlas_libs = ""  # List of numerics libs required for LOOS

            if env["ATLAS_LIBS"]:
                atlas_libs = env.Split(env["ATLAS_LIBS"])
            else:

                numerics = {
                    "openblas": 0,
                    "satlas": 0,
                    "atlas": 0,
                    "lapack": 0,
                    "f77blas": 0,
                    "cblas": 0,
                    "blas": 0,
                }

                if use_threads:
                    numerics["tatlas"] = 0
                    numerics["ptcblas"] = 0
                    numerics["ptf77blas"] = 0

                for libname in numerics:
                    if conf.CheckLib(libname, autoadd=0):
                        numerics[libname] = 1

                atlas_libs = []
                atlas_name = ""

                has_gfortran = 0
                if conf.CheckLib("gfortran", autoadd=0):
                    has_gfortran = 1

                if use_threads and numerics["tatlas"]:
                    atlas_libs.append("tatlas")
                    atlas_name = "tatlas"
                elif numerics["satlas"]:
                    atlas_libs.append("satlas")
                    atlas_name = "satlas"
                else:

                    if numerics["lapack"]:
                        atlas_libs.append("lapack")

                    if use_threads and (numerics["ptf77blas"] and numerics["ptcblas"]):
                        atlas_libs.extend(["ptf77blas", "ptcblas"])
                    elif numerics["f77blas"] and numerics["cblas"]:
                        atlas_libs.extend(["f77blas", "cblas"])
                    elif numerics["blas"]:
                        atlas_libs.append("blas")
                    else:
                        print("Error- you must have some kind of blas installed")
                        conf.env.Exit(1)

                    if numerics["atlas"]:
                        atlas_libs.append("atlas")
                        atlas_name = "atlas"

                # Try to figure out how to build with ATLAS...
                # We need these functions, so find a combination of libs and
                # libpaths will work...
                for funcname in ("dgesvd_", "dgemm_", "dtrmm_", "dsyev_"):
                    (ok, requires_gfortran) = checkForFunction(
                        conf, funcname, atlas_libs, has_gfortran
                    )
                    if requires_gfortran:
                        print("Build Requires gfortran")
                        atlas_libs.append("gfortran")

                    if not ok:
                        lib = checkLibsForFunction(
                            conf, funcname, list(numerics.keys()), atlas_libs
                        )
                        if lib:
                            atlas_libs.insert(0, lib)
                        else:
                            # Try putting scanning default_lib_path
                            # first...SUSE requires
                            # the lapack in /usr/lib first...
                            print(
                                (
                                    "Searching %s first for libraries..."
                                    % default_lib_path
                                )
                            )
                            # Remove the default_lib_path from the list and
                            # prepend...
                            libpaths = list(conf.env["LIBPATH"])
                            libpaths.remove(default_lib_path)
                            libpaths.insert(0, default_lib_path)
                            conf.env["LIBPATH"] = libpaths
                            (ok, requires_gfortran) = checkForFunction(
                                conf, funcname, atlas_libs, has_gfortran
                            )
                            if requires_gfortran:
                                print("Build requires gfortran")
                                atlas_libs.append("gfortran")

                            if not ok:
                                lib = checkLibsForFunction(
                                    conf, funcname, list(numerics.keys()), atlas_libs
                                )
                                if lib:
                                    atlas_libs.insert(0, lib)
                                else:
                                    print(
                                        "Error- could not figure out where ", funcname, " is located.",
                                    )
                                    print(
                                        "Try manually specifying ATLAS_LIBS and ATLAS_LIBPATH"
                                    )
                                    conf.env.Exit(1)

            # Hack to extend list rather than append a list into a list
            for lib in atlas_libs:
                conf.env.Append(LIBS=lib)
        elif env.USING_CONDA:
            conf.env.Append(LIBS="openblas")

        # Suppress those annoying maybe used unitialized warnings that -Wall
        # gives us...
        ccflags = conf.env["CCFLAGS"]
        conf.env.Append(
            CCFLAGS=["-Wno-maybe-uninitialized", "-Werror"]
        )  # Try suppressing, make bad flags an error
        ok = conf.TryCompile("", ".c")
        conf.env["CCFLAGS"] = ccflags
        if ok:
            conf.env.Append(CCFLAGS=["-Wno-maybe-uninitialized"])

        environOverride(conf)
        if "LIBS" in conf.env:
            print(
                "Autoconfigure will use these libraries to build LOOS:\n\t",
                conf.env["LIBS"],
            )
        if "LIBPATH" in conf.env:
            print(
                "Autoconfigure will add the following directories to find libs:\n\t",
                conf.env["LIBPATH"],
            )
        env = conf.Finish()


#########################################################################################3


def addDeprecatedOptions(opt):
    from SCons.Variables import PathVariable

    opt.Add(PathVariable("LAPACK", "Path to LAPACK", "",
PathVariable.PathAccept))
    opt.Add(PathVariable("ATLAS", "Path to ATLAS", "", PathVariable.PathAccept))
    opt.Add(
        PathVariable("ATLASINC", "Path to ATLAS includes", "",
PathVariable.PathAccept)
    )
    opt.Add(
        PathVariable("BOOSTLIB", "Path to BOOST libraries", "",
PathVariable.PathAccept)
    )
    opt.Add(
        PathVariable("BOOSTINC", "Path to BOOST includes", "",
PathVariable.PathAccept)
    )
    opt.Add("BOOSTREGEX", "Boost regex library name", "")
    opt.Add("BOOSTPO", "Boost program options library name", "")
    opt.Add(
        PathVariable(
            "LIBXTRA", "Path to additional libraries", "",
PathVariable.PathAccept
        )
    )
    opt.Add(
        PathVariable(
            "NETCDFINC", "Path to netcdf include files", "",
PathVariable.PathAccept
        )
    )
    opt.Add(
        PathVariable(
            "NETCDFLIB", "Path to netcdf library files", "",
PathVariable.PathAccept
        )
    )
    opt.Add(
        PathVariable(
            "ALTPATH", "Additional path to commands", "",
PathVariable.PathAccept
        )
    )
    opt.Add(
        PathVariable(
            "LIBS_OVERRIDE", "Override linked libs", "", PathVariable.PathAccept
        )
    )
    opt.Add(
        PathVariable(
            "LIBS_PATHS_OVERRIDE", "Override paths to libs", "",
PathVariable.PathAccept
        )
    )


def makeDeprecatedVariableWarning():
    state = {"warned": 0}

    def warning(what, mapto):
        if not state["warned"]:
            state["warned"] = 1
            print(
                """
***WARNING***
You are using old-style (deprecated) variables either
on the command line or in your custom.py file.  These
will be ignored.  The following deprecated variables
are set,
"""
            )
        print("\t%s: %s" % (what, mapto))

    return warning


def checkForDeprecatedOptions(env):
    mapping = {
        "LAPACK": "use ATLAS_LIBPATH",
        "ATLAS": "use ATLAS_LIBPATH",
        "ATLASINC": "no replacement",
        "BOOSTLIB": "use BOOST_LIBPATH or BOOST",
        "BOOSTINC": "use BOOST_INCLUDE or BOOST",
        "BOOSTREGEX": "use BOOST_LIBS",
        "BOOSTPO": "use BOOST_LIBS",
        "LIBXTRA": "use ATLAS_LIBS, BOOST_LIBS, or NETCDF_LIBS",
        "NETCDFINC": "use NETCDF_INCLUDE or NETCDF",
        "NETCDFLIB": "use NETCDF_LIBPATH or NETCDF",
        "ALTPATH": "Set your shell PATH instead",
        "LIBS_OVERRIDE": "use ATLAS_LIBS, BOOST_LIBS, or NETCDF_LIBS",
        "LIBS_PATHS_OVERRIDE": "use ATLAS_LIBPATH, BOOST_LIBPATH, or NETCDF_LIBPATH",
    }

    warner = makeDeprecatedVariableWarning()
    for name in mapping:
        if name in env:
            if env[name]:
                warner(name, mapping[name])
