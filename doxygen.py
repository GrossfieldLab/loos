# vim: set et sw=3 tw=0 fo=awqorc ft=python:
# -*- mode:python; coding:utf-8; -*-
#
# Astxx, the Asterisk C++ API and Utility Library.
# Copyright © 2005, 2006  Matthew A. Nicholson
# Copyright © 2006  Tim Blechmann
#
#  Copyright © 2007 Christoph Boehme
#
#  Copyright © 2012 Dirk Baechle
#
#  Copyright © 2013 Russel Winder
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License version 2.1 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os
import os.path
import glob
from fnmatch import fnmatch

# Currently supported output formats and their default
# values and output locations.
# From left to right:
#   1. default setting YES|NO
#   2. default output folder for this format
#   3. name of the (main) output file
#   4. default extension "
#   5. field for overriding the output file extension
output_formats = {
    "HTML": ("YES", "html", "index", ".html", "HTML_FILE_EXTENSION"),
    "LATEX": ("YES", "latex", "refman", ".tex", ""),
    "RTF": ("NO", "rtf", "refman", ".rtf", ""),
    "MAN": ("NO", "man", "", ".3", "MAN_EXTENSION"),
    "XML": ("NO", "xml", "index", ".xml", ""),
}


def DoxyfileParse(file_contents, conf_dir, data=None):
    """
    Parse a Doxygen source file and return a dictionary of all the values.
    Values will be strings and lists of strings.
    """
    if data is None:
        data = {}

    import shlex
    lex = shlex.shlex(instream=file_contents, posix=True)
    lex.wordchars += "*+./-:@"
    lex.whitespace = lex.whitespace.replace("\n", "")
    lex.escape = ""

    lineno = lex.lineno
    token = lex.get_token()
    key = None
    last_token = ""
    key_token = True  # The first token should be a key.
    next_key = False
    new_data = True

    def append_data(data, key, new_data, token):
        if new_data or len(data[key]) == 0:
            data[key].append(token)
        else:
            data[key][-1] += token

    while token:
        if token in ['\n']:
            if last_token not in ['\\']:
                key_token = True
        elif token in ['\\']:
            pass
        elif key_token:
            key = token
            key_token = False
        else:
            if token == "+=":
                if key not in data:
                    data[key] = []
            elif token == "=":
                if key == "TAGFILES" and key in data:
                    append_data(data, key, False, "=")
                    new_data = False
                elif key == "@INCLUDE" and key in data:
                    # don't reset the @INCLUDE list when we see a new @INCLUDE line.
                    pass
                else:
                    data[key] = []
            elif key == "@INCLUDE":
                # special case for @INCLUDE key: read the referenced
                # file as a doxyfile too.
                nextfile = token
                if not os.path.isabs(nextfile):
                    nextfile = os.path.join(conf_dir, nextfile)
                if nextfile in data[key]:
                    raise Exception("recursive @INCLUDE in Doxygen config: " + nextfile)
                data[key].append(nextfile)
                fh = open(nextfile, 'r')
                DoxyfileParse(fh.read(), conf_dir, data)
                fh.close()
            else:
                append_data(data, key, new_data, token)
                new_data = True

        last_token = token
        token = lex.get_token()

        if last_token == '\\' and token != '\n':
            new_data = False
            append_data(data, key, new_data, '\\')

    # compress lists of len 1 into single strings
    for (k, v) in data.items():
        if len(v) == 0:
            data.pop(k)

        # items in the following list will be kept as lists and not converted to strings
        if k in ["INPUT", "FILE_PATTERNS", "EXCLUDE_PATTERNS", "TAGFILES", "@INCLUDE"]:
            continue

        if len(v) == 1:
            data[k] = v[0]

    return data


def DoxySourceFiles(node, env):
    """
    Scan the given node's contents (a Doxygen file) and add
    any files used to generate docs to the list of source files.
    """
    default_file_patterns = [
        '*.c', '*.cc', '*.cxx', '*.cpp', '*.c++', '*.java', '*.ii', '*.ixx',
        '*.ipp', '*.i++', '*.inl', '*.h', '*.hh ', '*.hxx', '*.hpp', '*.h++',
        '*.idl', '*.odl', '*.cs', '*.php', '*.php3', '*.inc', '*.m', '*.mm',
        '*.py',
    ]

    default_exclude_patterns = [
        '*~',
    ]

    sources = []

    # We're running in the top-level directory, but the doxygen
    # configuration file is in the same directory as node; this means
    # that relative pathnames in node must be adjusted before they can
    # go onto the sources list
    conf_dir = os.path.dirname(str(node))

    data = DoxyfileParse(node.get_contents(), conf_dir)

    if data.get("RECURSIVE", "NO") == "YES":
        recursive = True
    else:
        recursive = False

    file_patterns = data.get("FILE_PATTERNS", default_file_patterns)
    exclude_patterns = data.get("EXCLUDE_PATTERNS", default_exclude_patterns)

    input = data.get("INPUT")
    if input:
        for node in data.get("INPUT", []):
            if not os.path.isabs(node):
                node = os.path.join(conf_dir, node)
            if os.path.isfile(node):
                sources.append(node)
            elif os.path.isdir(node):
                if recursive:
                    for root, dirs, files in os.walk(node):
                        for f in files:
                            filename = os.path.join(root, f)

                            pattern_check = reduce(lambda x, y: x or bool(fnmatch(filename, y)), file_patterns, False)
                            exclude_check = reduce(lambda x, y: x and fnmatch(filename, y), exclude_patterns, True)

                            if pattern_check and not exclude_check:
                                sources.append(filename)
                else:
                    for pattern in file_patterns:
                        sources.extend(glob.glob("/".join([node, pattern])))
    else:
        # No INPUT specified, so apply plain patterns only
        if recursive:
            for root, dirs, files in os.walk('.'):
                for f in files:
                    filename = os.path.join(root, f)

                    pattern_check = reduce(lambda x, y: x or bool(fnmatch(filename, y)), file_patterns, False)
                    exclude_check = reduce(lambda x, y: x and fnmatch(filename, y), exclude_patterns, True)

                    if pattern_check and not exclude_check:
                        sources.append(filename)
        else:
            for pattern in file_patterns:
                sources.extend(glob.glob(pattern))

    # Add @INCLUDEd files to the list of source files:
    for node in data.get("@INCLUDE", []):
        sources.append(node)

    # Add tagfiles to the list of source files:
    for node in data.get("TAGFILES", []):
        file = node.split("=")[0]
        if not os.path.isabs(file):
            file = os.path.join(conf_dir, file)
        sources.append(file)

    # Add additional files to the list of source files:
    def append_additional_source(option, formats):
        for f in formats:
            if data.get('GENERATE_' + f, output_formats[f][0]) == "YES":
                file = data.get(option, "")
                if file != "":
                    if not os.path.isabs(file):
                        file = os.path.join(conf_dir, file)
                    if os.path.isfile(file):
                        sources.append(file)
                break

    append_additional_source("HTML_STYLESHEET", ['HTML'])
    append_additional_source("HTML_HEADER", ['HTML'])
    append_additional_source("HTML_FOOTER", ['HTML'])

    return sources


def DoxySourceScan(node, env, path):
    """
    Doxygen Doxyfile source scanner.  This should scan the Doxygen file and add
    any files used to generate docs to the list of source files.
    """
    filepaths = DoxySourceFiles(node, env)
    sources = map(lambda path: env.File(path), filepaths)
    return sources


def DoxySourceScanCheck(node, env):
    """Check if we should scan this file"""
    return os.path.isfile(node.path)


def DoxyEmitter(target, source, env):
    """Doxygen Doxyfile emitter"""
    doxy_fpath = str(source[0])
    conf_dir = os.path.dirname(doxy_fpath)
    data = DoxyfileParse(source[0].get_contents(), conf_dir)

    targets = []
    out_dir = data.get("OUTPUT_DIRECTORY", ".")
    if not os.path.isabs(out_dir):
        out_dir = os.path.join(conf_dir, out_dir)

    # add our output locations
    for (k, v) in output_formats.items():
        if data.get("GENERATE_" + k, v[0]) == "YES":
            # Initialize output file extension for MAN pages
            if k == 'MAN':
                # Is the given extension valid?
                manext = v[3]
                if v[4] and v[4] in data:
                    manext = data.get(v[4])
                # Try to strip off dots
                manext = manext.replace('.', '')
                # Can we convert it to an int?
                try:
                    e = int(manext)
                except:
                    # No, so set back to default
                    manext = "3"

                od = env.Dir(os.path.join(out_dir, data.get(k + "_OUTPUT", v[1]), "man" + manext))
            else:
                od = env.Dir(os.path.join(out_dir, data.get(k + "_OUTPUT", v[1])))
            # don't clobber target folders
            env.Precious(od)
            # set up cleaning stuff
            env.Clean(od, od)

            # Add target files
            if k != "MAN":
                # Is an extension override var given?
                if v[4] and v[4] in data:
                    fname = v[2] + data.get(v[4])
                else:
                    fname = v[2] + v[3]
                of = env.File(os.path.join(out_dir, data.get(k + "_OUTPUT", v[1]), fname))
                targets.append(of)
                # don't clean single files, we remove the complete output folders (see above)
                env.NoClean(of)
            else:
                # Special case: MAN pages
                # We have to add a target file docs/man/man3/foo.h.3
                # for each input file foo.h, so we scan the config file
                # a second time... :(
                filepaths = DoxySourceFiles(source[0], env)
                for f in filepaths:
                    if os.path.isfile(f) and f != doxy_fpath:
                        of = env.File(os.path.join(out_dir,
                                                   data.get(k + "_OUTPUT", v[1]),
                                                   "man" + manext,
                                                   f + "." + manext))
                        targets.append(of)
                        # don't clean single files, we remove the complete output folders (see above)
                        env.NoClean(of)

    # add the tag file if neccessary:
    tagfile = data.get("GENERATE_TAGFILE", "")
    if tagfile != "":
        if not os.path.isabs(tagfile):
            tagfile = os.path.join(conf_dir, tagfile)
        targets.append(env.File(tagfile))

    return (targets, source)


def generate(env):
    """
    Add builders and construction variables for the
    Doxygen tool.  This is currently for Doxygen 1.4.6.
    """
    doxyfile_scanner = env.Scanner(
        DoxySourceScan,
        "DoxySourceScan",
        scan_check=DoxySourceScanCheck,
    )

    import SCons.Builder
    doxyfile_builder = SCons.Builder.Builder(
        action="cd ${SOURCE.dir}  &&  ${DOXYGEN} ${SOURCE.file}",
        emitter=DoxyEmitter,
        target_factory=env.fs.Entry,
        single_source=True,
        source_scanner=doxyfile_scanner,
    )

    env.Append(BUILDERS={
        'Doxygen': doxyfile_builder,
    })

    env.AppendUnique(
        DOXYGEN='doxygen',
    )


def exists(env):
    """
    Make sure doxygen exists.
    """
    return env.Detect("doxygen")
