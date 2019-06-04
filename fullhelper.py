#!/usr/bin/env python3
from sys import argv
import textwrap

usage = "path/to/loos/fullhelper.py thefullhelpmessage.txt thetool.cpp '\"the locator string\"' [desired_num_cols]"
fullhelper_fullhelp_message = """
This tool helps you wrap fullhelp messages for loos C++ source so that 
they will be formatted correctly in the WYSIWYG loos style, where all the 
needed punctuation and whitespace are specified in source. The objective is to
provide a way for users to write their desired fullhelp message in softwrapping
text editors, fairly common these days, as a plain text file. The correct
wrapping and line-decorating are then done by this script, and the adjusted
text is then written into the source at a spot specified by a locator or signal
string. This is done by assembling the entire updated source code as a list, 
then substituting it into a copy of the original source code (in case something
bad happens). You can check this file, which will have the same name as your
tool's source but prefixed with 'fullhelp-test', and if you like what you see
renaming the fullhelp-test file to be the source file. The only difference
should be the replacement of your signal string with the fullhelp string.

The intention of this script is for you to write your source code as you would,
including a placeholder string instead of your fullhelpstring in the definition
of 'string fullHelpMessage(void)' such as:
string fullHelpMessage(void) {
  string msg =
    "XXX";

  return(msg);
}

In this way you can continue to work on your code, able to build it and even
test the appearance of the options framework help phrases, without having to
come up with and correctly format a handsome and useful full help until the 
end. Obviously you can use any string you want, but make good choices. This
script will find and replace only the first instance of the signal string.  
You need to include its surrounding quotes (easily done in bash by 
including the signal string inside of doublequotes nested in single quotes, see
usage string). You can also make tweaks to the fullhelp source (the soft-wrapped
text file), then see what they look like by up-arrow rerunning your fullhelper
command at your terminal repeatedly.

The three obligatory arguments are the name of .txt containing your full help
message, the source code file you'd like to sub your message into, and the
locator string you are using as a target for replacement by the full help 
message. The fourth, optional, argument is the column count to wrap to. The
default column count is 80.
"""
# number of columns to wrap to
if '-h' in argv:
    print(fullhelper_fullhelp_message)
    print(usage)
    exit(0)
if '--fullhelp' in argv:
    print(fullhelper_fullhelp_message)
    print(usage)
    exit(0)
arglen = len(argv)
if arglen < 4 or arglen > 5:
    print(usage)
    exit(0)
# will be true if arglen == 3
cols = 80 # set default value of cols
# otherwise take user input on value of cols
if arglen == 5:
    cols = int(argv[4])
# signal text within sourcecode to locate where fullhelp should be put
signal = argv[3]
print(signal)
# textwrapper class instance, for wrapping the paragraphs in fullhelp msg.
wrapper = textwrap.TextWrapper(width=cols, initial_indent='\"', subsequent_indent='\"',drop_whitespace=False)
# the string to append to the end of each wrapped line
tail = '\\n\"\n'

# using wrapper, build fullhelpstring line by line from fullhelp text file (softwrapped)
fullhelpstr = ''
with open(argv[1], 'r') as fullhelp_text_file:
    for line in fullhelp_text_file:
        fullhelpstr += tail.join(wrapper.wrap(line)) + tail # must add tail to terminate last elt

# also write the string to stdout, as a first pass at courtesy to the user
print('This is the fullhelp string you\'ll print:')
print(fullhelpstr)

# Modify the source and write it to fullhelp-test prefixed file
with open(argv[2],'r') as loos_source_fo, open('fullhelp-test-'+argv[2], 'w') as test_fo:
    found_already = False
    for source_line in loos_source_fo:
        if not found_already and signal in source_line:
            test_fo.write(source_line.replace(signal, fullhelpstr))
            found_already = True
        else:
            test_fo.write(source_line)
