#!/usr/bin/env python3
from sys import argv
import textwrap

usage = "path/to/loos/utils/fullhelper.py thefullhelpmessage.txt thetool.cpp '\"the locator string\"' [desired_num_cols]"
fullhelper_fullhelp_message = """
This tool helps you wrap fullhelp messages for loos C++ source so that 
they will be formatted correctly in the whizzywhig loos style, where all the 
needed punctuation and whitespace are specified in source. The objective was to
provide a way for users to write their desired fullhelp message in softwrapping
text editors, fairly common these days, as a plain text file. The correct
wrapping and line-decorating are then done by this script, and the adjusted
text is then written into the source at a spot specified by a locator or signal
string. This is done by assembling the entire updated source code as a list,
writing a backup of the original source code (in case something bad happens),
then overwriting the source file with the newly assembled list of lines of
source code, which SHOULD be the same in every way except the signal string 
will be replaced by the correctly formatted lines of strings that comprise the
fullhelp message.

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
script will find all instances of the signal string and replace them. This 
means you need to include its surrounding quotes (easily done in bash by 
including the signal string inside of doublequotes nested in single quotes, see
usage string). It also means that if you use a signal string you include 
elsewhere, you will get additional copies of your fullhelp message inserted in
those places. This is one reason the backup gets written. The other being that
having the backup makes it easier to go back to editing your full help message
if you decide it needs to be changed after it's been inserted. 

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

# Modify the source storing it as a list of lines.
# Also creates a .backup file in case you wanted one of those...
newsource = []
with open(argv[2],'r') as loos_source_file, open(argv[2]+'.backup', 'w') as backup:
    for source_line in loos_source_file:
        backup.write(source_line)
        newsource.append(source_line.replace(signal, fullhelpstr))
# overwrite the old source code with the new source code
with open(argv[2],'w') as loos_source_file:
    for line in newsource:
        loos_source_file.write(line)