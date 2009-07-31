##############################################################
# Skeleton makefile for use with LOOS
#
# (c) 2008-2009 Tod D. Romo
#     Grossfield Lab
#     Department of Biochemistry & Biophysics
#     University of Rochester Medical School
#
##############################################################


### Where is LOOS?  This assumes you've installed LOOS somewhere...
### If you haven't installed LOOS, you will also need to change the
### LOOSINC and LOOSLIB appropriately...
LOOS=/opt/loos
LOOSINC=$(LOOS)/include
LOOSLIB=$(LOOS)/lib
#LOOSINC=$(LOOS)
#LOOSLIB=$(LOOS)

### Which compiler?
CXX=/usr/bin/g++

### Platform specific options
# These would normally be set by SCons for you...
#
# The following is required for MacOS X
#POPTS=-DREQUIRES_ULONG

### Release/debugging options
OPTS=-g -Wall -fno-inline $(POPTS)
#OPTS=-O3 -DNDEBUG -Wall $(POPTS)

### Use the following for Linux
BOOSTINC=-I/usr/include
BOOSTLIB=
ATLASLIB=-L/usr/lib64/atlas
ATLASINC=-I/usr/include/atlas
LDLIBS=-lloos -lboost_regex -lboost_program_options -llapack -latlas
LDFLAGS=-L$(LOOSLIB) $(BOOSTLIB) $(ATLASLIB)
CXXFLAGS=$(OPTS) -I. -I$(LOOSINC) $(BOOSTINC) $(ATLASINC)

### Use the following for MacOS X
#BOOSTINC=-I/usr/local/include      # Note: this will probably need to
                                    # change depending on your installation
#BOOSTLIB=-L/usr/local/lib
#LDLIBS=-lloos -lboost_regex -lboost_program_options   # Note: the library names
                                                       # will most likely need to change too...
#LDFLAGS=-L$(LOOSLIB) $(BOOSTLIB)
#CXXFLAGS=$(OPTS) -I. -I$(LOOSINC) $(BOOSTINC) -framework vecLib



###

SRCS=center-pdb.cpp
OBJS=$(SRCS:.cpp=.o)
TARGETS=center-pdb




all: $(TARGETS)

clean:
	-rm $(OBJS) $(TARGETS)



center-pdb: center-pdb.cpp




