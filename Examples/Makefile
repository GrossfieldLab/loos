##############################################################
# Skeleton makefile for use with LOOS
#
# (c) 2008 Tod D. Romo
#     Grossfield Lab
#     Department of Biochemistry & Biophysics
#     University of Rochester Medical School
#
##############################################################


### Where is LOOS?

# An absolute path
#LOOS=$(HOME)/Packages/LOOS-1.0

# Or a relative path...
LOOS=..


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
LDFLAGS=-L$(LOOS) $(BOOSTLIB) $(ATLASLIB)
CXXFLAGS=$(OPTS) -I. -I$(LOOS) $(BOOSTINC) $(ATLASINC)

### Use the following for MacOS X
#BOOSTINC=-I/usr/local/include
#BOOSTLIB=-L/usr/local/include
#LDLIBS=-lloos -lboost_regex -lboost_program_options
#LDFLAGS=-L$(LOOS) $(BOOSTLIB)
#CXXFLAGS=$(OPTS) -I. -I$(LOOS) $(BOOSTINC) -framework vecLib



###

SRCS=center-pdb.cpp
OBJS=$(SRCS:.cpp=.o)
TARGETS=center-pdb




all: $(TARGETS)

clean:
	-rm $(OBJS) $(TARGETS)



center-pdb: center-pdb.cpp




