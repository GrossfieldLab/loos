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
LOOS=$(HOME)/Packages/LOOS-1.0

### Which compiler?
CXX=/usr/bin/g++

### Release/debugging options
OPTS=-g -Wall -fno-inline
#OPTS=-O3 -DNDEBUG -Wall

### Use the following for Linux
BOOSTINC=-I/usr/include
BOOSTLIB=
ATLASLIB=-L/usr/lib64/atlas
ATLASINC=-I/usr/include/atlas
LDLIBS=-lloos -lboost_regex -llapack -latlas
LDFLAGS=-L$(LOOS) $(BOOSTLIB) $(ATLASLIB)
CXXFLAGS=$(OPTS) -I. -I$(LOOS) $(BOOSTINC) $(ATLASINC)

### Use the following for MacOS X
#BOOSTINC=-I/usr/local/include
#BOOSTLIB=-L/usr/local/include
#LDLIBS=-lloos -lboost_regex
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




