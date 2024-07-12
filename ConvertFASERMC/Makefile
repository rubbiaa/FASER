EXE = Convert.exe

COREUTILSDIR := ../CoreUtils

SRC = Convert.cc $(COREUTILSDIR)/TPOEvent.cc
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)
OBJ += dict.o

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)
INCS += $(COREUTILSDIR)/TPOEvent.hh

# Pythia8 library
PYTHIA8INC = $(PYTHIA8)/include
PYTHIA8LIB = $(PYTHIA8)/lib

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(COREUTILSDIR) -I$(PYTHIA8INC) -D_INCLUDE_PYTHIA_

LDLIBS = `root-config --glibs` -lEG -lEGPythia8 -lGeom -lPhysics
LDLIBS += -L$(PYTHIA8LIB) -lpythia8
LDLIBS += -Wl,-rpath,$(PYTHIA8LIB)

all : $(EXE)

dict.cxx : $(INCS) $(SRC)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR) -I$(PYTHIA8INC)

$(EXE) : $(OBJ) Makefile
	g++ $(ROOTCFLAGS) $(OBJ) $(LDLIBS) -o $(EXE)

# Compile source files to object files
%.o: %.C
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cxx
	$(CXX) $(ROOTCFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJ) $(EXE) dict.cxx dict_rdict.pcm

.PHONY : all clean

