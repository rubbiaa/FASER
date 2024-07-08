EXE = Convert.exe

COREUTILSDIR := ../CoreUtils

SRC = Convert.cc dict.cxx $(COREUTILSDIR)/TPOEvent.cc
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)
INCS += $(COREUTILSDIR)/TPOEvent.hh

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(COREUTILSDIR)

all : $(EXE)

dict.cxx : $(INCS) $(SRC)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR)

$(EXE) : $(OBJ)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom -o $(EXE)

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

