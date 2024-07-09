EXE = batchreco.exe

COREUTILSDIR := ../CoreUtils

SRC = BatchReco.cc $(COREUTILSDIR)/TPORecoEvent.cc dict.cxx $(COREUTILSDIR)/TcalEvent.cc $(COREUTILSDIR)/TPOEvent.cc
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(COREUTILSDIR)

all : $(EXE)

dict.cxx : $(INCS) BatchReco.cc
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

