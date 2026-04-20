EXE = ConvertGENIE.exe

COREUTILSDIR := ../CoreUtils

SRC = ConvertGENIE.cc
OBJ = ConvertGENIE.o TPOEvent.o dict.o

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

all : $(EXE) $(EXE_ECAL) $(EXE_AHCAL)

dict.cxx : $(INCS)
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

# Compile TPOEvent from CoreUtils with our flags
TPOEvent.o: $(COREUTILSDIR)/TPOEvent.cc $(COREUTILSDIR)/TPOEvent.hh
	$(CXX) $(ROOTCFLAGS) -c $(COREUTILSDIR)/TPOEvent.cc -o TPOEvent.o

# Clean target
clean:
	rm -f $(OBJ) $(EXE) dict.cxx dict_rdict.pcm ConvertGENIE.o TPOEvent.o

.PHONY : all clean

