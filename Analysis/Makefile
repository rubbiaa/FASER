EXE = AnalyReco.exe
EXE2 = plots.exe

COREUTILSDIR := ../CoreUtils

SRC = AnalyReco.cc
SRC2 = plots.C

SRCCORE = $(wildcard $(COREUTILSDIR)/*.cc)
OBJCORE := $(SRCCORE:.cc=.o)

SRC += $(SRCCORE)
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)
OBJ += dict.o

OBJ2 = $(SRC2:.C=.o)

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)
INCS += $(wildcard $(COREUTILSDIR)/*.hh)

ROOTCFLAGS = -fPIC -g $(shell root-config --cflags) 

ROOTCFLAGS += -I$(COREUTILSDIR)

all : $(EXE) $(EXE2)

dict.cxx : $(INCS) $(SRC)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR)

$(EXE) : $(OBJ) Makefile
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom -o $@
	g++ -shared -o libTPORec.so dict.o $(OBJCORE) -L$(ROOTSYS)/lib `root-config --glibs`

$(EXE2): $(OBJ2) Makefile
	g++ -o $(EXE2) $(OBJ2) `root-config --glibs` -lRooFit -lRooFitCore

# Compile source files to object files
%.o: %.C $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cxx $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJ) $(EXE) dict.cxx dict_rdict.pcm libTPORec.so

.PHONY : all clean

