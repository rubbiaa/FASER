EXE = batchreco.exe
EXE2 = dumphits.exe

include ../Makefile.common

SRC = BatchReco.cc
SRC2 = DumpHits.cc

SRCCORE = $(wildcard $(COREUTILSDIR)/*.cc)
OBJCORE := $(SRCCORE:.cc=.o)

SRC += $(SRCCORE)
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)
OBJ += dict.o

SRC2 += $(SRCCORE)
OBJ2 = $(SRC2:.C=.o)
OBJ2 := $(OBJ2:.cxx=.o)
OBJ2 := $(OBJ2:.cc=.o)
OBJ2 += dict.o

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)
INCS += $(wildcard $(COREUTILSDIR)/*.hh)

# ROOTCFLAGS += -fsanitize=address 

all : $(EXE) $(EXE2)

dict.cxx : $(INCS) $(SRC)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR) -I$(GENFITINCDIR)

$(EXE) : $(OBJ) Makefile
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom $(GENFITLIBS) -o $@
	g++ -shared -o libTPORec.so dict.o $(OBJCORE) -L$(ROOTSYS)/lib `root-config --glibs`

$(EXE2) : $(OBJ2) Makefile
	g++ $(ROOTCFLAGS) $(OBJ2) `root-config --glibs` -lEG -lGeom $(GENFITLIBS) -o $@
	g++ -shared -o libTPORec.so dict.o $(OBJCORE) -L$(ROOTSYS)/lib `root-config --glibs`

# Compile source files to object files
%.o: %.C $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cxx $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJ) $(EXE) $(EXE2) $(OBJ2) dict.cxx dict_rdict.pcm libTPORec.so

.PHONY : all clean

