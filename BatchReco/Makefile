EXE = batchreco.exe

FASERPSINCLUDEDIR := /home/rubbiaa/faserps/FASERPS/include/
FASERPSSRCDIR := /home/rubbiaa/faserps/FASERPS/src/
# FASERPSLIB := -L/home/rubbiaa/faserps/FASERPS/ -lFaserPS

COREUTILSDIR := ../CoreUtils

SRC = BatchReco.cc $(COREUTILSDIR)/TPORecoEvent.cc dict.cxx $(FASERPSSRCDIR)/TcalEvent.cc $(COREUTILSDIR)/TPOEvent.cc
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(FASERPSINCLUDEDIR) -I$(COREUTILSDIR)

all : $(EXE)

dict.cxx : $(INCS) BatchReco.cc
	rootcling -f dict.cxx -c LinkDef.h -I $(FASERPSINCLUDEDIR) -I$(COREUTILSDIR)

$(EXE) : $(OBJ)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom $(FASERPSLIB) -o $(EXE)

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

