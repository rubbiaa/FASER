EXE = evDisplay.exe

include ../Makefile.common

PANDORA_DIR = /data/sw/products/PandoraPFA/install
FASTJET_DIR = /data/sw/products/fastjet-install
FASTJET_INCLUDE = -I$(FASTJET_DIR)/include
FASTJET_LIB = -L$(FASTJET_DIR)/lib -lfastjet

ROOTCFLAGS +=  -fsanitize=address $(FASTJET_INCLUDE)

SRC = evDisplay.C MyMainFrame.cc 

##SRCCORE = $(wildcard $(COREUTILSDIR)/*.cc)
SRCCORE = $(COREUTILSDIR)/TPOEvent.cc $(COREUTILSDIR)/TPORecoEvent.cc $(COREUTILSDIR)/TcalEvent.cc $(COREUTILSDIR)/TTKTrack.cc $(COREUTILSDIR)/DBScan.cc $(COREUTILSDIR)/TPSTrack.cc $(COREUTILSDIR)/TPSCluster.cc 
OBJCORE := $(SRCCORE:.cc=.o)

SRC += $(SRCCORE)

OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)
OBJ += dict.o

INCS = $(wildcard *.h)
INCS += $(wildcard *.hh)
INCS += $(wildcard $(COREUTILSDIR)/*.hh)

all : $(EXE)

dict.cxx : $(INCS) $(SRC)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR) -I$(GENFITINCDIR) $(FASTJET_INCLUDE)

$(EXE) : $(OBJ) Makefile
#	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs --evelibs` -lEG -lGeom -o $(EXE)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom $(GENFITLIBS) $(FASTJET_LIB) -o $(EXE)

# Compile source files to object files
%.o: %.C $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cxx $(INCS)
	$(CXX) $(ROOTCFLAGS) -c $< -o $@

# Clean target
clean:
	rm -f $(OBJ) $(EXE) dict.cxx dict_rdict.pcm

.PHONY : all clean

