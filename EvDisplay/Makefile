EXE = evDisplay.exe

FASERPSINCLUDEDIR := /home/rubbiaa/faserps/FASERPS/include/
FASERPSSRCDIR := /home/rubbiaa/faserps/FASERPS/src/
# FASERPSLIB := -L/home/rubbiaa/faserps/FASERPS/ -lFaserPS

SRC = evDisplay.C MyMainFrame.cc dict.cxx $(FASERPSSRCDIR)/TcalEvent.cc $(FASERPSSRCDIR)/TPOEvent.cc $(FASERPSSRCDIR)/../G__PS.cxx
OBJ = $(SRC:.C=.o)
OBJ := $(OBJ:.cxx=.o)
OBJ := $(OBJ:.cc=.o)

INCS = $(wildcard *.h)

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(FASERPSINCLUDEDIR)

all : $(EXE) $(EXE2)

dict.cxx : $(INCS) evDisplay.C
	rootcling -f dict.cxx -c LinkDef.h -I $(FASERPSINCLUDEDIR)

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
	rm -f $(OBJS) $(EXEC) dict.cxx dict_rdict.pcm
