EXE = evDisplay.exe

FASERPSINCLUDEDIR := /home/rubbiaa/faserps/FASERPS/include/
FASERPSSRCDIR := /home/rubbiaa/faserps/FASERPS/src/
# FASERPSLIB := -L/home/rubbiaa/faserps/FASERPS/ -lFaserPS

SRC = evDisplay.C $(FASERPSSRCDIR)/TcalEvent.cc $(FASERPSSRCDIR)/../G__PS.cxx
OBJ = $(SRC:.C=.o)

ROOTCFLAGS = -g $(shell root-config --cflags)

ROOTCFLAGS += -I$(FASERPSINCLUDEDIR)

all : $(EXE) $(EXE2)

$(EXE) : $(OBJ)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -lGeom $(FASERPSLIB) -o $(EXE)

# Compile source files to object files
%.o: %.C
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cxx
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
