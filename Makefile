EXE = t.exe

SRC = t.C faserntuplib.o
OBJ = $(SRC:.C=.o)

ROOTCFLAGS = -g $(shell root-config --cflags)

all : $(EXE)

$(EXE) : $(OBJ)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -o $(EXE)

# Compile source files to object files
%.o: %.C
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
