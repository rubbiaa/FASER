EXE = t.exe
EXE2 = s.exe

SRC = t.C faserntuplib.o
OBJ = $(SRC:.C=.o)

SRC2 = s.C TKinFitter.cc TFitParticleESpher.cc TAbsFitParticle.cc TFitConstraintEp.cc TAbsFitConstraint.cc
OBJ2 = s.o TKinFitter.o TFitParticleESpher.o TAbsFitParticle.o TFitConstraintEp.o TAbsFitConstraint.o

ROOTCFLAGS = -g $(shell root-config --cflags)

all : $(EXE) $(EXE2)

$(EXE) : $(OBJ)
	g++ $(ROOTCFLAGS) $(OBJ) `root-config --glibs` -lEG -o $(EXE)

$(EXE2) : $(OBJ2)
	g++ $(ROOTCFLAGS) $(OBJ2) `root-config --glibs` -lEG -o $(EXE2)

# Compile source files to object files
%.o: %.C
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
%.o: %.cc
	$(CXX) $(ROOTCFLAGS) -c $< -o $@
