EXE = t.exe

all : $(EXE)

$(EXE) : t.C
	g++ `root-config --cflags` t.C `root-config --glibs` -lEG -o $(EXE)	
