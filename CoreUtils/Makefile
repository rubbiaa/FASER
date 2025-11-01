include ../Makefile.common

.PHONY : dict
dict : LinkDef.h $(INCS) $(SRCCORE)
	rootcling -f dict.cxx -c LinkDef.h -I$(COREUTILSDIR) -I$(GENFITINCDIR)
	g++ $(ROOTCFLAGS) -fPIC -shared -o libCoreUtilsDict.so dict.cxx *.cc $(SRCCORE:.cc=.o) -L$(ROOTSYS)/lib `root-config --glibs` $(GENFITLIBS)
