g++ $(root-config --cflags) -O2 CombineFluxes.C -o CombineFluxes $(root-config --glibs)

g++ -o countInteractions countInteractions.cc `root-config --cflags --libs`

g++ -o plotInteractionVertices plotInteractionVertices.cc `root-config --cflags --libs`

g++ -o etaCoverage etaCoverage.cc `root-config --cflags --libs`


rootcling -f dict.cxx -c LinkDef.h -I../CoreUtils -I${PYTHIA8}/include

g++ -g `root-config --cflags` \
    -I../CoreUtils -I${PYTHIA8}/include -D_INCLUDE_PYTHIA_ \
    -c dict.cxx -o dict.o

g++ -g `root-config --cflags` \
    -I../CoreUtils -I${PYTHIA8}/include -D_INCLUDE_PYTHIA_ \
    countinteractions_tpoevent.cc ../CoreUtils/TPOEvent.cc dict.o \
    `root-config --glibs` -lEG -lEGPythia8 -lGeom -lPhysics \
    -L${PYTHIA8}/lib -lpythia8 -Wl,-rpath,${PYTHIA8}/lib \
    -o countinteractions_tpoevent.exe