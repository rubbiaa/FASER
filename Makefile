TOPDIR = $(shell pwd)

PATCHFILE = $(TOPDIR)/genfit.patch
FILETOAPPLY = $(TOPDIR)/GenFit/test/vertexingTest/main.cc

.PHONY: clhep rave

clhep: clhep_git
	if [ ! -d CLHEP-install ]; then \
		mkdir -p CLHEP-build; \
		mkdir -p CLHEP-install; \
		cd CLHEP-build && cmake -DCMAKE_INSTALL_PREFIX=../CLHEP-install ../CLHEP; \
		make -j; \
		make install; \
	fi

clhep_git:
	if [ ! -f CLHEP/README.md ]; then \
		git clone https://gitlab.cern.ch/CLHEP/CLHEP.git; \
	fi

rave_tar: 
	if [ ! -d rave-0.6.25 ]; then \
		tar -xvf rave-0.6.25.tar.gz; \
	fi

rave: rave_tar
	if [ ! -d rave-install ]; then \
		mkdir -p rave-install; \
		cd rave-0.6.25 && ./configure --prefix=$(TOPDIR)/rave-install --disable-java \
		--with-clhep=$(TOPDIR)/CLHEP-install; \
		make CXXFLAGS="-g -std=c++11" -j; \
		make install; \
	fi

genfit_git:
	if [ ! -d GenFit ]; then \
		git clone https://github.com/GenFit/GenFit.git; \
	fi

genfit: genfit_git
	if [ ! -d GenFit-build ]; then \
		mkdir -p GenFit-build; \
		mkdir -p GenFit-install; \
		cd GenFit-build && cmake \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCMAKE_INSTALL_PREFIX=$(TOPDIR)/GenFit-install \
		-DRave_CFLAGS="-DRaveDllExport= -DWITH_FLAVORTAGGING -DWITH_KINEMATICS" \
		-DRave_INCLUDE_DIRS=$(TOPDIR)/rave-install/include/ \
		-DRave_LDFLAGS="-Wl,-rpath-link,$(TOPDIR)/rave-install/lib/ -L$(TOPDIR)/rave-install/lib/ -lRaveBase -L$(TOPDIR)/CLHEP-install/lib/ -lCLHEP" \
		../GenFit; \
		make CXXFLAFS="-g" -j; \
		sh CMakeFiles/gtests.dir/link.txt; \
		make -j; \
		make install; \
	fi
	# if make fails, run sh CMakeFiles/gtests.dir/link.txt and then make again

.PHONY clean:
	rm -rf CLHEP-build CLHEP-install
	rm -rf CLHEP
	rm -rf rave-0.6.25 rave-install
	rm -rf GenFit-build GenFit-install
	rm -rf GenFit

