TOPDIR = $(shell pwd)

PYTHIA8_DIR = $(TOPDIR)/pythia8312
PYTHIA8_INCLUDE_DIR = $(PYTHIA8_DIR)/include
PYTHIA8_LIBRARY = $(PYTHIA8_DIR)/lib/libpythia8.a

OPENGL = /usr/include/GL
OPEN_gl_LIBRARY=/usr/

PATCHFILE = $(TOPDIR)/genfit.patch
FILETOAPPLY = $(TOPDIR)/GenFit/test/vertexingTest/main.cc

.PHONY: root pythia8 clhep rave genfit

root_src:
	if [ ! -d root-6.32.02 ]; then \
		wget https://cernbox.cern.ch/s/5LgzZYg141mQ03c/download -O root_v6.32.02.source.tar.gz; \
		tar xvfz root_v6.32.02.source.tar.gz;\
	fi

root_cmake:
		cd root-build && cmake \
		-DCMAKE_INSTALL_PREFIX=$(TOPDIR)/root-install \
		-DPYTHIA8_INCLUDE_DIR=$(PYTHIA8_INCLUDE_DIR) \
		-DPYTHIA8_LIBRARY=$(PYTHIA8_LIBRARY) \
		-Dpythia8=ON -Dopengl=ON -Dbuiltin_openui5=OFF \
		../root-6.32.02/

root_make:
	cd root-build && make -j4

root_install:
	cd root-build && make install

root:
	if [ ! -d root-build ]; then \
		mkdir -p root-build; \
		mkdir -p root-install; \
		cd root-build && cmake \
		-DCMAKE_INSTALL_PREFIX=$(TOPDIR)/root-install \
		-DPYTHIA8_INCLUDE_DIR=$(PYTHIA8_INCLUDE_DIR) \
		-DPYTHIA8_LIBRARY=$(PYTHIA8_LIBRARY) \
		-Dpythia8=ON -Dopengl=ON -Dbuiltin_openui5=OFF \
		../root-6.32.02/; \
	fi

pythia8_tar:
	if [ ! -d pythia8312 ]; then \
		wget https://cernbox.cern.ch/s/lOhu8P3H0beVnb0/download -O  pythia8312.tgz; \
		tar xvfz pythia8312.tgz; \
	fi

pythia8: pythia8_tar
	if [ ! -d pythia8312/lib ]; then \
		cd pythia8312; \
		./configure --prefix=$(TOPDIR)/pythia-install; \
		make; \
	fi

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
		cd Genfit; \
		patch -p0 -u -i ../genfit.patch; \
	fi

genfit: genfit_git
	if [ ! -d GenFit-build ]; then \
		mkdir -p GenFit-build; \
		mkdir -p GenFit-install; \
		cd GenFit-build && cmake \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCMAKE_INSTALL_PREFIX=$(TOPDIR)/GenFit-install \
		-DGTEST_LIBRARY=$(TOPDIR)/googletest-install/lib64/libgtest.a -DGTEST_INCLUDE_DIR=$(TOPDIR)/googletest-install/include \
		-DGTEST_MAIN_LIBRARY=$(TOPDIR)/googletest-install/lib64/libgtest_main.a \
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

googletest_git:
	if [ ! -d googletest ]; then \
		git clone https://github.com/google/googletest.git; \
	fi

googletest: googletest_git
	if [ ! -d googletest-build ]; then \
		mkdir -p googletest-build; \
		mkdir -p googletest-install; \
		cd googletest-build && cmake \
		-DCMAKE_INSTALL_PREFIX=$(TOPDIR)/googletest-install \
		../googletest; \
	fi


.PHONY clean:
	rm -rf pythia8312 pythia8312.tgz
	rm -rf CLHEP-build CLHEP-install
	rm -rf CLHEP
	rm -rf rave-0.6.25 rave-install
	rm -rf GenFit-build GenFit-install
	rm -rf GenFit
	rm -rf googletest googletest-install googletest-build


