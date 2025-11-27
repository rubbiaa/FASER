cmake_minimum_required(VERSION 3.12)
project(FaserCalDisplay)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# ======== ROOT Configuration ========
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS
    Core Hist Matrix Physics GenVector Eve EG Geom Gpad RGL Rint Gui GuiHtml Thread Gdml)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(-Wno-deprecated-declarations)

# ======== GenFit Configuration ========
set(GENFITINCDIR /Users/ukose/sw/kits/GenFit-install/include)
set(GENFITLIBDIR /Users/ukose/sw/kits/GenFit-install/lib)

include_directories(${GENFITINCDIR})
link_directories(${GENFITLIBDIR})

# ======== Sources ========
set(SOURCES
    FaserCalDisplay.cxx
    ../CoreUtils/TPOEvent.cc
    ../CoreUtils/TPORecoEvent.cc
    ../CoreUtils/TcalEvent.cc
    ../CoreUtils/TTKTrack.cc
    ../CoreUtils/TPSCluster.cc
    ../CoreUtils/DBScan.cc
    ../CoreUtils/TPSTrack.cc
    ../CoreUtils/GenMagneticField.cc
    ../CoreUtils/TMuTrack.cc
    )

# ======== ROOT Dictionary ========
root_generate_dictionary(G__FaserCalDisplay FaserCalDisplay.h LINKDEF disp_linkdef.h)

# ======== Shared Library ========
add_library(FaserCalDisplay SHARED ${SOURCES} G__FaserCalDisplay.cxx)
target_include_directories(FaserCalDisplay PRIVATE
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ../CoreUtils
    ${GENFITINCDIR}
)
target_link_libraries(FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

# ======== Executables ========
add_executable(FaserCalDisplayApp main.cxx)
target_include_directories(FaserCalDisplayApp PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(FaserCalDisplayApp FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

add_executable(FaserCalAnalyzer FaserCalAnalyzer.cxx)
target_include_directories(FaserCalAnalyzer PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(FaserCalAnalyzer FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

add_executable(FaserCalAnalyzer2nd FaserCalAnalyzer2nd.cxx)
target_include_directories(FaserCalAnalyzer2nd PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(FaserCalAnalyzer2nd FaserCalDisplay ${ROOT_LIBRARIES} genfit2)


add_executable(DumpMuonHits DumpMuonHits.cxx)
target_include_directories(DumpMuonHits PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(DumpMuonHits FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

add_executable(DumpEventList DumpEventList.cxx)
target_include_directories(DumpEventList PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(DumpEventList FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

#add_executable(FaserCalMuonSpectrometer FaserCalMuonSpectrometer.cxx)
#target_include_directories(FaserCalMuonSpectrometer PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
#target_link_libraries(FaserCalMuonSpectrometer FaserCalDisplay ${ROOT_LIBRARIES} genfit2)


add_executable(DumpMuonSpectrometer DumpMuonSpectrometer.cxx)
target_include_directories(DumpMuonSpectrometer PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_link_libraries(DumpMuonSpectrometer FaserCalDisplay ${ROOT_LIBRARIES} genfit2)

# ======== RPATH for macOS Runtime Discovery ========
set(CMAKE_INSTALL_RPATH "@loader_path/../lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
list(APPEND CMAKE_INSTALL_RPATH ${GENFITLIBDIR} $ENV{ROOTSYS}/lib)
##in linux
#set(CMAKE_INSTALL_RPATH "${ROOT_LIBRARY_DIR};${GENFITLIBDIR};$ENV{ROOTSYS}/lib")
#set(CMAKE_BUILD_RPATH "${ROOT_LIBRARY_DIR};${GENFITLIBDIR};$ENV{ROOTSYS}/lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# Debug
# Debug output
message(STATUS "ROOT Library Directory: ${ROOT_LIBRARY_DIR}")
message(STATUS "GenFit Library Directory: ${GENFITLIBDIR}")
message(STATUS "GenFit Library: ${GENFIT_LIB}")

