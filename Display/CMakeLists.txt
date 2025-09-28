cmake_minimum_required(VERSION 3.12)
project(FaserCalDisplay)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# Find ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS GenVector Eve EG Geom Gpad RGL Rint Gui GuiHtml Thread Gdml)

# Include ROOT directories
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

add_definitions(-Wno-deprecated-declarations)

# Include GenFit directories
set(GENFITINCDIR ../GenFit-install/include)
include_directories(${GENFITINCDIR})
set(GENFITLIBDIR ../GenFit-install/lib)
link_directories(${GENFITLIBDIR})

# Define the source files
set(SOURCES
    FaserCalDisplay.cxx
    ../CoreUtils/TPOEvent.cc
    ../CoreUtils/TPORecoEvent.cc
    ../CoreUtils/TcalEvent.cc
    ../CoreUtils/TTKTrack.cc
    ../CoreUtils/TPSCluster.cc
    ../CoreUtils/DBScan.cc
    ../CoreUtils/TPSTrack.cc
)

# Generate dictionary
root_generate_dictionary(G__FaserCalDisplay FaserCalDisplay.h LINKDEF disp_linkdef.h)

# Create shared library for FaserCalDisplay
add_library(FaserCalDisplay SHARED ${SOURCES} G__FaserCalDisplay.cxx)

# Diplay app
# Add executable
add_executable(FaserCalDisplayApp main.cxx)

# Include directories
target_include_directories(FaserCalDisplay PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})
target_include_directories(FaserCalDisplayApp PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})

# Link libraries
find_library(GENFIT_LIB libgenfit2.so PATHS ${GENFITLIBDIR})
if (NOT GENFIT_LIB)
    message(FATAL_ERROR "GenFit library not found!")
endif()

target_link_libraries(FaserCalDisplayApp FaserCalDisplay stdc++fs ${ROOT_LIBRARIES} ${GENFIT_LIB})

########################
## Add analyzer executable
add_executable(FaserCalAnalyzer FaserCalAnalyzer.cxx)

target_include_directories(FaserCalAnalyzer PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})

target_link_libraries(FaserCalAnalyzer
    FaserCalDisplay
    stdc++fs
    ${ROOT_LIBRARIES}
    ${GENFIT_LIB}
)

########################
## Add muonchecker executable
add_executable(DumpMuonHits DumpMuonHits.cxx)

target_include_directories(DumpMuonHits PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})

target_link_libraries(DumpMuonHits
    FaserCalDisplay
    stdc++fs
    ${ROOT_LIBRARIES}
    ${GENFIT_LIB}
)
########################
## Add EventList executable
add_executable(DumpEventList DumpEventList.cxx)

target_include_directories(DumpEventList PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ../CoreUtils ${GENFITINCDIR})

target_link_libraries(DumpEventList
    FaserCalDisplay
    stdc++fs
    ${ROOT_LIBRARIES}
    ${GENFIT_LIB}
)



# Add RPATH for runtime library discovery
set(CMAKE_INSTALL_RPATH "${ROOT_LIBRARY_DIR};${GENFITLIBDIR};$ENV{ROOTSYS}/lib")
set(CMAKE_BUILD_RPATH "${ROOT_LIBRARY_DIR};${GENFITLIBDIR};$ENV{ROOTSYS}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Debug output
message(STATUS "ROOT Library Directory: ${ROOT_LIBRARY_DIR}")
message(STATUS "GenFit Library Directory: ${GENFITLIBDIR}")
message(STATUS "GenFit Library: ${GENFIT_LIB}")
