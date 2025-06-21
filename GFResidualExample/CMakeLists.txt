cmake_minimum_required(VERSION 2.8)
project(GFResidualExample LANGUAGES CXX)


# Find ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist)
include(${ROOT_USE_FILE})

# Use environment variables for GenFit2
# It is expected you have something like
#   export GENFIT=...$USER/genfit/lib/
#   export GENFIT_INCDIR=${GENFIT}/../include/
#   export LD_LIBRARY_PATH=${GENFIT}:${LD_LIBRARY_PATH}
# in your ~/.bashrc
find_path(GENFIT2_INCLUDE_DIR
    NAMES AbsTrackRep.h
    PATHS $ENV{GENFIT_INCDIR}
    DOC "GenFit2 include directory"
)

find_library(GENFIT2_LIBRARY
    NAMES genfit2
    PATHS $ENV{GENFIT}
    DOC "GenFit2 library"
)
if(NOT GENFIT2_INCLUDE_DIR OR NOT GENFIT2_LIBRARY)
    message(FATAL_ERROR "Could not find GenFit2. Please ensure GENFIT and GENFIT_INCDIR are set in your environment.")
endif()

message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")
message(STATUS "GenFit include dir: ${GENFIT2_INCLUDE_DIR}")
message(STATUS "GenFit library: ${GENFIT2_LIBRARY}")


include_directories(${GENFIT2_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

add_executable(GFResidualExample minimal_genfit_residual_example.cpp)
target_link_libraries(GFResidualExample ${GENFIT2_LIBRARY} ${ROOT_LIBRARIES})

