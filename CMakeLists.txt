cmake_minimum_required( VERSION 2.8 )
project ( FermiOwnHMC )
include(CTest)

set( FermiOwnHMC_VERSION_MAJOR 0 )
set( FermiOwnHMC_VERSION_MINOR 1 )

set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11 -lm -Wall -g -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -std=c++11 -lm -Wall -O0 -g")

include_directories($ENV{EIGEN3_INCLUDE_DIR})
include_directories("src/Actions")
include_directories("src/UpdateAlgorithms")
include_directories("src/UpdateAlgorithms/HMC_Support")
include_directories("FermiOwnCommon/src/Math")
include_directories("FermiOwnCommon/src/Data")

add_subdirectory("FermiOwnCommon")
add_subdirectory("src")
add_subdirectory("test")
#add_subdirectory("apps")