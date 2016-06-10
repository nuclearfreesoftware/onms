#NMSAlphaNData Library

message(STATUS "Building NMSAlphaNData")
message(STATUS "Source folder: ${CMAKE_SOURCE_DIR}")

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(NMSAlphaNData)

find_package(Geant4)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# GEANT
find_package(Geant4 REQUIRED vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# SOURCES AND HEADERS of PROJECT
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(NMSAlphaNData ${sources} ${headers})
target_link_libraries(NMSAlphaNData ${Geant4_LIBRARIES})

# add_executable (test_AlphaNData test_AlphaNData.cc)
# target_link_libraries (test_AlphaNData NMSAlphaNData)
# target_link_libraries (test_AlphaNData ${Geant4_LIBRARIES})

# add_executable (test_AlphaNSpectrum test_AlphaNSpectrum.cc)
# target_link_libraries (test_AlphaNSpectrum NMSAlphaNData)
# target_link_libraries (test_AlphaNSpectrum ${Geant4_LIBRARIES})

# add_executable (test_stoppingpower test_stoppingpower.cc)
# target_link_libraries (test_stoppingpower NMSAlphaNData)
# target_link_libraries (test_stoppingpower ${Geant4_LIBRARIES})

install (FILES ${headers} DESTINATION include)
install (TARGETS NMSAlphaNData DESTINATION lib)
# install (TARGETS test_AlphaNData DESTINATION bin )

