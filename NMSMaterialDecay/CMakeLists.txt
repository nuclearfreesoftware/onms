#NMSMaterialDecay Library

message(STATUS "Building NMSMaterialDecay")
message(STATUS "Source folder: ${CMAKE_SOURCE_DIR}")

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(NMSMaterialDecay)

include_directories(/geant/fission)

find_package(Geant4)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/../NMSAlphaNData/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(NMSMaterialDecay ${sources} ${headers})
target_link_libraries(NMSMaterialDecay NMSAlphaNData)

install (FILES ${headers} DESTINATION include)
install (TARGETS NMSMaterialDecay DESTINATION lib)

