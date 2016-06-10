cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(nms)


message(STATUS "Building NMS (main cmake folder)")
message(STATUS "Source folder: ${CMAKE_SOURCE_DIR}")

# FISSION (in source tree)
add_subdirectory(fission)
include_directories(${PROJECT_SOURCE_DIR}/fission/include)

# NMS libraries in project tree
add_subdirectory(NMSAlphaNData)
include_directories(${PROJECT_SOURCE_DIR}/NMSAlphaNData/include)
add_subdirectory(NMSMaterialDecay)
include_directories(${PROJECT_SOURCE_DIR}/NMSMaterialDecay/include)
add_subdirectory(NMSPulsetrainAnalysis)
include_directories(${PROJECT_SOURCE_DIR}/NMSPulsetrainAnalysis/include)


# GEANT
find_package(Geant4 REQUIRED vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# SOURCES AND HEADERS of PROJECT
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# EXECUTABLE
add_executable(nms nms.cc ${sources} ${headers})

###### LINK TO LIBRARIES

#----- Geant
target_link_libraries(nms ${Geant4_LIBRARIES})

#----- NMSAlphaNData
target_link_libraries(nms NMSAlphaNData)

#----- NMSMaterialDecay
target_link_libraries(nms NMSMaterialDecay)

#----- NMSPulsetrainAnalysis
target_link_libraries(nms NMSPulsetrainAnalysis)

#----- FISSION Lib from LLNL (FIXME: add cmake for compilation of sources later)
target_link_libraries(nms Fission)


# There was a good reason for linking NMSMaterialDecay before Fission!

file(GLOB ALL_MAC *.mac)


foreach(_macfile ${ALL_MAC})
  message(STATUS "Copy .mac script '${_macfile}'")
  configure_file(${_macfile} ${PROJECT_BINARY_DIR} COPYONLY)
endforeach()

# copy GDML files, schema files and examples
file(COPY ${PROJECT_SOURCE_DIR}/gdml DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${PROJECT_BINARY_DIR})
