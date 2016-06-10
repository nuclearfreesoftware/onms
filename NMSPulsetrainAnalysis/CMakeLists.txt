cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(NMSPulsetrainAnalysis)

find_package( PkgConfig REQUIRED)
pkg_check_modules( gsl REQUIRED gsl )

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(NMSPulsetrainAnalysis ${sources} ${headers})
target_link_libraries( NMSPulsetrainAnalysis ${gsl_LIBRARIES} )

add_executable (nmspa_example example/nmspa_example.cc)
target_link_libraries (nmspa_example NMSPulsetrainAnalysis)

install (FILES ${headers} DESTINATION include)
install (TARGETS NMSPulsetrainAnalysis DESTINATION lib)
install (TARGETS nmspa_example DESTINATION example )
install (FILES example/cf252_example_pulsetrain_list DESTINATION example)
