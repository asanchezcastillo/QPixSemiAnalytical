cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

#set(CMAKE_C_COMPILER "/usr/local/bin/gcc-12")
#set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-12")

project(Semi)

# Locate the ROOT package and define a number of useful targets and variables.

find_package(ROOT REQUIRED COMPONENTS RIO Net MathCore Hist Physics Core ROOTDataFrame)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

find_package(nlohmann_json 3.2.0 REQUIRED)
add_library(SemiLib SHARED PropagationTimeModel.cc SemiAnalyticalModel.cc)

target_link_libraries(SemiLib PUBLIC ROOT::RIO ROOT::Net ROOT::MathCore ROOT::Hist ROOT::Physics ROOT::Core ROOT::ROOTDataFrame)
target_link_libraries(SemiLib PRIVATE nlohmann_json::nlohmann_json)
# Create the main program using the library.
   
add_executable(OpticalSimulation main.cc PropagationTimeModel.cc SemiAnalyticalModel.cc)
target_link_libraries(OpticalSimulation SemiLib)
