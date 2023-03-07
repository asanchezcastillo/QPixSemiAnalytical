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
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
add_library(SemiLib SHARED PropagationTimeModel.cc SemiAnalyticalModel.cc)

target_link_libraries(SemiLib PUBLIC ROOT::RIO ROOT::Net ROOT::MathCore ROOT::Hist ROOT::Physics ROOT::Core ROOT::ROOTDataFrame)
target_link_libraries(SemiLib PRIVATE nlohmann_json::nlohmann_json)
# Create the main program using the library.
   
add_executable(main main.cc PropagationTimeModel.cc SemiAnalyticalModel.cc)
target_link_libraries(main SemiLib)
