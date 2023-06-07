cmake_minimum_required(VERSION 3.7)
project(MicMac)

set(CMAKE_CXX_STANDARD 14)

# Compilation options
option(DEBUG_MODE "Debug mode (with asserts and such) " OFF) #OFF by default
set(ALWAYSON_CXX_FLAGS "--std=c++14")
if(DEBUG_MODE)
    set(CMAKE_CXX_FLAGS "-O0 -g ${ALWAYSON_CXX_FLAGS}") # debug mode
    message("-- INFO: Compiling in debug mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
else()
    set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG ${ALWAYSON_CXX_FLAGS}") # release mode
    message("-- INFO: Compiling in release mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
endif(DEBUG_MODE)
unset(DEBUG_MODE CACHE)

include_directories("src")
include_directories("src/lib")
include_directories("utils")

add_executable(neutral src/neutral.cpp)
add_executable(selection src/selection.cpp)