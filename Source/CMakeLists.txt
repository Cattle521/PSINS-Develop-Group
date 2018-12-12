#Templates are created by Cattle 2018.12.2
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(PSINS)

#SET(CMAKE_C_COMPILER "/.....-gcc ")
#SET(CMAKE_CXX_COMPILER "/.....-g++ ")

# Compiler setting for gcc
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPSINS_RMEMORY")
    # Run cmake -DCMAKE_BUILD_TYPE=Debug ../
    # Set g++ debug flags
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")
    # Run cmake -DCMAKE_BUILD_TYPE=Release ../
    # Set g++ release flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
endif()

#find_package(Threads REQUIRED)

# Set exectable binary file directory
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add header files directory
INCLUDE_DIRECTORIES(Algorithm)

# Add subdirectory contain "CMakesLists.txt"
ADD_SUBDIRECTORY(demo)
ADD_SUBDIRECTORY(Algorithm)
