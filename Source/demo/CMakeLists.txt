# Add executable demo file
ADD_EXECUTABLE(demo demo.cpp)
# Link dynamic PSINS library to demo
TARGET_LINK_LIBRARIES(demo PSINS_static)

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    # Add a demo test
    # Run executable binary file "demo" under direcotry "CMAKE"/demo
    ADD_CUSTOM_TARGET(test demo WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/demo)
endif()
