# Preferable to avoid warnings during execution
cmake_minimum_required(VERSION 3.9)
project("Bilateral Teleoperation Ariadne")

# Set to search the root source directory for the FindAriadne.cmake module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Set binary directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set target name
set(TARGET_NAME main)

# Set to look for includes in usr/local directory, required for GMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/include -g")

# Set the C++ standard to C++17, as used by Ariadne
set(CMAKE_CXX_STANDARD 17)

# Find Ariadne and Yaml packages
find_package(Ariadne REQUIRED)
find_package(YAML REQUIRED)

# Add the Ariadne headers and Yaml headers to the includes
include_directories(SYSTEM ${ARIADNE_INCLUDE_DIRS} ${YAML_INCLUDES})

# Set the executable along with the required source files
add_executable(${TARGET_NAME} main.cpp)

# Instruct to link against the ariadne and yaml libraries
target_link_libraries(${TARGET_NAME} ${ARIADNE_LIBRARIES} ${YAML_LIBRARIES})

# Create plots directory after building the target
add_custom_command(TARGET ${TARGET_NAME} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/build/bin/plots")