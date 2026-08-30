cmake_minimum_required(VERSION 3.16)
project(stopping)

# Find Geant4 package
find_package(Geant4 REQUIRED)

# Include current directory for headers
include_directories(${PROJECT_SOURCE_DIR})

# Gather source and header files
file(GLOB SOURCES "*.cc")
file(GLOB HEADERS "*.hh")

# Create executable
add_executable(stopping ${SOURCES} ${HEADERS})

# Link Geant4 libraries
target_link_libraries(stopping ${Geant4_LIBRARIES})
