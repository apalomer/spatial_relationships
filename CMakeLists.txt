#####################################################################################
#                                                                                   #
#                      CMake Configuration                                          #
#                                                                                   #
#####################################################################################

# Minimum version required
cmake_minimum_required(VERSION 2.8)

# Project Name
project(spatial_relationship)

# Set bin output
set(CMAKE_BIN "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BIN})

# c++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")

# Ceres
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIR})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
add_definitions(${EIGEN_DEFINITIONS})

#Doxygen
option(BUILD_DOC "Build documentation" ON)
if(BUILD_DOC)

    # Find package
    find_package(Doxygen REQUIRED)

    # Doxygen options
    option(DOXYGEN_SHOW_FILES "Display files in doxygen" ON)
    option(DOXYGEN_CALL_GRAPH "Display call graph for each function" ON)
    option(DOXYGEN_CALLER_GRAPH "Display caller graph for each function" ON)
    if(DOXYGEN_SHOW_FILES)
        set(DOXYGEN_SHOW_FILES YES)
    else()
        set(DOXYGEN_SHOW_FILES NO)
    endif()
    if(DOXYGEN_CALL_GRAPH)
        set(DOXYGEN_CALL_GRAPH YES)
    else()
        set(DOXYGEN_CALL_GRAPH NO)
    endif()
    if(DOXYGEN_CALLER_GRAPH)
        set(DOXYGEN_CALLER_GRAPH YES)
    else()
        set(DOXYGEN_CALLER_GRAPH NO)
    endif()

    # Set input and output files
    set(DOXYGEN_CONFIG_IN ${CMAKE_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_CONFIG_OUT ${CMAKE_BINARY_DIR}/Doxyfile)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/doc)

    # Configure
    configure_file(${DOXYGEN_CONFIG_IN} ${DOXYGEN_CONFIG_OUT} @ONLY)

    # Run doxygen
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONFIG_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen in ${DOXYGEN_OUTPUT_DIR}"
        VERBATIM
    )

endif()

# QtCreator add all files
file(GLOB_RECURSE EXTRA_FILES src/*)
file(GLOB_RECURSE EXTRA_INCLUDES include/*)
add_custom_target(${PROJECT_NAME}_OTHER_FILES ALL WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} SOURCES ${EXTRA_FILES} ${EXTRA_INCLUDES})

# Include directories
include_directories(include)

# Add sources
add_subdirectory(src)

