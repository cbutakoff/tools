cmake_minimum_required(VERSION 2.8)

PROJECT(tools)


FIND_PACKAGE(VTK REQUIRED)


set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Release
 RelWithDebInfo MinSizeRel."
       FORCE)

INCLUDE(${VTK_USE_FILE})


# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)



FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})


ADD_EXECUTABLE(FillSurfaceHoles  FillSurfaceHoles.cpp)
TARGET_LINK_LIBRARIES(FillSurfaceHoles ${VTK_LIBRARIES} ${Boost_LIBRARIES})

