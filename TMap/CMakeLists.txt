cmake_minimum_required(VERSION 3.0)

PROJECT(TMap)

FIND_PACKAGE(VTK REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

ADD_EXECUTABLE(TMap TMap.h TMap.cpp main.cpp)
TARGET_LINK_LIBRARIES(TMap ${VTK_LIBRARIES} )


