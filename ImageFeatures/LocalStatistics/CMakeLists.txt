cmake_minimum_required(VERSION 3.0)

PROJECT(LocalIntensity)

FIND_PACKAGE(ITK REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
FIND_PACKAGE(Armadillo REQUIRED)

include_directories(${ARMADILLO_INCLUDES})
INCLUDE(${ITK_USE_FILE})



ADD_EXECUTABLE(SampleLocalStatistics SampleLocalStatistics.cpp)
TARGET_LINK_LIBRARIES(SampleLocalStatistics ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ARMADILLO_LIBRARIES})


