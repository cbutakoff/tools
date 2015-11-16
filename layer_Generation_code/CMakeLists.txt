cmake_minimum_required(VERSION 2.8)

PROJECT(tools)

FIND_PACKAGE(ITK REQUIRED)

FIND_PACKAGE(VTK REQUIRED)


INCLUDE(${VTK_USE_FILE})
INCLUDE(${ITK_USE_FILE})

set(VMTK_INCLUDES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
include_directories(${VMTK_INCLUDES})


set(VMTK_LIBRARIES 
vtkvmtkCommon
vtkvmtkIO
vtkvmtkITK
vtkvmtkComputationalGeometry
vtkvmtkDifferentialGeometry
vtkvmtkSegmentation
vtkvmtkRendering
vtkvmtkContrib
vtkvmtkMisc
)



SET(MYCOMMON_FILES CommonTools.h CommonTools.cpp)


ADD_EXECUTABLE(SmoothMeshThroughImage ${MYCOMMON_FILES} SmoothMeshTrhoughImage.cpp)
TARGET_LINK_LIBRARIES(SmoothMeshThroughImage ${ITK_LIBRARIES} ${VTK_LIBRARIES}  ${VMTK_LIBRARIES})


ADD_EXECUTABLE(GenerateAtrialLayers ${MYCOMMON_FILES} GenerateAtrialLayers.cpp)
TARGET_LINK_LIBRARIES(GenerateAtrialLayers  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${VMTK_LIBRARIES})



