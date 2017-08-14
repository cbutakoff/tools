#cmake_minimum_required(VERSION 2.8)

#PROJECT(tools)

FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(VTK REQUIRED)


INCLUDE(${VTK_USE_FILE})
INCLUDE(${ITK_USE_FILE})

set(VMTK_INCLUDES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
set(VMTK_LIBRARIES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
include_directories(${VMTK_INCLUDES})

link_directories(${VMTK_LIBRARIES})

set(VMTK_LIBRARIES 
vtkvmtkCommon
vtkvmtkIO
vtkvmtkITK
vtkvmtkComputationalGeometry
vtkvmtkDifferentialGeometry
vtkvmtkSegmentation
#vtkvmtkRendering
vtkvmtkContrib
vtkvmtkMisc
)



#SET(MYCOMMON_FILES ../Common/CommonTools.h ../Common/CommonTools.cpp)


ADD_EXECUTABLE(SmoothMeshThroughImage ${MYCOMMON_FILES} SmoothMeshTrhoughImage.cpp)
TARGET_LINK_LIBRARIES(SmoothMeshThroughImage ${ITK_LIBRARIES} ${VTK_LIBRARIES}  ${VMTK_LIBRARIES})


ADD_EXECUTABLE(GenerateAtrialLayers ${MYCOMMON_FILES} GenerateAtrialLayers.cpp)
TARGET_LINK_LIBRARIES(GenerateAtrialLayers  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${VMTK_LIBRARIES})



