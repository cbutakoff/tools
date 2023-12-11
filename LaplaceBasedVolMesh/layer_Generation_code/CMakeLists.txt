#cmake_minimum_required(VERSION 3.0)

#PROJECT(tools)

FIND_PACKAGE(ITK REQUIRED)
FIND_PACKAGE(VTK REQUIRED)


INCLUDE(${ITK_USE_FILE})

set(VMTK_INCLUDES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
set(VMTK_LIBRARIES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
include_directories(${VMTK_INCLUDES})

link_directories(${VMTK_LIBRARIES})

set(VMTK_LIBS
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
TARGET_LINK_LIBRARIES(SmoothMeshThroughImage ${ITK_LIBRARIES} ${VTK_LIBRARIES}  ${VMTK_LIBS})


ADD_EXECUTABLE(GenerateAtrialLayers ${MYCOMMON_FILES} GenerateAtrialLayers.cpp)
TARGET_LINK_LIBRARIES(GenerateAtrialLayers  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${VMTK_LIBS})



