cmake_minimum_required(VERSION 2.8)

PROJECT(tools)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${tools_BINARY_DIR}/bin")
endif()


include_directories(Common)
SET(MYCOMMON_FILES ../Common/CommonTools.h ../Common/CommonTools.cpp)

set(MYPROJECTS
BasicMeshOperations
cgal
FillSurfaceHoles
IrisTools
#layer_Generation_code
LVMeshFromMRI
#VolumetricLVMesh
DistTransform
) 

foreach(program ${MYPROJECTS})
	option(BUILD_${program} "Build ${program}" ON)
	if(BUILD_${program})
		add_subdirectory(${program})
	endif(BUILD_${program})
endforeach(program)


