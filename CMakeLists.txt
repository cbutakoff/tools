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


option(BUILD_BASIC "Build Basic Tools" ON)
option(BUILD_CGAL "Build CGAL Tools" ON)
option(BUILD_FillSurfaceHoles "Build FillSurfaceHoles" ON)
option(BUILD_IrisTools "Build Iris Tools" ON)
option(BUILD_layer_Generation_code "Build Layer Generation" ON)
option(BUILD_LVMeshFromMRI "Build LVMeshFromMRI" ON)
option(BUILD_VolumetricLVMesh "Build VolumetricLVMesh" ON)

if(BUILD_BASIC)
	add_subdirectory(BasicMeshOperations)
endif(BUILD_BASIC)

if(BUILD_CGAL)
	add_subdirectory(cgal)
endif(BUILD_CGAL)

if(BUILD_FillSurfaceHoles)
	add_subdirectory(FillSurfaceHoles)
endif(BUILD_FillSurfaceHoles)

if(BUILD_IrisTools)
	add_subdirectory(IrisTools)
endif(BUILD_IrisTools)

if(BUILD_layer_Generation_code)
	add_subdirectory(layer_Generation_code)
endif(BUILD_layer_Generation_code)

if(BUILD_LVMeshFromMRI)
	add_subdirectory(LVMeshFromMRI)
endif(BUILD_LVMeshFromMRI)

if(BUILD_VolumetricLVMesh)
	add_subdirectory(VolumetricLVMesh)
endif(BUILD_VolumetricLVMesh)

