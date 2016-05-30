cmake_minimum_required(VERSION 2.8)
 
PROJECT(IrisTools)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(ExtractUGRegions ExtractUGRegions)
target_link_libraries(ExtractUGRegions ${VTK_LIBRARIES})

