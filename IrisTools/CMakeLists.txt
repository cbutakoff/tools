#cmake_minimum_required(VERSION 2.8)
 
#PROJECT(IrisTools)
 
find_package(VTK 7 REQUIRED)
include(${VTK_USE_FILE})
 

set(IRIS_TOOLS 
ExtractUGRegions
Reclassify
ReplaceLabel
)



foreach(program ${IRIS_TOOLS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${VTK_LIBRARIES})
endforeach(program)



