#cmake_minimum_required(VERSION 3.0)
 
#PROJECT(IrisTools)
 
find_package(VTK REQUIRED)
 

set(IRIS_TOOLS 
ExtractUGRegions
Reclassify
ReplaceLabel
)



foreach(program ${IRIS_TOOLS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${VTK_LIBRARIES})
endforeach(program)



