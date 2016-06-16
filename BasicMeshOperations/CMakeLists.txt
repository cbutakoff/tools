FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(VTK 7 REQUIRED)

FIND_PACKAGE(Armadillo)
if(ARMADILLO_FOUND)
	include_directories(${ARMADILLO_INCLUDES})
endif(ARMADILLO_FOUND)


INCLUDE(${VTK_USE_FILE})
INCLUDE(${ITK_USE_FILE})




# add a target to generate API documentation with Doxygen
#find_package(Doxygen)
#if(DOXYGEN_FOUND)
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#add_custom_target(doc
#${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
#WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#COMMENT "Generating API documentation with Doxygen" VERBATIM
#)
#endif(DOXYGEN_FOUND)

#these are all compiled the same way
set(ITK_VTK_PROGRAMS 
CreateImageMask
ExtractShapeRegion
ResampleImage
VTKConvert
SetScalars
PassScalars
PassScalarsReverse
PassScalarsInterp
chaste2vtk
Morphology3D
nrrd2vtk
S2SError
MarchingCubes
CheckTetraMesh
PassLabels2TetMesh
)

set(VTK_PROGRAMS 
ProjectImage2Mesh
ensight2vtk
)


foreach(program ${ITK_VTK_PROGRAMS})
  add_executable(${program} ${MYCOMMON_FILES} ${program}.cpp)
  target_link_libraries(${program} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
endforeach(program)

foreach(program ${VTK_PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${VTK_LIBRARIES})
endforeach(program)




