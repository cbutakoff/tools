FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(VTK REQUIRED)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

#FIND_PACKAGE(Armadillo)
#if(ARMADILLO_FOUND)
#	include_directories(${ARMADILLO_INCLUDES})
#endif(ARMADILLO_FOUND)


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
ExtractShapeRegion
ResampleImage
nrrd2vtk
Morphology3D
chaste2vtk
S2SError
)

set(VTK_PROGRAMS 
CreateImageMask
ProjectImage2Mesh
ensight2vtk
SurfaceMeshBoolean
MarkVolmeshCellsInsideSurface
VTKConvert
SetScalars
PassScalars
PassScalarsReverse
PassScalarsInterp
MarchingCubes
CheckTetraMesh
PassLabels2TetMesh
CreateSurfaceMeshEnvelope
PassVolmeshLabels2Image
CreatePixelatedTetMesh
AddFibers2Mesh
EstimateConductionVelocity
SampleImage2UG
)


foreach(program ${ITK_VTK_PROGRAMS})
  add_executable(${program} ${MYCOMMON_FILES} ${program}.cpp)
  target_link_libraries(${program} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
endforeach(program)

foreach(program ${VTK_PROGRAMS})
  add_executable(${program} ${MYVTKCOMMON_FILES} ${program}.cpp)
  target_link_libraries(${program} ${VTK_LIBRARIES})
endforeach(program)




