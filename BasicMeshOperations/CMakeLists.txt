FIND_PACKAGE(ITK 5)
FIND_PACKAGE(VTK 8.2 REQUIRED)

find_package(Eigen3)


#FIND_PACKAGE(Armadillo)
#if(ARMADILLO_FOUND)
#	include_directories(${ARMADILLO_INCLUDES})
#endif(ARMADILLO_FOUND)


INCLUDE(${VTK_USE_FILE})




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
SampleUG2Image
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
SampleImage2UG
InterpolateFibers
RelabelSurfaceMesh
ThresholdDS
CalculateEdgeLengths
UniformityIndex
)

set(VTKEIGEN_PROGRAMS 
EstimateConductionVelocity
CalculateFlatteningDistortion
)



find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


IF (ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})

    foreach(program ${ITK_VTK_PROGRAMS})
        add_executable(${program} ${MYCOMMON_FILES} ${program}.cpp)
        target_link_libraries(${program} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
    endforeach(program)
ENDIF(ITK_FOUND)

foreach(program ${VTK_PROGRAMS})
  add_executable(${program} ${MYVTKCOMMON_FILES} ${program}.cpp)
  target_link_libraries(${program} ${VTK_LIBRARIES})
endforeach(program)


IF (EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})


    foreach(program ${VTKEIGEN_PROGRAMS})
        add_executable(${program} ${MYVTKCOMMON_FILES} ${program}.cpp)
        target_link_libraries(${program} ${VTK_LIBRARIES})
    endforeach(program)
ENDIF (EIGEN3_FOUND)



