cmake_minimum_required(VERSION 2.8)

PROJECT(tools)

FIND_PACKAGE(ITK REQUIRED)

FIND_PACKAGE(VTK REQUIRED)


set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: Debug Release
 RelWithDebInfo MinSizeRel."
       FORCE)

FIND_PACKAGE(Armadillo)
if(ARMADILLO_FOUND)
	include_directories(${ARMADILLO_INCLUDES})
endif(ARMADILLO_FOUND)


INCLUDE(${VTK_USE_FILE})
INCLUDE(${ITK_USE_FILE})



#include_directories(imported) 
#include_directories(imported/pgm2itkvol) 



# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

SET(MYCOMMON_FILES CommonTools.h CommonTools.cpp)

ADD_EXECUTABLE(CreateImageMask ${MYCOMMON_FILES} CreateImageMask.cpp)
TARGET_LINK_LIBRARIES(CreateImageMask ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(MeshHeart ${MYCOMMON_FILES} MeshHeart.cpp)
TARGET_LINK_LIBRARIES(MeshHeart ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


ADD_EXECUTABLE(ExtractShapeRegion ${MYCOMMON_FILES} ExtractShapeRegion.cpp)
TARGET_LINK_LIBRARIES(ExtractShapeRegion ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(GenerateVolumetricLVMesh ${MYCOMMON_FILES} GenerateVolumetricLVMesh.cpp)
TARGET_LINK_LIBRARIES(GenerateVolumetricLVMesh ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(ResampleImage ${MYCOMMON_FILES} ResampleImage.cpp)
TARGET_LINK_LIBRARIES(ResampleImage ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(VTKConvert ${MYCOMMON_FILES} VTKConvert.cpp)
TARGET_LINK_LIBRARIES(VTKConvert ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(MRIRemesh ${MYCOMMON_FILES} MRIRemesh.cpp)
TARGET_LINK_LIBRARIES(MRIRemesh ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(MakeBiventricularMesh ${MYCOMMON_FILES} MakeBiventricularMesh.cpp)
TARGET_LINK_LIBRARIES(MakeBiventricularMesh ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(LabelBiventricularMesh ${MYCOMMON_FILES} LabelBiventricularMesh.cpp)
TARGET_LINK_LIBRARIES(LabelBiventricularMesh ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(CloseBVMesh ${MYCOMMON_FILES} CloseBVMesh.cpp)
TARGET_LINK_LIBRARIES(CloseBVMesh ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(SetScalars ${MYCOMMON_FILES} SetScalars.cpp)
TARGET_LINK_LIBRARIES(SetScalars ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(PassScalars ${MYCOMMON_FILES} PassScalars.cpp)
TARGET_LINK_LIBRARIES(PassScalars ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(PassScalarsReverse ${MYCOMMON_FILES} PassScalarsReverse.cpp)
TARGET_LINK_LIBRARIES(PassScalarsReverse ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


ADD_EXECUTABLE(PassScalarsInterp ${MYCOMMON_FILES} PassScalarsInterp.cpp)
TARGET_LINK_LIBRARIES(PassScalarsInterp ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


#ADD_EXECUTABLE(TransformPhilipsHeart ${MYCOMMON_FILES} TransformPhilipsHeart)
#TARGET_LINK_LIBRARIES(TransformPhilipsHeart ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(Blob ${MYCOMMON_FILES} Blob.cpp)
TARGET_LINK_LIBRARIES(Blob ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


#ADD_EXECUTABLE(LabelBranches3D imported/LabelBranches3D.cpp imported/itkImageToVTKImageFilter.h imported/itkImageToVTKImageFilter.txx imported/itkBinaryThinningImageFilter3D.h imported/itkBinaryThinningImageFilter3D.txx)
#TARGET_LINK_LIBRARIES(LabelBranches3D ${ITK_LIBRARIES} )

#gives an ITK error
#ADD_EXECUTABLE(itkvol2pgm  imported/pgm2itkvol/itkvol2pgm.cxx 
#imported/pgm2itkvol/mccodimage.h
#imported/pgm2itkvol/mcimage.c
#imported/pgm2itkvol/mcimage.h
#imported/pgm2itkvol/mcutil.h)
#c++TARGET_LINK_LIBRARIES(itkvol2pgm ${ITK_LIBRARIES}  )

#ADD_EXECUTABLE(pgm2itkvol imported/pgm2itkvol/pgm2itkvol.cxx 
#imported/pgm2itkvol/mccodimage.h
#imported/pgm2itkvol/mcimage.c
#imported/pgm2itkvol/mcimage.h
#imported/pgm2itkvol/mcutil.h)
#TARGET_LINK_LIBRARIES(pgm2itkvol ${ITK_LIBRARIES}  )


ADD_EXECUTABLE(MeshSegmentationLaplace ${MYCOMMON_FILES} MeshSegmentationLaplace.cpp)
TARGET_LINK_LIBRARIES(MeshSegmentationLaplace ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


ADD_EXECUTABLE(chaste2vtk ${MYCOMMON_FILES} chaste2vtk.cpp)
TARGET_LINK_LIBRARIES(chaste2vtk  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(Morphology3D ${MYCOMMON_FILES} Morphology3D.cpp)
TARGET_LINK_LIBRARIES(Morphology3D  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )



ADD_EXECUTABLE(nrrd2vtk ${MYCOMMON_FILES} nrrd2vtk.cpp)
TARGET_LINK_LIBRARIES(nrrd2vtk  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(ComputeLocalMeshDensity ${MYCOMMON_FILES} ComputeLocalMeshDensity.cpp)
TARGET_LINK_LIBRARIES(ComputeLocalMeshDensity  ${ITK_LIBRARIES} ${VTK_LIBRARIES} )





ADD_EXECUTABLE(S2SError ${MYCOMMON_FILES} S2SError.cpp)
TARGET_LINK_LIBRARIES(S2SError  ${ITK_LIBRARIES}  ${VTK_LIBRARIES})


if(ARMADILLO_FOUND)
	ADD_EXECUTABLE(zernike  zernike.cpp)
	TARGET_LINK_LIBRARIES(zernike  ${VTK_LIBRARIES} ${ARMADILLO_LIBRARIES})

	ADD_EXECUTABLE(GeneralizedRicciFlow ${MYCOMMON_FILES} GeneralizedRicciFlow.cpp)
	TARGET_LINK_LIBRARIES(GeneralizedRicciFlow  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ARMADILLO_LIBRARIES})

	include_directories("/usr/include/eigen3/")
	ADD_EXECUTABLE(ls_conformal ${MYCOMMON_FILES} ls_conformal.cpp)
		TARGET_LINK_LIBRARIES(ls_conformal ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${ARMADILLO_LIBRARIES} )

endif(ARMADILLO_FOUND)


ADD_EXECUTABLE(MarchingCubes ${MYCOMMON_FILES} MarchingCubes.cpp)
TARGET_LINK_LIBRARIES(MarchingCubes  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(CheckTetraMesh ${MYCOMMON_FILES} CheckTetraMesh.cpp)
TARGET_LINK_LIBRARIES(CheckTetraMesh  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(PassLabels2TetMesh ${MYCOMMON_FILES} PassLabels2TetMesh.cpp)
TARGET_LINK_LIBRARIES(PassLabels2TetMesh  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(tetmesh_uniform ${MYCOMMON_FILES} tetmesh_uniform.cpp)
TARGET_LINK_LIBRARIES(tetmesh_uniform  ${ITK_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(ProjectImage2Mesh ProjectImage2Mesh.cpp)
TARGET_LINK_LIBRARIES(ProjectImage2Mesh  ${VTK_LIBRARIES})

ADD_EXECUTABLE(GenerateMeshHoleCovers CommonTools.cpp CommonTools.h GenerateMeshHoleCovers.cpp)
TARGET_LINK_LIBRARIES(GenerateMeshHoleCovers ${ITK_LIBRARIES} ${VTK_LIBRARIES} )


ADD_EXECUTABLE(heap_test MinHeap/MinHeap.h MinHeap/MinHeap.txx MinHeap/heap_test.cpp)
TARGET_LINK_LIBRARIES(heap_test )

ADD_EXECUTABLE(AtriumSegment  MinHeap/MinHeap.h MinHeap/MinHeap.txx AtriumSegment.cpp)
TARGET_LINK_LIBRARIES(AtriumSegment ${VTK_LIBRARIES} )


option(USE_BOOST "Boost dependent programs will be compiled" ON)

IF(USE_BOOST)
	FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)

	include_directories(${Boost_INCLUDE_DIRS})

	ADD_EXECUTABLE(AdvFrontHoleFilling  AdvFrontHoleFilling.cpp)
	TARGET_LINK_LIBRARIES(AdvFrontHoleFilling ${VTK_LIBRARIES} ${Boost_LIBRARIES})
ENDIF(USE_BOOST)


option(USE_VMTK "Vmtk dependent programs will be compiled" ON)


if(USE_VMTK)
	set(VMTK_INCLUDES "/usr/include/vmtk" CACHE FILEPATH "VMTK Includes")
	include_directories(${VMTK_INCLUDES})


	set(VMTK_LIBRARIES 
	vtkvmtkCommon
	vtkvmtkIO
	vtkvmtkITK
	vtkvmtkComputationalGeometry
	vtkvmtkDifferentialGeometry
	vtkvmtkSegmentation
	vtkvmtkRendering
	vtkvmtkContrib
	vtkvmtkMisc
	)

	#for older vmtk
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

	ADD_EXECUTABLE(UniformRemesh ${MYCOMMON_FILES} UniformRemesh.cpp)
	TARGET_LINK_LIBRARIES(UniformRemesh  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${VMTK_LIBRARIES})

	ADD_EXECUTABLE(GenerateAtrialLayers ${MYCOMMON_FILES} GenerateAtrialLayers.cpp)
	TARGET_LINK_LIBRARIES(GenerateAtrialLayers  ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${VMTK_LIBRARIES})

	ADD_EXECUTABLE(SmoothMeshThroughImage ${MYCOMMON_FILES} SmoothMeshTrhoughImage.cpp)
	TARGET_LINK_LIBRARIES(SmoothMeshThroughImage ${ITK_LIBRARIES} ${VTK_LIBRARIES}  ${VMTK_LIBRARIES})
endif(USE_VMTK)


