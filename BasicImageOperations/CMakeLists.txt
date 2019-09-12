FIND_PACKAGE(ITK 5 REQUIRED)
FIND_PACKAGE(VTK 8.2)


INCLUDE(${ITK_USE_FILE})

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})



#these are all compiled the same way
set(ITK_PROGRAMS 
RotateImage3D
CaclulateEntropyFiberAngles
ExtractLabels
PixelCounter
CombineImageLabels
)

set(ITK_VTK_PROGRAMS 
ParticleDensityImage
)

foreach(program ${ITK_PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${ITK_LIBRARIES})
endforeach(program)

IF (VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})

  foreach(program ${ITK_VTK_PROGRAMS})
      add_executable(${program} ${MYCOMMON_FILES} ${program}.cpp)
      target_link_libraries(${program} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  endforeach(program)
ENDIF(VTK_FOUND)





