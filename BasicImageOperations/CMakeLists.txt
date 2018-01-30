FIND_PACKAGE(ITK 4 REQUIRED)


INCLUDE(${ITK_USE_FILE})

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})



#these are all compiled the same way
set(ITK_PROGRAMS 
RotateImage3D
CaclulateEntropyFiberAngles
)


foreach(program ${ITK_PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${ITK_LIBRARIES})
endforeach(program)





