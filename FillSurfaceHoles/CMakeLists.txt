FIND_PACKAGE(VTK 7 REQUIRED)
INCLUDE(${VTK_USE_FILE})


FIND_PACKAGE(Boost COMPONENTS serialization REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})


ADD_EXECUTABLE(FillSurfaceHoles SurfaceHoleFiller.h  SurfaceHoleFiller.cpp FillSurfaceHoles.cpp)
TARGET_LINK_LIBRARIES(FillSurfaceHoles ${VTK_LIBRARIES} ${Boost_LIBRARIES})

