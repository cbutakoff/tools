cmake_minimum_required(VERSION 3.0)
project (boost_zernike)

FIND_PACKAGE(Boost COMPONENTS serialization  REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

FIND_PACKAGE(Armadillo)
if(ARMADILLO_FOUND)
	include_directories(${ARMADILLO_INCLUDES})
endif(ARMADILLO_FOUND)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=500") 


ADD_EXECUTABLE(boost_zernike  boost_zernike.cxx defs.h)
TARGET_LINK_LIBRARIES(boost_zernike   ${Boost_LIBRARIES})

ADD_EXECUTABLE(generate_inputs generate_inputs.cxx)
TARGET_LINK_LIBRARIES(generate_inputs)

ADD_EXECUTABLE(decompose_reconstruct  decompose_reconstruct.cxx)
TARGET_LINK_LIBRARIES(decompose_reconstruct ${ITK_LIBRARIES}  ${Boost_LIBRARIES} ${ARMADILLO_LIBRARIES} )




