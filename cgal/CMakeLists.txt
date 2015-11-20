project( Mesh_3_applications )

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)



FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )


  include( CGAL_CreateSingleSourceCGALProgram )

# Compilable examples
  if( CGAL_ImageIO_USE_ZLIB )
    	create_single_source_cgal_program( "image2mesh.cpp" )
  else()
    message( STATUS "NOTICE: The example mesh_3D_image.cpp needs CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
  endif()
  create_single_source_cgal_program( "poly2mesh.cpp" )
  create_single_source_cgal_program( "segment_mesh_sdf.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

ADD_EXECUTABLE(vtk2inrimage  vtk2inrimage.cpp)
TARGET_LINK_LIBRARIES(vtk2inrimage  ${VTK_LIBRARIES} )


