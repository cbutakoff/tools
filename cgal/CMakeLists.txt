project( Mesh_3_applications )

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)



FIND_PACKAGE(VTK REQUIRED)
find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  include_directories(/usr/include/eigen3/) 
  include( CGAL_CreateSingleSourceCGALProgram )

# Compilable examples
  #if( CGAL_ImageIO_USE_ZLIB )
  #    	create_single_source_cgal_program( "image2mesh.cpp" )
  #else()
  #    message( STATUS "NOTICE: The example mesh_3D_image.cpp needs CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
  #endif()
  #create_single_source_cgal_program( "poly2mesh.cpp" )
  create_single_source_cgal_program( "segment_mesh_sdf.cpp" )
  #create_single_source_cgal_program( "ls_conformal.cpp" )
  #create_single_source_cgal_program( "poisson_recon.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

ADD_EXECUTABLE(vtk2inrimage  vtk2inrimage.cpp)
TARGET_LINK_LIBRARIES(vtk2inrimage  ${VTK_LIBRARIES} )


