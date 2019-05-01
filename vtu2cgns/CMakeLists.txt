cmake_minimum_required(VERSION 2.8)

# Gui
#-------------------------------------------------------


find_package(VTK REQUIRED)

include(${VTK_USE_FILE})

include_directories( "${CMAKE_SOURCE_DIR}/Mesher/src/" )

set (CGNS_INSTALLATION_DIR ${tools_SOURCE_DIR}/extern/CGNS/build CACHE PATH "Path to installed CGNS" )




add_executable ( vtu2cgns vtu2cgns.cpp      )

include_directories( "${CGNS_INSTALLATION_DIR}/include/" )
target_link_libraries ( vtu2cgns   ${CGNS_INSTALLATION_DIR}/lib/libcgns.a ${VTK_LIBRARIES}                  )
