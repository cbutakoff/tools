# Enable ExternalProject CMake module
include(ExternalProject)

################
# Define a test
add_executable(SurfaceHoleFiller_test SurfaceHoleFiller_test.cpp)

ExternalProject_Add(googletest
    GIT_REPOSITORY    https://github.com/google/googletest.git
    GIT_TAG           master
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/googletest-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/googletest-build"
    CMAKE_ARGS -DBUILD_GTEST=ON
	       -DBUILD_SHARED_LIBS=ON
               -Dgtest_force_shared_crt=ON
    # Disable install step
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
#    LOG_DOWNLOAD ON
#    LOG_CONFIGURE ON
#    LOG_BUILD ON
)


# Specify include dir for googlemock and googletest
ExternalProject_Get_Property(googletest source_dir binary_dir)
include_directories(${source_dir}/googletest/include)


######################################
# Configure the test to use GoogleTest
#
# If used often, could be made a macro.

add_dependencies(SurfaceHoleFiller_test googletest)

target_link_libraries(SurfaceHoleFiller_test ${binary_dir}/googlemock/gtest/libgtest.a)
target_link_libraries(SurfaceHoleFiller_test ${binary_dir}/googlemock/gtest/libgtest_main.a)
target_link_libraries(SurfaceHoleFiller_test pthread)

##################################
# Just make the test runnable with
#   $ make test

enable_testing()
add_test(NAME    SurfaceHoleFiller_test 
         COMMAND SurfaceHoleFiller_test)

