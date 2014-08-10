cmake_minimum_required( VERSION 2.8 )

set( CMAKE_CXX_FLAGS "-std=c++0x -u MAIN__" )

# PCL
find_package( PCL )
include_directories( ${PCL_INCLUDE_DIRS} )
link_directories( ${PCL_LIBRARY_DIRS} )
add_definitions( ${PCL_DEFINITIONS} )

# levmar
add_subdirectory( levmar )


# Ceres
find_package( Ceres )
include_directories( ${CERES_INCLUDE_DIRS} )

include_directories( ${CMAKE_SOURCE_DIR}/include )
include_directories( ${CMAKE_SOURCE_DIR} )
link_directories( ${CMAKE_SOURCE_DIR}/lib )

# SQ_utils library
add_library( SQ_utils src/SQ_utils.cpp src/analytic_equations.cpp )
target_link_libraries( SQ_utils levmar ${PCL_LIBRARIES} ${CERES_LIBRARIES} pthread )
set_target_properties( SQ_utils PROPERTIES RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_SOURCE_DIR}/bin ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )


# SQ_fitter library
add_library( SQ_fitter src/SQ_fitter.cpp )
target_link_libraries( SQ_fitter SQ_utils levmar ${PCL_LIBRARIES} ${CERES_LIBRARIES} pthread )
set_target_properties( SQ_fitter PROPERTIES RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_SOURCE_DIR}/bin ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )

# Add tests
add_subdirectory( tests )