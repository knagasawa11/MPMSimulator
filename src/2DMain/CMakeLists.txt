include( CMakeSourceFiles.txt )

# Track the Git commit within the source
include( GetGitRevisionDescription )
get_git_head_revision( GIT_REFSPEC GIT_SHA1 )

# Track the build time
string( TIMESTAMP BUILD_DATE_TIME "%Y-%m-%d %H:%M:%S" )

# Create a file with build information
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h" @ONLY )

# if( CMAKE_Fortran_COMPILER )
#   set_source_files_properties( CompileDefinitions.h CompileDefinitions.cpp PROPERTIES COMPILE_DEFINITIONS FORTRAN_FOUND )
# endif()

# if( USE_HDF5 )
#   set( Headers "${Headers};HDF5File.h" )
#   set( Sources "${Sources};HDF5File.cpp" )
# endif()

add_executable( PantaRhei3D ${Headers} ${Sources} )
# if( ENABLE_IWYU )
#   set_property( TARGET scisim PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path} )
# endif()

# Add a reliable indicator of the current compiler
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  target_compile_definitions( PantaRhei3D PUBLIC CMAKE_DETECTED_GNU_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  target_compile_definitions( PantaRhei3D PUBLIC CMAKE_DETECTED_INTEL_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  target_compile_definitions( PantaRhei3D PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  target_compile_definitions( PantaRhei3D PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()

# Eigen3 is header only and required for callers
find_package( Eigen3 3.2.9 REQUIRED )
target_include_directories( PantaRhei3D SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )
# Disable Eigen's built in multithreading as we will use OpenMP around Eigen
target_compile_definitions( PantaRhei3D PUBLIC EIGEN_DONT_PARALLELIZE )
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  # EIGEN_DONT_ALIGN also enables EIGEN_DONT_VECTORIZE
  target_compile_definitions( PantaRhei3D PUBLIC EIGEN_DONT_ALIGN EIGEN_INITIALIZE_MATRICES_BY_NAN )
endif()

if( USE_HDF5 )
  #find_package( HDF5 REQUIRED C )
  find_package( HDF5 REQUIRED )
  if( NOT ${HDF5_FOUND} )
    message( FATAL_ERROR "HDF5 not found." )
  endif()
  target_include_directories( PantaRhei3D SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS} )
  target_compile_definitions( PantaRhei3D PUBLIC USE_HDF5 )
  # target_link_libraries( PantaRhei3D INTERFACE ${HDF5_C_LIBRARIES} )
  target_link_libraries( PantaRhei3D ${HDF5_LIBRARIES} )
endif()

# if( USE_OPENMP )
#   find_package( OpenMP )
#   if( NOT OPENMP_FOUND )
#     message( FATAL_ERROR "Error, failed to locate OpenMP." )
#   endif()
#   # target_compile_options( PantaRhei3D PRIVATE ${OpenMP_CXX_FLAGS} )
#   set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp" )
#   # target_link_libraries( PantaRhei3D ${OpenMP_LIBRARIES} )
# endif()

find_package( RapidXML REQUIRED )
target_include_directories( PantaRhei3D SYSTEM PRIVATE ${RAPIDXML_INCLUDE_DIR} )
