# SPDX-License-Identifier: Apache-2.0

# Library setup

set ( LIB_ID SAM )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( USE_INTEL_FORTRAN )
  set ( BLA_VENDOR Intel10_64lp )
elseif ( LINUX )
  if ( FTENV_WARNINGS )
    string ( APPEND CMAKE_Fortran_FLAGS " -Wno-compare-reals -Wno-function-elimination" )
  endif ( FTENV_WARNINGS )
  if ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 10.0 )
    string ( APPEND CMAKE_Fortran_FLAGS " -fallow-argument-mismatch" )
  endif ( CMAKE_C_COMPILER_VERSION VERSION_GREATER 10.0 )
endif ( USE_INTEL_FORTRAN )

find_package ( BLAS REQUIRED )

add_subdirectory ( SPR )


# Build and install

file ( GLOB F77_SOURCE_FILES *.f )
add_library ( ${LIB_ID} ${F77_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} ${BLAS_LIBRARIES} )
if ( "${PACKAGE_ID}" STREQUAL "SAM" )
  do_install ( ${LIB_ID} "${CMAKE_INSTALL_PREFIX}" )
endif ( "${PACKAGE_ID}" STREQUAL "SAM" )
