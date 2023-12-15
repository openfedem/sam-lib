# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0

# Library setup

set ( LIB_ID SPR )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

option ( BUILD_SPR_INTEGER8 "Build the SPR library using 64-bit integer" ON )
mark_as_advanced ( BUILD_SPR_INTEGER8 )

if ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_Fortran_FLAGS " -Wno-unused-parameter" )
endif ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )

if ( BUILD_SPR_INTEGER8 AND PLATFORM_BITSIZE_IS_64_BIT )
  if ( WIN )
    string ( APPEND CMAKE_Fortran_FLAGS " -4I8" )
  elseif ( LINUX )
    if ( USE_INTEL_FORTRAN )
      string ( APPEND CMAKE_Fortran_FLAGS " -i8" )
    else ( USE_INTEL_FORTRAN )
      string ( APPEND CMAKE_Fortran_FLAGS " -fdefault-integer-8" )
    endif ( USE_INTEL_FORTRAN )
  endif ( WIN )
  set ( LIB_ID SPR_I8 )
endif ( BUILD_SPR_INTEGER8 AND PLATFORM_BITSIZE_IS_64_BIT )


# Build and install

file ( GLOB F77_SOURCE_FILES *.f )
add_library ( ${LIB_ID} ${F77_SOURCE_FILES} )
target_link_libraries ( ${LIB_ID} SAM )
if ( "${PACKAGE_ID}" STREQUAL "SAM" )
  install ( TARGETS ${LIB_ID}
            LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" )
endif ( "${PACKAGE_ID}" STREQUAL "SAM" )
