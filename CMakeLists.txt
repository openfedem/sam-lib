# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required ( VERSION 2.8...3.5 )

# Project setup

set ( DOMAIN_ID FEDEM )
set ( PACKAGE_ID SAM )

project ( SAM Fortran )
message ( STATUS "Generating build project for ${PROJECT_SOURCE_DIR}" )

find_path ( _MODULES FedemConfig.cmake
            PATHS $ENV{CMAKE_MODULES}
           "${PROJECT_SOURCE_DIR}/cmake/Modules/"
          )
if ( _MODULES )
  message ( STATUS "NOTE : Using ${_MODULES}" )
  list ( APPEND CMAKE_MODULE_PATH ${_MODULES} )
else ( _MODULES )
  message ( STATUS "ERROR : Missing path to FedemConfig.cmake" )
  message ( FATAL_ERROR "Set environment variable CMAKE_MODULES and try again" )
endif ( _MODULES )
unset ( _MODULES CACHE )

include ( FedemConfig )

add_subdirectory ( src )
