# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required ( VERSION 2.8...3.10 )

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

include ( CheckPFUnit )
include ( FedemConfig )

if ( BUILD_TESTS )

  enable_testing ()
  set ( CTEST_OPTIONS --force-new-ctest-process --output-on-failure )
  if ( CMAKE_CONFIGURATION_TYPES )
    list ( APPEND CTEST_OPTIONS --build-config \"$<CONFIGURATION>\" )
  endif ( CMAKE_CONFIGURATION_TYPES )

  add_custom_target ( check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS} )
  include ( pFUnit ) # Using the pFUnit test framework for Fortran unit tests

endif ( BUILD_TESTS )

add_subdirectory ( src )
if ( BUILD_TESTS )
  add_subdirectory ( test )
endif ( BUILD_TESTS )
