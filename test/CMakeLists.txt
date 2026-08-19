# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0

# Build setup

set ( LIB_ID testSAM )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

if ( pFUnit_FOUND )
# Add some unit tests using pFUnit (executed via ctest)
  enable_fortran_tests ()
  add_fortran_test ( test_SAM testSAM.pf SAM )
endif ( pFUnit_FOUND )
