C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRDOF (MADOF,MINEX,NANOD,IDOF,LDOF,INOD)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRDOF                  GROUP 9 / PUBLIC
C
C     T A S K :  To determine the node identifier and local (node) dof
C                number corresponding to (global) dof number IDOF
C
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-01-11 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IDOF,INOD,LDOF,NANOD,
     &                  MADOF(*),MINEX(*)
C                                                ! Local variables
      INTEGER           I
C ----------------------------------------------------------------------
C
      INOD = 0
      LDOF = 0
      DO 10 I=1,NANOD
         IF (IDOF .LT. MADOF(I+1)) THEN
            INOD = MINEX(I)
            LDOF = IDOF - MADOF(I) + 1
            GO TO 100
         ENDIF
   10 CONTINUE
C
  100 RETURN
      END
