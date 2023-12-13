C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRMXD ( MADOF , MPMNPC,
     $                    MMNPC , MPAR    )
C
C @(#)sprmxd.f 1.1 95/12/29
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRMXD                  GROUP 9 / PUBLIC
C
C     T A S K :  To determine the order of the largest finite element
C                in the structure model and place it in MPAR(20).
C
C
C     ROUTINES CALLED/REFERENCED :    MAX    (FORTRAN library)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-02-19 / 1.0
C
C **********************************************************************
C
C     IMPLICIT          NONE
      INTEGER           MADOF(*), MMNPC(*), MPMNPC(*), MPAR(50)
C
C                                                ! local variables
      INTEGER           IEL,INOD,J,MAXDOF,NEDOF
C ----------------------------------------------------------------------
C
      MAXDOF = 0
      DO 20 IEL=1,MPAR(2)
         NEDOF = 0
         DO 10 J=MPMNPC(IEL),MPMNPC(IEL+1)-1
            INOD  = MMNPC(J)
            NEDOF = NEDOF + MADOF(INOD+1) - MADOF(INOD)
   10    CONTINUE
         MAXDOF = MAX(MAXDOF,NEDOF)
   20 CONTINUE
C
      MPAR(20) = MAXDOF
C
      RETURN
      END

