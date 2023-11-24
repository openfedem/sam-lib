C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRTRA (A,B,C,WA,MSPAR,MTREES,MSIFA,M,N,KSA,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRTRA                  GROUP 9 / PUBLIC
C
C     TASK:  to perform the matrix multiplication
C          C = BT*A*B    (BT = transpose of B)
C     where  A  is a symmetric 'sparse' matrix (KSA .gt. 0) or a
C     diagonal matrix (KSA .le. 0), and  B  is a rectangular matrix
C
C     ROUTINES CALLED/REFERENCED :  DDOT/SDOT (BLAS)
C                                   SPRPRM    (SAM-9)
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   97-02-04 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,KSA,LPU,M,N,  MSPAR(50),MTREES(*),MSIFA(*)
      DOUBLE PRECISION  A(*),B(M,N),C(N,N),WA(M)
C
      INTEGER           I,J
      DOUBLE PRECISION  DDOT
C
      EXTERNAL          DDOT
C ----------------------------------------------------------------------
      IERR = 0
C
      DO 50 J=1,N
         CALL SPRPRM (A,B(1,J),C,WA,
     &                MSPAR,MTREES,MSIFA,M,1,KSA,88,LPU,IERR)
         IF (IERR .LT. 0)  GOTO 90
C
         DO 40 I=J,N
            C(I,J) = DDOT(M,B(1,I),1,WA(1),1)
            C(J,I) = C(I,J)
   40    CONTINUE
   50 CONTINUE
      GO TO 100
C                                                ! error exit
   90 WRITE (LPU,690)
C
  100 RETURN
C
  690 FORMAT(///' *** ERROR RETURN from  SAM  library routine SPRTRA')
C
      END
