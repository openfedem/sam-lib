C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLP3 (ALPHA,BETA,EVL,EVERR,MANCUR,
     +                   EPS,NLV,NGDEV,LPU,IFLAG)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLP3                 GROUP 9 / PRIVATE
C
C     TASK :  To print computed eigenvalues and their corresponding
C             relative errors after a given Lanczos step (no. NLV).
C             The step at which the eigenvalues were accepted is also
C             printed.
C             Optionally (IFLAG > 1) the tridiagonal form, stored in
C             ALPHA and BETA, is also printed
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-30 / 1.0
C
C **********************************************************************
C     CONDITIONALS  :   S - SINGLE PRECISION
C     (COLUMN 2)        D - DOUBLE PRECISION
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IFLAG,LPU,NGDEV,NLV,       MANCUR(NLV)
      DOUBLE PRECISION  EPS,  ALPHA(NLV),BETA(NLV),EVERR(NLV),EVL(NLV)
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (IFLAG .GT. 0) THEN
            WRITE (LPU,600)  NLV
            IF (IFLAG .GT. 1) THEN
               WRITE (LPU,610)
               DO 25 I=1,NLV
                  WRITE (LPU,620) I,BETA(I),ALPHA(I)
   25          CONTINUE
            ENDIF
            WRITE (LPU,630) NGDEV,EPS
         ENDIF
         WRITE (LPU,640)
         DO 50 I=1,NLV
            WRITE (LPU,650) I,EVL(I),EVERR(I),MANCUR(I)
   50    CONTINUE
      ENDIF
C
      RETURN
C ----------------------------------------------------------------------
C
  600 FORMAT(////5X,21('*')/5X,'Lanczos step no.',I5 / 5X,21('*'))
  610 FORMAT(//  5X,'Tridiagonal form :' ///
     &          12X,'Subdiagonal',6X,'diagonal'/)
  620 FORMAT(I7,1P2D16.5)
  630 FORMAT(//  5X,'Number of "good" eigenvalues :',I5,'   (tol =',
     &           1PE11.3,' )')
  640 FORMAT(// 15X,'COMPUTED',13X,'RELATIVE',6X,'ACCEPTED'/
     &          15X,'EIGENVALUES',10X,'ERRORS',8X,'AT STEP'/)
  650 FORMAT(I7,1PD24.15,D15.4,I9)
C
      END
