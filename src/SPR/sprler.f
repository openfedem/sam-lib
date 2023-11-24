C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLER (N,I1,I2,I3,LPU,IERR)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLER                  GROUP 9 / PRIVATE
C
C     TASK :  To print error/warning messages and set the error flag
C             for eigenproblem routine SPRLAN
C
C     ROUTINES CALLED/REFERENCED :  None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   96-12-30 / 1.0
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           IERR,I1,I2,I3,LPU,N
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         IF (N .LT. 10) THEN
C                                                ! warnings
            WRITE (LPU,600)
            IF (N .EQ. 1) THEN
               WRITE (LPU,601) I1,I2,I3
            ELSEIF (N .EQ. 2) THEN
               WRITE (LPU,602)
            ELSEIF (N .EQ. 3) THEN
               WRITE (LPU,603) I1
            ENDIF
         ELSE
C                                                ! errors
            WRITE (LPU,690)
            IF (N.GT.9 .AND. N.LT.21) THEN
               WRITE (LPU,610)
               IF (N .EQ. 11)  WRITE (LPU,611) I1,I2
               IF (N .EQ. 12)  WRITE (LPU,612) I1
               IF (N .EQ. 13)  WRITE (LPU,613) I1
               IF (N .EQ. 14)  WRITE (LPU,614) I1
               IF (N .EQ. 15)  WRITE (LPU,615) I1,I2,I3
               IF (N .EQ. 16)  WRITE (LPU,616) I1,I2
               IF (N .EQ. 17)  WRITE (LPU,617) I1,I2
               IF (N .EQ. 18)  WRITE (LPU,618) I1,I2,I3
               IF (N .EQ. 19)  WRITE (LPU,619) I1
               IF (N .EQ. 20)  WRITE (LPU,620)
            ELSEIF (N .EQ. 21) THEN
               WRITE (LPU,621)
            ELSEIF (N .EQ. 22) THEN
               WRITE (LPU,622) I1
            ELSEIF (N .EQ. 23) THEN
               WRITE (LPU,623)
            ELSEIF (N .EQ. 24) THEN
               WRITE (LPU,624) I1
            ELSEIF (N .EQ. 25) THEN
               WRITE (LPU,625) I1
            ELSEIF (N .EQ. 26) THEN
               WRITE (LPU,626) I1
            ELSEIF (N .EQ. 27) THEN
               WRITE (LPU,627)
            ELSEIF (N .EQ. 31) THEN
               WRITE (LPU,631) I1
            ELSEIF (N .EQ. 32) THEN
               WRITE (LPU,632) I1
            ELSEIF (N .EQ. 33) THEN
               WRITE (LPU,610)
               WRITE (LPU,633) I1,I2
            ELSEIF (N .EQ. 34) THEN
               WRITE (LPU,610)
               WRITE (LPU,634) I1,I2
            ENDIF
         ENDIF
      ENDIF
C
      IF (N.LT.10) THEN
         IERR = 1
      ELSE IF (N.LT.21 .OR. N.EQ.33 .OR. N.EQ.34) THEN
         IERR = -1
      ELSE IF (N.EQ.21 .OR. N.EQ.22) THEN
         IF (IERR .EQ. -3) THEN
            IERR = 19-N
         ELSE
            IERR = -10
         ENDIF
      ELSE IF (N.EQ.26 .OR. N.EQ.27) THEN
         IERR = -7
      ELSE IF (N.LT.31) THEN
         IERR = 19-N
      ELSE
         IERR = 23-N
      ENDIF
C
      RETURN
C ----------------------------------------------------------------------
  600 FORMAT(///' *** WARNING from  S A M  library routine SPRLAN')
  601 FORMAT(5X,'Only',I5,' eigenvalues accepted'/
     &       5X,'in',I5,' Lanczos steps'/
     &       5X,'(the required number was',I5,')')
  602 FORMAT(5X,'Nonzero shift may render A non-positive definite'/
     &       5X,'and cause an error during Cholesky factorization' )
  603 FORMAT(5X,'Problem has only',I5,' finite eigenvalues')
  610 FORMAT(5X,'Illegal or inconsistent input parameters')
  611 FORMAT(5X,'Storage codes :  KSA =',I3,5X,'KSB =',I3)
  612 FORMAT(5X,'Reorthogonalization code =',I5)
  613 FORMAT(5X,'Parameter KEX (=',I5,') out of range')
  614 FORMAT(5X,'Parameter KPD (=',I5,') out of range')
  615 FORMAT(5X,'Check parameters KPD, KSA and KSB =',3I6)
  616 FORMAT(5X,'Check parameters KSB and KSR =',2I6)
  617 FORMAT(5X,'Check parameters KSA and NEVAL =',2I6)
  618 FORMAT(5X,'Check parameters NEVAL, NEVEC and MAXLAN =',3I6)
  619 FORMAT(5X,'KSA =',I3,' is inconsistent with non-zero shift')
  620 FORMAT(5X,'Inconsistent sparse control information')
C
  621 FORMAT(5X,'Error during factorization of matrix A')
  622 FORMAT(5X,'Matrix B (KSB =',I3,') cannot be factorized')
  623 FORMAT(5X,'Problem has no finite eigenvalues (N=1)')
  624 FORMAT(5X,'Error during premultiplication (by matrix H)'/
     &       5X,'Lanczos step no.',I5)
  625 FORMAT(5X,'Error during solution of small, tridiag. eigenproblem'
     &      /5X,'(of dimension',I5,')')
  626 FORMAT(5X,'Error during recovery of eigenvector no.',I5)
  627 FORMAT(5X,'Error during recovery of eigenvectors')
  631 FORMAT(5X,'Breakdown! - 11th restart required in step no.',I5)
  632 FORMAT(5X,'No eigenvalues accepted (after',I5,' Lanczos steps)')
  633 FORMAT(5X,'Matrix A has a different dimension than matrix B'/
     &       5X,'NEQA   =',I10,5X,'NEQB   =',I10)
  634 FORMAT(5X,'Matrix A has a different sparsity pattern than B'/
     &       5X,'NZEROA =',I10,5X,'NZEROB =',I10/
     &       5X,'Non-zero shift is then not allowed')
C
  690 FORMAT(///' *** ERROR return from  S A M  library routine SPRLAN')
C
      END
