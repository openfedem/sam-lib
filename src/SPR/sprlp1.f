C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLP1 (NNZERO,MIP,SHIFT,N,NEVAL,NEVEC,MAXLAN,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLP1                  GROUP 9 / PRIVATE
C
C     TASK :  To print key input information to the eigensolution in
C             in subroutine SPRLAN
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
      INTEGER           LPU,MAXLAN,N,NEVAL,NEVEC,NNZERO,MIP(6)
      DOUBLE PRECISION  SHIFT
C
      INTEGER           I
C ----------------------------------------------------------------------
      IF (LPU .GT. 0) THEN
         WRITE (LPU,600)
         IF (NNZERO .GT. 0) THEN
            WRITE (LPU,610) N,NNZERO
         ELSE
            WRITE (LPU,611) N
         ENDIF
         I = ABS(MIP(2))
         IF (I .EQ. 1)  WRITE (LPU,620)
         IF (I .EQ. 2)  WRITE (LPU,630)
         WRITE (LPU,640) MIP(3),MIP(6),MIP(5),MIP(4),
     &                   NEVAL,NEVEC,MAXLAN,SHIFT
      ENDIF
C
      RETURN
C ------------------------------------------------ FORMATS
C
  600 FORMAT('1'/5X,43('*')//
     &           5X,'Print from  S A M  library routine SPRLAN'/
     &           5X,'(General eigenproblem by truncated Lanczos)'//
     &           5X,43('*')//)
  610 FORMAT(5X,'Matrix dimension     :',I9/
     &       5X,'Stored elements in A :',I9)
  611 FORMAT(5X,'Matrix dimension     :',I9)
  620 FORMAT(5X,'Form of B-matrix     :  sparse')
  630 FORMAT(5X,'Form of B-matrix     :  diagonal')
  640 FORMAT(/5X,'Start vector code        :',I5/
     &        5X,'Factorization indicator  :',I5/
     &        5X,'Parameter KEX            :',I5/
     &        5X,'Reorthogonalization code :',I5/
     &        5X,'Requested eigenvalues    :',I5/
     &        5X,'Requested eigenvectors   :',I5/
     &        5X,'Max. number of steps     :',I5/
     &        5X,'Shift value              :',1PE11.3)
C
      END
