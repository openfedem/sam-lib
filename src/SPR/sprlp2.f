C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRLP2 (EVL,EVERR,TOL,MANCUR,MOP,NLV,LPU)
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRLP2                 GROUP 9 / PRIVATE
C
C     TASK :  To print key information on completion of eigensolution
C             in LANCZ2, including eigenvalues and their relative errors
C
C
C     ROUTINES CALLED/REFERENCED :  SPRLP3     (SAM-9)
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
      INTEGER           LPU,NLV,       MANCUR(NLV),MOP(10)
      DOUBLE PRECISION  EVL(NLV),EVERR(NLV),TOL(2)
C
C ----------------------------------------------------------------------
      IF (LPU.GT.0) THEN
         WRITE (LPU,600) MOP(1),MOP(2),MOP(3),MOP(4),MOP(5),
     &                   MOP(6),MOP(7),MOP(8),MOP(9),MOP(10)
         WRITE (LPU,610) TOL(1),TOL(2)
         CALL SPRLP3 (EVL,EVL,EVL,EVERR,MANCUR,TOL(1),NLV,MOP(2),LPU,0)
      ENDIF
C
      RETURN
C ----------------------------------------------------------------------
C
  600 FORMAT(////5X,14('*')/5X,'Leaving SPRLAN'/5X,14('*')//
     &           5X,'No. of "good" eigenvalues             :',I5/
     &           5X,'No. of acceptable eigenvalues         :',I5/
     &           5X,'No. of Lanczos steps (vectors)        :',I5/
     &           5X,'No. of eigenvalues smaller than shift :',I5/
     &           5X,'No. of eigenvectors returned          :',I5/
     &           5X,'No. of zero elements in diag.(B)      :',I5/
     &           5X,'No. of tridiag. eigv.problems solved  :',I5/
     &           5X,'No. of Lanczos vectors for which'/
     &           5X,'reorthogonalization was carried out   :',I5/
     &           5X,'No. of vector dot products performed'/
     &           5X,'during reorthogonalization            :',I5/
     &           5X,'No. of restarts                       :',I5)
  610 FORMAT(/   5X,'Relative eigenvalue tolerance    :',1PE11.3 /
     &           5X,'Diagonal decay tolerance         :',1PE11.3 )
C
      END
