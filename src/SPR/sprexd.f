C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPREXD (MSPAR,MTREES,MSIFA,SM,DIAG,LPU,IERR)

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           MSPAR(*), MTREES(*), MSIFA(*), LPU, IERR
      DOUBLE PRECISION  SM(*), DIAG(*)

C ======================================================================
C  S A M  library routine :  SPREXD                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To extract the diagonal from the factor matrix stored in SM.
C
C     Created   : Mar. 23, 2004 (kmo)
C     Revisions : Mnt. xx, 200x (kmo)
C
C     MSPAR - INTEGER(*)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Not changed.
C     SM - DOUBLE PRECISION(NSM)
C      Entry : As output from SPRADM.
C      Exit  : The L and D factors, A is destroyed.
C     LPU - INTEGER
C      Entry : Output unit.
C      Exit  : Not changed.
C     IERR - INTEGER
C      Entry : Not defined
C      Exit  : Is set to zero if a successful call, and to -1 otherwise
C
C     Working arrays
C     --------------
C
C     None
C
C     Procedures
C     ----------
C
C     SPRXD1
C
C     Intrinsic
C     ---------
C
C     None
C
C     Include Blocks
C     --------------
C
C     None
C
C     Common Blocks
C     -------------
C
C     None
C
C     ------------------------------------------------------------------

      INTEGER   LNZ,NEQ,NSUPER,NZEROL,XLINDX,XLNZ,XSUPER

      EXTERNAL  SPRXD1

C     ==================================================================

C     -------------------
C     GET INFO FROM MSPAR
C     -------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         IERR = -1
         CALL SPRER1 ( 30, 'SPREXD', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

      IERR   = 0
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NZEROL = MSPAR(16)

C     ------------------
C     SET ARRAY POINTERS
C     ------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      XLNZ   = MSPAR(48)
      LNZ    = 1 + NEQ

C     -------------------------------------------------------
C     EXTRACT D FROM L BY A CALL TO SPRXD1 AND STORE IN DIAG.
C     -------------------------------------------------------
      CALL SPRXD1 ( NSUPER, NEQ, NZEROL, MTREES(XSUPER), MTREES(XLINDX),
     $              MSIFA(XLNZ), SM(LNZ), DIAG )

      RETURN
      END


      SUBROUTINE SPRXD1 (NSUPER,NEQ,NZEROL,XSUPER,XLINDX,XLNZ,LNZ,DIAG)

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER,NEQ,NZEROL
      INTEGER           XSUPER(NSUPER+1),XLINDX(NSUPER+1),XLNZ(NSUPER+1)
      DOUBLE PRECISION  LNZ(NZEROL),DIAG(NEQ)

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRXD1
C
C --- TO EXTRACT THE DIAGONAL FROM THE FACTOR MATRIX STORED IN LNZ.
C
C
C     CREATED   : MAR. 23, 2004 (KMO)
C     REVISIONS : MNT. XX, 200X (KMO)
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
C     NEQ    : INTEGER
C              NUMBER OF EQUATIONS.
C     NZEROL : INTEGER
C              SIZE OF THE FACTOR MATRIX INCLUDING THE DIAGONAL.
C     XSUPER : INTEGER(NSUPER+1)
C              THE FUNDAMENTAL SUPERNODE PARTITION.
C     XLINDX : INTEGER(NSUPER+1)
C              THE NODE BLOCK NON-ZERO STRUCTURE OF L.
C     XLNZ   : INTEGER(NSUPER+1)
C              POINTER TO START OF SUBMATRICES IN LNZ FOR EACH
C              SUPERNODE.
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              THE CHOLESKY SUPERNODE SUBMATRICES OF L STORED CON-
C              SEQUTIVELY ON DENSE FORMAT INCLUDED THE DIAGONAL.
C
C     ON EXIT
C     -------
C
C     DIAG   : DOUBLE PRECISION(NEQ)
C              THE DIAGONAL OF L.
C
C     WORKING ARRAYS
C     --------------
C
C     NONE
C
C     PROCEDURES
C     ----------
C
C     NONE
C
C     INTRINSIC
C     ---------
C
C     NONE
C
C     ------------------------------------------------------------------

      INTEGER   CHDLEN,IPJS,J,JS

      DO 200 JS = 1, NSUPER
         CHDLEN = XLINDX(JS+1) - XLINDX(JS)
         IPJS = XLNZ(JS)
         DO 100 J = XSUPER(JS), XSUPER(JS+1)-1
            DIAG(J) = LNZ(IPJS)
            IPJS = IPJS + CHDLEN
            CHDLEN = CHDLEN - 1
  100    CONTINUE
  200 CONTINUE

      RETURN
      END
