C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRDAD (MSPAR,MTREES,MSIFA,SM,SIGMA,LPU,IERR)

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           MSPAR(*), MTREES(*), MSIFA(*), LPU, IERR
      DOUBLE PRECISION  SM(*), SIGMA

C ======================================================================
C  S A M  library routine :  SPRDAD                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To add a scalar to the diagonal of the factor matrix stored in SM.
C
C     Created   : May. 25, 2009 (kmo)
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
C      Exit  : The same matrix, but with SIGMA added to the diagonal.
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
C     SPRDA1
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

      INTEGER   LNZ,NSUPER,NZEROL,XLINDX,XLNZ,XSUPER

      EXTERNAL  SPRDA1

C     ==================================================================

C     -------------------
C     GET INFO FROM MSPAR
C     -------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         IERR = -1
         CALL SPRER1 ( 30, 'SPRDAD', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

      IERR   = 0
      NSUPER = MSPAR(11)
      NZEROL = MSPAR(16)

C     ------------------
C     SET ARRAY POINTERS
C     ------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      XLNZ   = MSPAR(48)
      LNZ    = MSPAR( 8) + 1

C     ---------------------------------------------------
C     ADD SIGMA TO THE DIAGONAL OF L BY A CALL TO SPRDA1.
C     ---------------------------------------------------
      CALL SPRDA1 ( NSUPER, NZEROL, MTREES(XSUPER), MTREES(XLINDX),
     $              MSIFA(XLNZ), SM(LNZ), SIGMA )

      RETURN
      END


      SUBROUTINE SPRDA1 (NSUPER,NZEROL,XSUPER,XLINDX,XLNZ,LNZ,SIGMA)

C     ------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER           NSUPER,NZEROL
      INTEGER           XSUPER(NSUPER+1),XLINDX(NSUPER+1),XLNZ(NSUPER+1)
      DOUBLE PRECISION  LNZ(NZEROL),SIGMA

C     ------------------------------------------------------------------
C
C     PURPOSE
C     -------
C
C     SPRDA1
C
C --- TO ADD A SCALAR TO THE DIAGONAL OF THE FACTOR MATRIX STORED IN LNZ.
C
C
C     CREATED   : MAY. 25, 2009 (KMO)
C     REVISIONS : MNT. XX, 200X (KMO)
C
C
C     ON ENTRY
C     --------
C
C     NSUPER : INTEGER
C              NUMBER OF SUPERNODES.
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
C     SIGMA  : DOUBLE PRECISION
C              SCALAR VALUE TO ADD TO ALL DIAGONAL TERMS OF L
C
C     ON EXIT
C     -------
C
C     LNZ    : DOUBLE PRECISION(NZEROL)
C              THE CHOLESKY SUPERNODE SUBMATRICES OF L
C              WITH UPDATED THE DIAGONAL.
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
            LNZ(IPJS) = LNZ(IPJS) + SIGMA
            IPJS = IPJS + CHDLEN
            CHDLEN = CHDLEN - 1
  100    CONTINUE
  200 CONTINUE

      RETURN
      END
