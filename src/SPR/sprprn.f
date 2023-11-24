C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRPRN (L     , MSPAR , MTREES, MSIFA , LPU   , IERR  )

C     ------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION   L(*)
      INTEGER            MSPAR(*), MTREES(*), MSIFA(*), LPU, IERR

C ======================================================================
C  S A M  library routine :  SPRPRN                   GROUP 9 / PUBLIC
C ======================================================================
C
C     Purpose
C     -------
C
C --- To print out the compressed stored symmetric matrix L.
C
C     Created   : Mar. 11, 2003 (kmo)
C
C     L - DOUBLE PRECISION(*)
C      Entry : A symmetric matrix where the lower triangle is stored
C              in compressed format.
C      Exit  : Unchanged.
C     MSPAR - INTEGER(*)
C      Entry : Matrix of sparse parameters as output from SPRSMB.
C      Exit  : Unchanged.
C     MTREES - INTEGER(NTREES)
C      Entry : As output from SPRSMB.
C      Exit  : Unchanged.
C     MSIFA - INTEGER(NMSIFA)
C      Entry : As output from SPRSMB.
C      Exit  : Unchanged.
C     LPU - INTEGER
C      Entry : Output for error messages.
C      Exit  : Unchanged.
C     IERR - INTEGER
C      Entry : Need not be set.
C      Exit  : Error flag, it is set to zero in case of a normal return.
C
C     Working arrays
C     --------------
C
C     None
C
C     Procedures
C     ----------
C
C     SPRPR1
C     SPRER1
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

      INTEGER           LINDX , NELLIN, NOFSUB, NSUPER,
     $                  XLINDX, XLNZ  , XSUPER

      PARAMETER       ( NELLIN = 10 )

      EXTERNAL          SPRPR1, SPRER1

C     ==================================================================

C     --------------------
C     GET INFO FROM MSPAR.
C     --------------------
      IF ( MSPAR(1) .LT. 0 ) THEN
         MSPAR(1) = -11
         IERR = -1
         CALL SPRER1 ( 30, 'SPRPRN', 0, 0, 0, LPU, IERR )
         RETURN
      ENDIF

      IERR   = 0
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)

C     -----------------------
C     SET POINTERS TO MTREES.
C     -----------------------
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)

C     ----------------------
C     SET POINTERS TO MSIFA.
C     ----------------------
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

C     -----------------------------------
C     PRINT L TO LPU BY A CALL TO SPRPR1.
C     -----------------------------------
      CALL SPRPR1 ( LPU, NELLIN, NSUPER, NOFSUB, MTREES(XSUPER),
     $              MSIFA(XLNZ), MTREES(XLINDX), MSIFA(LINDX), L )

C     ==================================================================

      RETURN

C     ------------------------------------------------------------------
C     END OF MODULE SPRPRN
      END


      SUBROUTINE SPRPR1 (LPU   , NELLIN, NSUPER, NOFSUB, XSUPER, XLNZ  ,
     $                   XLINDX, LINDX , LNZ   )

      IMPLICIT NONE

      INTEGER            LPU   , NELLIN, NSUPER, NOFSUB
      INTEGER            XSUPER(NSUPER+1)      , XLNZ(NSUPER+1)        ,
     $                   XLINDX(NSUPER+1)      , LINDX(NOFSUB)
      DOUBLE PRECISION   LNZ(*)

C     -------------------------------------------------------
C --- TO PRINT OUT A LOWER SYMMETRIC MATRIX COLUMN BY COLUMN.
C     -------------------------------------------------------
      INTEGER            I     , II    , IPLNZ , ISTRT , ISTOP ,
     $                   J     , JS    , JSTRT , LENLN , NLINE

      DO 400 JS = 1, NSUPER
         ISTRT  = XLINDX(JS)
         ISTOP  = XLINDX(JS+1) - 1
         IPLNZ  = XLNZ(JS)
         DO 300 J = XSUPER(JS), XSUPER(JS+1)-1
            NLINE = 1 + (ISTOP-ISTRT)/NELLIN
            JSTRT = ISTRT
            WRITE(LPU,600) J
            DO 200 II = 1, NLINE
               LENLN = MIN(NELLIN,ISTOP-JSTRT+1)
               WRITE(LPU,610) (LINDX(JSTRT+I),I=0,LENLN-1)
               WRITE(LPU,620) (LNZ(IPLNZ+I),I=0,LENLN-1)
               JSTRT = JSTRT + LENLN
               IPLNZ = IPLNZ + LENLN
  200       CONTINUE
            ISTRT = ISTRT + 1
  300    CONTINUE
  400 CONTINUE

C     ==================================================================
      RETURN

  600 FORMAT(/' +++++ Equation',I8,' +++++')
  610 FORMAT(I9,19I13)
  620 FORMAT(1P20E13.5)

C     ------------------------------------------------------------------
C     END OF MODULE SPRPR1
      END
