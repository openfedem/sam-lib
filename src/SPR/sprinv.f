C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRINV
     $( MSPAR , MTREES, MSIFA , SV    , SEV   , SINDEX )

      IMPLICIT NONE


      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), SINDEX(*)
      DOUBLE PRECISION
     $  SV(*), SEV(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To write back a parent solution to SPR data structure for a child.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*
*     MSPAR - INTEGER(*)
*      Entry : As output from SPRSMB.
*      Exit  : Not changed.
*     MTREES - INTEGER(NTREES)
*      Entry : As output from SPRSMB.
*      Exit  : Not changed.
*     MSIFA - INTEGER(NMSIFA)
*      Entry : As output from SPRSMB.
*      Exit  : Not changed.
*     SV - DOUBLE PRECISION(NSM)
*      Entry : Not defined.
*      Exit  : The solution vector of the child superelement updated
*              with the parent solution for retracking.
*     SEV - DOUBLE PRECISION(NEQ2)
*      Entry : The full superelement solution extracted from the parent
*              to retrack the child internal superelement degrees of
*              freedom.
*              NEQ2 = MSPAR(54)
*      Exit  : Not changed
*
*     Working arrays
*     --------------
*
*     SINDEX - INTEGER(NEQ)
*      Entry : Not defined.
*      Exit  : Need not be saved. Contains for each equation in the
*              system that is retained, the pointer to the location
*              in the superelement matrix in SEV.
*
* --- ------------------------------------------------------------------

      INTEGER
     $  NEQ   , NSUPER, NOFSUB, NEQ2

      INTEGER
     $  XLINDX, LINDX , XSUPER, SUPSUP, XSIBL , SIBL

* --- MSPAR
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)
      NEQ2   = MSPAR(54)

* --- MTREES
      XSUPER = MSPAR(45)
      XSIBL  = XSUPER + NSUPER + 1
      SIBL   = XSIBL  + NSUPER + 1
      XLINDX = MSPAR(46)
      SUPSUP = XLINDX + NSUPER + 1

* --- MSIFA
      LINDX  = MSPAR(47)

* --- Extract the superelement load vector.
      CALL SPRIN1
     $( NEQ,NSUPER,NEQ2,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $  MSIFA(LINDX),SV,SEV,SINDEX )

      RETURN
      END


      SUBROUTINE SPRIN1
     $( NEQ   , NSUPER, NEQ2  , XSUPER, XLINDX, SUPSUP, LINDX , SV    ,
     $  SEV   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  NEQ, NSUPER, NEQ2,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  SINDEX(NEQ)

      DOUBLE PRECISION
     $  SV(NEQ), SEV(NEQ2)

* --- ------------------------------------------------------------------

      INTEGER
     $  I, ISTRT, ISTOP, JS, ROWIND, SUPSZE

* --- ------------------------------------------------------------------

      DO I = 1, NEQ
         SINDEX(I) = 0
      END DO

      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1
            DO I = ISTRT, ISTOP
               ROWIND = LINDX(I)
               SINDEX(ROWIND) = ROWIND
            END DO
         ENDIF
      END DO

      ROWIND = 0
      DO I = 1, NEQ
         IF ( SINDEX(I).GT.0 ) THEN
            ROWIND = ROWIND + 1
            SINDEX(I) = ROWIND
*ACD             SV(I) = 0.0D0
         ENDIF
      END DO

      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            I = XLINDX(JS)
            SUPSZE = XSUPER(JS+1)-XSUPER(JS)
            CALL SPRIN2
     $    ( NEQ, NEQ2, SUPSZE, LINDX(I), SV, SEV, SINDEX )
         ENDIF
      END DO

      RETURN
      END


      SUBROUTINE SPRIN2
     $( NEQ   , LDSEV , M     , LINDX , SV    , SEV   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  NEQ   , LDSEV , M     , LINDX(M), SINDEX(NEQ)

      DOUBLE PRECISION
     $  SV(NEQ), SEV(LDSEV)

      INTEGER
     $  I, II, IND

      DO I = 1, M
         IND = LINDX(I)
         II = SINDEX(IND)
*ACD          SV(IND) = SV(IND) + SEV(II)
         SV(IND) = SEV(II)
      END DO

      RETURN
      END
