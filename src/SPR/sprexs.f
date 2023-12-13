C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPREXS
     $( MSPAR , MTREES, MSIFA , NRHS  , SV    , SEV   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), NRHS, SINDEX(*)

      DOUBLE PRECISION
     $  SV(*), SEV(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To extract a superelement load vector from SV.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 Apr. 07, 2008 (kmo)
*                 Added NRHS as argument. Extract all RHSs in one go.
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
*     NRHS - INTEGER
*      Entry : Number of right-hand sides stored in SV and SEV.
*      Exit  : Not changed.
*     SV - DOUBLE PRECISION(NEQ,NRHS)
*      Entry : As output from SPRFS1.
*      Exit  : Not changed.
*     SEV - DOUBLE PRECISION(NEQ2,NRHS)
*      Entry : Not defined.
*      Exit  : The full superelement load extracted from the partially
*              reduced substructure vector.
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
      CALL SPREX3
     $( NSUPER,NEQ,NEQ2,NRHS,MTREES(XSUPER),MTREES(XLINDX),
     $  MTREES(SUPSUP),MSIFA(LINDX),SV,SEV,SINDEX )

      RETURN
      END


      SUBROUTINE SPREX3
     $( NSUPER, NEQ   , NEQ2  , NRHS  , XSUPER, XLINDX, SUPSUP, LINDX ,
     $  SV    , SEV   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  NSUPER, NEQ, NEQ2, NRHS,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  SINDEX(NEQ)

      DOUBLE PRECISION
     $  SV(NEQ,NRHS), SEV(NEQ2,NRHS)

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
         ENDIF
      END DO

      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            I = XLINDX(JS)
            SUPSZE = XSUPER(JS+1)-XSUPER(JS)
            CALL SPREX4
     $    ( NEQ, NEQ2, NRHS, SUPSZE, LINDX(I), SV, SEV, SINDEX )
         ENDIF
      END DO

      RETURN
      END


      SUBROUTINE SPREX4
     $( NEQ   , LDSEV , NRHS  , M     , LINDX , SV    , SEV   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  NEQ   , LDSEV , NRHS  , M     , LINDX(M), SINDEX(NEQ)

      DOUBLE PRECISION
     $  SV(NEQ,NRHS), SEV(LDSEV,NRHS)

      INTEGER
     $  I, J, IND

      DO I = 1, M
         IND = LINDX(I)
         DO J = 1, NRHS
            SEV(SINDEX(IND),J) = SV(IND,J)
         END DO
      END DO

      RETURN
      END
