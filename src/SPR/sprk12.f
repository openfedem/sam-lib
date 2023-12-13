C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRK12
     $( MSPAR , MTREES, MSIFA , SM    , K12   , RINDEX )

      IMPLICIT NONE

      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), RINDEX(*)

      DOUBLE PRECISION
     $  SM(*), K12(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To extract the lower triangular matrix associated with the
*     external dofs from SM. If this is done prior to factorization,
*     but after assembly the results is the coefficient matrix A12.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 Mar. 14, 2003 (kmo)
*                 Added internal permutation array PERE2I.
*                 Removed some unused local variables and argument LPU.
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
*     SM - DOUBLE PRECISION(NSM)
*      Entry : As output from SPRFC1.
*      Exit  : Not changed.
*     K12 - DOUBLE PRECISION(NEQ1,NEQ2)
*      Entry : Not defined.
*      Exit  : The matrix associated with the external dofs, extracted
*              from the substructure matrix.
*
*     Working arrays
*     --------------
*
*     RINDEX - INTEGER(NEQ)
*      Entry : Not defined.
*      Exit  : Need not be saved. Contains for each equation in the
*              system that is external the row index, ie the equation.
*
* --- ------------------------------------------------------------------

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ   , NEQ1  , NEQ2  , NSUPER

      INTEGER
     $  XLINDX, PERE2I, LINDX , XLNZ  , XSUPER, SUPSUP

* --- MSPAR
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NEQ2   = MSPAR(54)
      NEQ1   = NEQ - NEQ2
      LIPERM = MSPAR(58) .GT. 0

* --- MTREES
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      SUPSUP = XLINDX + NSUPER + 1

* --- MSIFA
      PERE2I = 1 + NEQ
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

* --- Extract K12.
      CALL SPRKC1
     $( NEQ,NEQ1,NEQ2,NSUPER,MTREES(XSUPER),MTREES(XLINDX),
     $  MTREES(SUPSUP),MSIFA(LINDX),MSIFA(XLNZ),SM,K12,RINDEX,
     $  LIPERM,MSIFA(PERE2I) )

      RETURN
      END


      SUBROUTINE SPRKC1
     $( NEQ   , NEQ1  , NEQ2  , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX ,
     $  XLNZ  , SM    , K12   , RINDEX, LIPERM, PERE2I )

      IMPLICIT NONE

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ, NEQ1, NEQ2, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ), PERE2I(NEQ)

      DOUBLE PRECISION
     $  SM(*), K12(NEQ1,NEQ2)

* --- ------------------------------------------------------------------

      INTEGER
     $  COLIND, I, IPSM, ISTRT, ISTOP, J, JS, JSTRT, JSTOP, ROWIND

      DOUBLE PRECISION
     $  ZERO
      PARAMETER
     $( ZERO = 0.0D+00 )

* --- ------------------------------------------------------------------

      DO J = 1, NEQ2
         DO I = 1, NEQ1
            K12(I,J) = ZERO
         END DO
      END DO

* --- Compute RINDEX.
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER

         IF ( SUPSUP(JS).GT.0 ) THEN

* --------- This supernode is in the retained
*           set and we mark the indices.
            DO I = XSUPER(JS), XSUPER(JS+1)-1
               RINDEX(I) = I
            END DO

         ENDIF

      END DO

      COLIND = 0
      ROWIND = 0
      DO J = 1, NEQ
         if ( LIPERM ) then
            I = PERE2I(J)
         else
            I = J
         end if

         IF ( RINDEX(I).GT.0 ) THEN

* --------- This index belongs to retained set and is a
*           row index in K21.
            ROWIND = ROWIND + 1
            RINDEX(I) = ROWIND

         ELSE

* --------- This index belongs to internal set and is a
*           column index in K21.
            COLIND = COLIND + 1
            RINDEX(I) = -COLIND

         ENDIF

      END DO

      DO JS = 1, NSUPER

         IF ( SUPSUP(JS).LE.0 ) THEN

* --------- This supernode is in the internal set
*           and we need to scatter its coefficients
*           that belong to K21 to K12.
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Next position in SM.
            IPSM = XLNZ(JS)

* --------- Insert the coefficient set for the current supernode.
            DO J = JSTRT, JSTOP

* ------------ Get the column of K21. We know that it is negative
*              and we need to negate.
               COLIND = -RINDEX(J)

* ------------ Now, traverse down the column of J while picking
*              the row indices that belong to K21.
               DO I = ISTRT, ISTOP

* --------------- Get the row index.
                  ROWIND = RINDEX(LINDX(I))
                  IF ( ROWIND.GT.0 ) THEN

* ------------------ This is a row index that belongs to K21.
*                    That is we have a column of K12 and COLIND
*                    is the corresponding row.
                     K12(COLIND,ROWIND) = SM(IPSM)
                  ENDIF

                  IPSM = IPSM + 1

               END DO

* ------------ Update index pointer.
               ISTRT = ISTRT + 1

            END DO

         ENDIF

      END DO

      RETURN
      END
