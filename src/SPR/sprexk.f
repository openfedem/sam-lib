C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPREXK
     $( MSPAR , MTREES, MSIFA , SM    , SEM   , SINDEX )

      IMPLICIT NONE

      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), SINDEX(*)

      DOUBLE PRECISION
     $  SM(*), SEM(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To extract a superelement stiffness matrix from SM.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 May. 12, 2003 (kmo)
*                 Added internal permutation array PERE2I.
*                 Use SPRKB3 instead of SPREX1 which has been removed.
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
*     SEM - DOUBLE PRECISION(NEQ2*NEQ2)
*      Entry : Not defined.
*      Exit  : The full, but symmetric superelement matrix, extracted
*              from the partially factorized substructure matrix.
*              NEQ2 = MSPAR(54)
*
*     Working arrays
*     --------------
*
*     SINDEX - INTEGER(NEQ)
*      Entry : Not defined.
*      Exit  : Need not be saved. Contains for each equation in the
*              system that is retained, the pointer to the location
*              in the superelement matrix in SEM.
*
* --- ------------------------------------------------------------------

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ   , NEQ2  , NSUPER, NOFSUB

      INTEGER
     $  XLINDX, PERE2I, LINDX , XLNZ  , XSUPER, SUPSUP

* --- MSPAR
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)
      NEQ2   = MSPAR(54)
      LIPERM = MSPAR(58) .GT. 0

* --- MTREES
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      SUPSUP = XLINDX + NSUPER + 1

* --- MSIFA
      PERE2I = 1 + NEQ
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

* --- Extract the superelement stiffness matrix.
      CALL SPRKB3
     $( NEQ,NEQ2,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $  MSIFA(LINDX),MSIFA(XLNZ),SM,SEM,SINDEX,LIPERM,MSIFA(PERE2I) )

      RETURN
      END
