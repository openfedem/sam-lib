C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRK22
     $( MSPAR , MTREES, MSIFA , SM    , IFLAG , XK22  , K22   , RINDEX)

      IMPLICIT NONE

      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), IFLAG, XK22(*), RINDEX(*)

      DOUBLE PRECISION
     $  SM(*), K22(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To extract the lower triangular matrix associated with the
*     external dofs from SM. If this is done prior to factorization,
*     but after assembly the results is the coefficient matrix A22.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 Mar. 14, 2003 (kmo)
*                 Added internal permutation array PERE2I.
*                 Removed some unused local variables and argument LPU.
*                 May. 06, 2003 (kmo & acd)
*                 Added subroutine SPRKB3 with correct use of PERE2I.
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
*     IFLAG - INTEGER
*      Entry : Option for storage of K22.
*              =  0 : Compressed storage.
*              =  1 : Full storage, lower triangle only
*              =  2 : Full storage, symmetric rectangular matrix
*              We have always a full K22 if MSPAR(53) = 1.
*     XK22 - INTEGER(SZER22)
*      Entry : Not defined.
*      Exit  : The index information of K22.
*              SZER22 = MSPAR(56)
* ***          Not changed when SPRKB3 is used.
*     K22 - DOUBLE PRECISION(SZEL22)
*      Entry : Not defined.
*      Exit  : The matrix associated with the external dofs, extracted
*              from the substructure matrix. The matrix is stored on
*              compressed form.
*              SZEL22 = MSPAR(55)
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
     $  NEQ   , NEQ2  , NSUPER, NOFSUB

      INTEGER
     $  XLINDX, PERE2I, LINDX , XLNZ  , XSUPER, SUPSUP

      INTEGER
     $  NSUP22, XSUP22, XSUB22, XLNZ22, SUB22

* --- MSPAR
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)
      NSUP22 = MSPAR(53)
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

* --- XK22
      XSUP22 = 1
      XSUB22 = XSUP22 + NSUP22 + 1
      XLNZ22 = XSUB22 + NSUP22 + 1
      SUB22  = XLNZ22 + NSUP22 + 1

* --- Extract K22 (XK22)
      IF ( IFLAG.EQ.0 ) THEN
         CALL SPRKB1
     $   ( NEQ,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $     MSIFA(LINDX),MSIFA(XLNZ),NSUP22,XK22(XSUP22),XK22(XSUB22),
     $     XK22(XLNZ22),XK22(SUB22),SM,K22,RINDEX )
      ELSEIF ( IFLAG.EQ.1 ) THEN
         CALL SPRKB2
     $   ( NEQ,NEQ2,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $     MSIFA(LINDX),MSIFA(XLNZ),NSUP22,XK22(XSUP22),XK22(XSUB22),
     $     XK22(XLNZ22),XK22(SUB22),SM,K22,RINDEX )
      ELSE
         CALL SPRKB3
     $   ( NEQ,NEQ2,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $     MSIFA(LINDX),MSIFA(XLNZ),SM,K22,RINDEX,LIPERM,MSIFA(PERE2I) )
      ENDIF

      RETURN
      END


      SUBROUTINE SPRKB1
     $( NEQ   , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX , XLNZ  , NSUP22,
     $  XSUP22, XSUB22, XLNZ22, SUB22 , SM    , K22   , RINDEX  )

      IMPLICIT NONE

      INTEGER
     $  NEQ, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ)

      INTEGER
     $  NSUP22, XSUP22(NSUP22+1), XSUB22(NSUP22+1), XLNZ22(NSUP22+1),
     $  SUB22(*)

      DOUBLE PRECISION
     $  SM(*), K22(*)

* --- ------------------------------------------------------------------

      INTEGER
     $  I,IPK22,IPSM,IPXK22,IPXSUP,ISTRT,ISTOP,J,JS,JSTRT,JSTOP,ROWIND,
     $  SUP22T

* --- ------------------------------------------------------------------

* --- Compute RINDEX
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1
            DO I = ISTRT, ISTOP
               ROWIND = LINDX(I)
               RINDEX(ROWIND) = ROWIND
            END DO
         ENDIF
      END DO
      ROWIND = 0
      DO I = 1, NEQ
         IF ( RINDEX(I).GT.0 ) THEN
            ROWIND = ROWIND + 1
            RINDEX(I) = ROWIND
         ENDIF
      END DO

* --- Start NSUP22
      SUP22T = 0

* --- Start XSUP22
      IPXSUP = 1
      XSUP22(SUP22T+1) = IPXSUP

* --- Start XLNZ22
      IPK22 = 1
      XLNZ22(SUP22T+1) = IPK22

* --- Start XSUB22
      IPXK22 = 1
      XSUB22(SUP22T+1) = IPXK22

      DO JS = 1, NSUPER

         IF ( SUPSUP(JS).GT.0 ) THEN

            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1

* --------- Insert the row index set for current supernode
            DO I = ISTRT, ISTOP
               ROWIND = LINDX(I)
               IF ( RINDEX(ROWIND).GT.0 ) THEN
                  SUB22(IPXK22) = RINDEX(ROWIND)
                  IPXK22 = IPXK22 + 1
               ENDIF
            END DO

            IPSM = XLNZ(JS)
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Insert the coefficient set for current supernode
            DO J = JSTRT, JSTOP
               DO I = ISTRT, ISTOP
                  ROWIND = LINDX(I)
                  IF ( RINDEX(ROWIND).GT.0 ) THEN
                     K22(IPK22) = SM(IPSM)
                     IPK22 = IPK22 + 1
                  ENDIF
                  IPSM = IPSM + 1
               END DO
               ISTRT = ISTRT + 1
               IPXSUP = IPXSUP + 1
            END DO

* --------- Update the number of supernodes in K22
            SUP22T = SUP22T + 1

* --------- Update XSUP22, XSUB22 and XLNZ22
            XSUP22(SUP22T+1) = IPXSUP
            XSUB22(SUP22T+1) = IPXK22
            XLNZ22(SUP22T+1) = IPK22

         ENDIF

      END DO

      RETURN
      END


      SUBROUTINE SPRKB2
     $( NEQ   , NEQ2  , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX , XLNZ  ,
     $  NSUP22, XSUP22, XSUB22, XLNZ22, SUB22 , SM    , K22   , RINDEX )

      IMPLICIT NONE

      INTEGER
     $  NEQ, NEQ2, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ)

      INTEGER
     $  NSUP22, XSUP22(NSUP22+1), XSUB22(NSUP22+1), XLNZ22(NSUP22+1),
     $  SUB22(*)

      DOUBLE PRECISION
     $  SM(*), K22(*)

* --- ------------------------------------------------------------------

      INTEGER
     $  I,IPK22,IPSM,ISTRT,ISTOP,J,JS,JSTRT,JSTOP,NCOL,NROW,ROWIND

      DOUBLE PRECISION
     $  ZERO
      PARAMETER
     $( ZERO = 0.0D+00 )

* --- ------------------------------------------------------------------

      DO I = 1, NEQ2*(NEQ2+1)/2
         K22(I) = ZERO
      END DO

* --- Compute RINDEX and SUB22
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1
            DO I = ISTRT, ISTOP
               ROWIND = LINDX(I)
               RINDEX(ROWIND) = ROWIND
            END DO
         ENDIF
      END DO
      ROWIND = 0
      DO I = 1, NEQ
         IF ( RINDEX(I).GT.0 ) THEN
            ROWIND = ROWIND + 1
            RINDEX(I) = ROWIND
            SUB22(ROWIND) = ROWIND
         ENDIF
      END DO

* --- Start XSUP22
      XSUP22(1) = 1

* --- Start XLNZ22
      IPK22 = 1
      XLNZ22(1) = IPK22

* --- Start XSUB22
      XSUB22(1) = 1

      NCOL = 0
      DO JS = 1, NSUPER

         IF ( SUPSUP(JS).GT.0 ) THEN

            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1

            IPSM  = XLNZ(JS)
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Insert the coefficient set for current supernode
            DO J = JSTRT, JSTOP
               DO I = ISTRT, ISTOP
                  ROWIND = LINDX(I)
                  NROW = RINDEX(ROWIND)
                  IF ( NROW.GT.0 ) THEN
                     IPK22 = NCOL*NEQ2-(NCOL*(NCOL-1)/2)+NROW-NCOL
                     K22(IPK22) = SM(IPSM)
                  ENDIF
                  IPSM = IPSM + 1
               END DO
               ISTRT = ISTRT + 1
               NCOL = NCOL + 1
            END DO

         ENDIF

      END DO

* --- Update XSUP22, XSUB22 and XLNZ22
      XSUP22(2) = XSUP22(1) + NEQ2
      XSUB22(2) = XSUB22(1) + NEQ2
      XLNZ22(2) = XLNZ22(1) + NEQ2*(NEQ2+1)/2

      RETURN
      END


      SUBROUTINE SPRKB3
     $( NEQ   , NEQ2  , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX , XLNZ  ,
     $  SM    , K22   , RINDEX, LIPERM, PERE2I  )

      IMPLICIT NONE

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ, NEQ2, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ), PERE2I(NEQ)

      DOUBLE PRECISION
     $  SM(*), K22(NEQ2,NEQ2)

* --- ------------------------------------------------------------------

      INTEGER
     $  I,IPSM,ISTRT,ISTOP,J,JS,JSTRT,JSTOP,NCOL,NROW

      DOUBLE PRECISION
     $  ZERO
      PARAMETER
     $( ZERO = 0.0D+00 )

* --- ------------------------------------------------------------------

      DO J = 1, NEQ2
         DO I = 1, NEQ2
            K22(I,J) = ZERO
         END DO
      END DO

* --- Compute RINDEX
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN
            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1
            DO I = ISTRT, ISTOP
               NROW = LINDX(I)
               RINDEX(NROW) = NROW
            END DO
         ENDIF
      END DO
      NROW = 0
      DO J = 1, NEQ
         if ( LIPERM ) then
            I = PERE2I(J)
         else
            I = J
         end if
         IF ( RINDEX(I).GT.0 ) THEN
            NROW = NROW + 1
            RINDEX(I) = NROW
         ENDIF
      END DO

      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).GT.0 ) THEN

            ISTRT = XLINDX(JS)+1
            ISTOP = XLINDX(JS+1)-1
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Insert the coefficient set for current supernode
            IPSM = XLNZ(JS)
            DO J = JSTRT, JSTOP
               NCOL = RINDEX(J)
               K22(NCOL,NCOL) = SM(IPSM)
               IPSM = IPSM + 1
               DO I = ISTRT, ISTOP
                  NROW = RINDEX(LINDX(I))
                  IF ( NROW.GT.0 ) THEN
                     K22(NROW,NCOL) = SM(IPSM)
                     K22(NCOL,NROW) = SM(IPSM)
                  ENDIF
                  IPSM = IPSM + 1
               END DO
               ISTRT = ISTRT + 1
            END DO

         ENDIF
      END DO

      RETURN
      END
