C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRK11
     $( MSPAR , MTREES, MSIFA , SM    , IFLAG , XK11   , K11   , RINDEX)

      IMPLICIT NONE

      INTEGER
     $  MSPAR(*), MTREES(*), MSIFA(*), IFLAG, XK11(*), RINDEX(*)

      DOUBLE PRECISION
     $  SM(*), K11(*)

* --- ------------------------------------------------------------------
*
*     Purpose
*     -------
*
* --- To extract the lower triangular matrix associated with the
*     internal dofs from SM. If this is done prior to factorization,
*     but after assembly the results is the coefficient matrix A11.
*
*     Created   : Jul. 01, 1998 (acd)
*     Revisions : Feb. 27, 2003 (kmo)
*                 Get pointers from MSPAR instead of recalculating them.
*                 Apr. 11, 2011 (kmo)
*                 Removed some unused local variables and argument LPU.
*                 Added subroutine SPRKA3.
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
*      Entry : Option for storage of K11.
*              =  0 : Compressed storage.
*              =  2 : Full storage, symmetric rectangular matrix
*     XK11 - INTEGER(SZER11)
*      Entry : Not defined.
*      Exit  : The index information of K11.
*              SZER11 = MSPAR(60)
*     K11 - DOUBLE PRECISION(SZEL11)
*      Entry : Not defined.
*      Exit  : The matrix associated with the internal dofs, extracted
*              from the substructure matrix. The matrix is stored on
*              compressed form.
*              SZEL11 = MSPAR(59)
*
*     Working arrays
*     --------------
*
*     RINDEX - INTEGER(NEQ)
*      Entry : Not defined.
*      Exit  : Need not be saved. Contains for each equation in the
*              system that is internal the row index, ie the equation.
*
* --- ------------------------------------------------------------------

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ   , NEQ1  , NSUPER, NOFSUB

      INTEGER
     $  XLINDX, PERE2I, LINDX , XLNZ  , XSUPER, SUPSUP

      INTEGER
     $  NSUP11, SUPRET, XSUP11, XSUB11, XLNZ11, SUB11

* --- MSPAR
      NEQ    = MSPAR( 8)
      NSUPER = MSPAR(11)
      NOFSUB = MSPAR(15)
      SUPRET = MSPAR(53)
      NEQ1   = NEQ - MSPAR(54)
      LIPERM = MSPAR(58) .GT. 0

* --- MTREES
      XSUPER = MSPAR(45)
      XLINDX = MSPAR(46)
      SUPSUP = XLINDX + NSUPER + 1

* --- MSIFA
      PERE2I = 1 + NEQ
      LINDX  = MSPAR(47)
      XLNZ   = MSPAR(48)

* --- NSUP11 - the number of internal supernodes.
      NSUP11 = NSUPER - SUPRET

* --- XK11
      XSUP11 = 1
      XSUB11 = XSUP11 + NSUP11 + 1
      XLNZ11 = XSUB11 + NSUP11 + 1
      SUB11  = XLNZ11 + NSUP11 + 1

* --- Extract K11 (XK11).
      IF ( IFLAG.EQ.0 ) THEN
         CALL SPRKA1
     $   ( NEQ,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $     MSIFA(LINDX),MSIFA(XLNZ),NSUP11,XK11(XSUP11),XK11(XSUB11),
     $     XK11(XLNZ11),XK11(SUB11),SM,K11,RINDEX )
      ELSE
         CALL SPRKA3
     $   ( NEQ,NEQ1,NSUPER,MTREES(XSUPER),MTREES(XLINDX),MTREES(SUPSUP),
     $     MSIFA(LINDX),MSIFA(XLNZ),SM,K11,RINDEX,LIPERM,MSIFA(PERE2I) )
      ENDIF

      RETURN
      END


      SUBROUTINE SPRKA1
     $( NEQ   , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX , XLNZ  , NSUP11,
     $  XSUP11, XSUB11, XLNZ11, SUB11 , SM    , K11   , RINDEX  )

      IMPLICIT NONE

      INTEGER
     $  NEQ, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ)

      INTEGER
     $  NSUP11, XSUP11(NSUP11+1), XSUB11(NSUP11+1), XLNZ11(NSUP11+1),
     $  SUB11(*)

      DOUBLE PRECISION
     $  SM(*), K11(*)

* --- ------------------------------------------------------------------

      INTEGER
     $  I,IPK11,IPSM,IPXK11,IPXSUP,ISTRT,ISTOP,J,JS,JSTRT,JSTOP,
     $  ROWIND,SUP11T

* --- ------------------------------------------------------------------

* --- Compute RINDEX.
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).LE.0 ) THEN
            ISTRT = XSUPER(JS)
            ISTOP = XSUPER(JS+1)-1
            DO I = ISTRT, ISTOP
               RINDEX(I) = I
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

* --- Start NSUP11
      SUP11T = 0

* --- Start XSUP11
      IPXSUP = 1
      XSUP11(SUP11T+1) = IPXSUP

* --- Start XLNZ11
      IPK11 = 1
      XLNZ11(SUP11T+1) = IPK11

* --- Start XSUB11
      IPXK11 = 1
      XSUB11(SUP11T+1) = IPXK11

      DO JS = 1, NSUPER

         IF ( SUPSUP(JS).LE.0 ) THEN

            ISTRT = XLINDX(JS)
            ISTOP = XLINDX(JS+1)-1

* --------- Insert the row index set for the current supernode.
            DO I = 1, ISTRT, ISTOP
               ROWIND = LINDX(I)
               IF ( RINDEX(ROWIND).GT.0 ) THEN
                  SUB11(IPXK11) = RINDEX(ROWIND)
                  IPXK11 = IPXK11 + 1
               ENDIF
            END DO

            IPSM = XLNZ(JS)
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Insert the coefficient set for the current supernode.
            DO J = JSTRT, JSTOP
               DO I = ISTRT, ISTOP
                  ROWIND = LINDX(I)
                  IF ( RINDEX(ROWIND).GT.0 ) THEN
                     K11(IPK11) = SM(IPSM)
                     IPK11 = IPK11 + 1
                  ENDIF
                  IPSM = IPSM + 1
               END DO
               ISTRT = ISTRT + 1
               IPXSUP = IPXSUP + 1
            END DO

* --------- Update the number of supernodes in K11.
            SUP11T = SUP11T + 1

* --------- Update XSUP11, XSUB11 and XLNZ11
            XSUP11(SUP11T+1) = IPXSUP
            XSUB11(SUP11T+1) = IPXK11
            XLNZ11(SUP11T+1) = IPK11

         ENDIF

      END DO

      RETURN
      END


      SUBROUTINE SPRKA3
     $( NEQ   , NEQ1  , NSUPER, XSUPER, XLINDX, SUPSUP, LINDX , XLNZ  ,
     $  SM    , K11   , RINDEX, LIPERM, PERE2I  )

      IMPLICIT NONE

      LOGICAL
     $  LIPERM

      INTEGER
     $  NEQ, NEQ1, NSUPER,
     $  XSUPER(NSUPER+1), XLINDX(NSUPER+1), SUPSUP(NSUPER), LINDX(*),
     $  XLNZ(NSUPER+1), RINDEX(NEQ), PERE2I(NEQ)

      DOUBLE PRECISION
     $  SM(*), K11(NEQ1,NEQ1)

* --- ------------------------------------------------------------------

      INTEGER
     $  I,IPSM,ISTRT,ISTOP,J,JS,JSTRT,JSTOP,NCOL,NROW

      DOUBLE PRECISION
     $  ZERO
      PARAMETER
     $( ZERO = 0.0D+00 )

* --- ------------------------------------------------------------------

      DO J = 1, NEQ1
         DO I = 1, NEQ1
            K11(I,J) = ZERO
         END DO
      END DO

* --- Compute RINDEX
      DO I = 1, NEQ
         RINDEX(I) = 0
      END DO
      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).LE.0 ) THEN
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
         IF ( LIPERM ) THEN
            I = PERE2I(J)
         ELSE
            I = J
         END IF
         IF ( RINDEX(I).GT.0 ) THEN
            NROW = NROW + 1
            RINDEX(I) = NROW
         ENDIF
      END DO

      DO JS = 1, NSUPER
         IF ( SUPSUP(JS).LE.0 ) THEN

            ISTRT = XLINDX(JS)+1
            ISTOP = XLINDX(JS+1)-1
            JSTRT = XSUPER(JS)
            JSTOP = XSUPER(JS+1)-1

* --------- Insert the coefficient set for current supernode
            IPSM = XLNZ(JS)
            DO J = JSTRT, JSTOP
               NCOL = RINDEX(J)
               IF ( NCOL.GT.0 .AND. NCOL.LE.NEQ1 ) THEN
                  K11(NCOL,NCOL) = SM(IPSM)
               END IF
               IPSM = IPSM + 1
               DO I = ISTRT, ISTOP
                  NROW = RINDEX(LINDX(I))
                  IF ( NCOL.GT.0 .AND. NCOL.LE.NEQ1 .AND.
     $                 NROW.GT.0 .AND. NROW.LE.NEQ1 ) THEN
                     K11(NROW,NCOL) = SM(IPSM)
                     K11(NCOL,NROW) = SM(IPSM)
                  ENDIF
                  IPSM = IPSM + 1
               END DO
               ISTRT = ISTRT + 1
            END DO

         ENDIF
      END DO

      RETURN
      END
