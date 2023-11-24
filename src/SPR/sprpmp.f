C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE SPRPMP ( MSICA , MTREES,
     $                    MEQN  , MPAR  ,
     $                    MSPAR , MDOFNC  )
C
C **********************************************************************
C
C   S A M  LIBRARY ROUTINE :  SPRPMP                  GROUP 9 / PUBLIC
C
C     T A S K :  To determine some pointers in MSICA, MTREES and MSIFA,
C                and place them in MSPAR, and to establish control
C                arrays NODMAP and SUPMAP and store them in MDOFNC.
C
C
C     ROUTINES CALLED/REFERENCED :   None
C
C     PROGRAMMED BY :   Kolbein Bell
C     DATE/VERSION  :   94-02-19 / 1.0
*
*     ACD - 98-08-24 - MSPAR(50) -> MSPAR(*), and
*                      added superelement option.
C
C     KMO - 03-02-27 - MSPAR(47) = 1 + MSPAR(58)*MSPAR(8), and
C                      MSPAR(39) = MPAR(20) + MSPAR(24) + MSPAR(8) + ...
C
C **********************************************************************
C
      IMPLICIT          NONE
      INTEGER           MDOFNC(*),MEQN(*),MPAR(*),MSICA(*),MTREES(*),
     $                  MSPAR(*)
C                                                ! local variables
      INTEGER           I,IDOF,IEQ,IP,J,JPN,JPX
C ----------------------------------------------------------------------
C
C                                                ! pointers in MSICA
      MSPAR(41) = 1
      MSPAR(43) = MSPAR(41) + MSPAR(5) + MSPAR(6) + 2
      MSPAR(44) = MSPAR(43) + MSPAR(6) + 1
      MSPAR(42) = MSPAR(44) + MSPAR(6) + MSPAR(8)
C                                                ! pointers in MTREES
      MSPAR(45) = 1
      MSPAR(46) = MSPAR(45) + 3*(MSPAR(11) + 1)
C                                                ! pointers in MSIFA
      MSPAR(47) = 1 + MSPAR(58)*MSPAR(8)
      MSPAR(48) = MSPAR(47) + MSPAR(15)
C
C --- Control matrices NODMAP AND SUPMAP
C
      JPX = MSPAR(43)
      JPN = MSPAR(44)
C
      DO 40 I=1,MSPAR(6)
         DO 20 J=MSICA(JPX+I-1),MSICA(JPX+I)-1
            IDOF = MSICA(JPN+J-1)
            IEQ  = MEQN(IDOF)
            MDOFNC(IEQ) = I
   20    CONTINUE
   40 CONTINUE
C
      IP = MSPAR(8) + 1
      DO 80 I=1,MSPAR(11)
         DO 60 J=MTREES(I),MTREES(I+1)-1
            MDOFNC(IP+J-1) = I
   60    CONTINUE
   80 CONTINUE
C
      MSPAR(39) = MPAR(20) + MSPAR(24) + MSPAR(8) +
     $            MSPAR(24)*(MSPAR(24)+1)/2 +
     $            MSPAR(25)*(MSPAR(25)+1)/2
C
      RETURN
      END
