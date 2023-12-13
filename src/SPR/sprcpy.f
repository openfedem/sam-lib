C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE  SPRCPY  (N,SX,INCX,SY,INCY)
C @(#)sprcpy.f 1.1 95/12/29
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C
      INTEGER SX(*),SY(*)
      INTEGER I,INCX,INCY,IX,IY,N
C
      IF(N.LE.0)RETURN
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END




