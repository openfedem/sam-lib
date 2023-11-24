C     SPDX-License-Identifier: Apache-2.0
C
      SUBROUTINE I7COPY (N, IX, INCX, IY)
      INTEGER IX(*), IY(*)
      if (INCX .EQ. 0) then
         do I = 1, N
            IY(I) = IX(1)
         end do
      else if (INCX .EQ. 1) then
         M = MOD(N,7)
         do I = 1, M
            IY(I) = IX(I)
         end do
         do I = M+1, N, 7
            IY(I)   = IX(I)
            IY(I+1) = IX(I+1)
            IY(I+2) = IX(I+2)
            IY(I+3) = IX(I+3)
            IY(I+4) = IX(I+4)
            IY(I+5) = IX(I+5)
            IY(I+6) = IX(I+6)
         end do
      end if
      RETURN
      END
