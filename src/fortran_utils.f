      subroutine zero_double_array(x,n)
      integer n
      double precision x(n)

      integer i

      do i = 1, n
         x(i) = 0d0
      end do
      end subroutine
