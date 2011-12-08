module monotonecubic_module
  use kind_module, only: i4b, dp
  implicit none

! History:
!   2009-5-15 T. Enomoto <eno@jamstec.go.jp>
!     * first version

! References:
!  Williamson, D. L. and P. J. Rasch, 1989: Two-dimensional semi-Lagrangian transport
!    with shape-preserving interpolation
!  Fritsch, F. N. and R. E. Carlson, 1980: Monotone piecewise cubic interpolation.
!    SIAM J. Numer. Anal., 17, 238--246

contains

  subroutine monotonecubic_init(f, dxr, d)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: f, dxr ! dxr = 1/(x_i+1-x_i)
    real(kind=dp), dimension(:), intent(inout) :: d

    real(kind=dp) :: delta, alpha, beta, tau
    integer(kind=i4b) :: i, n

    n = size(f)
    do i=1,  n-1
      delta = (f(i+1)-f(i))*dxr(i)
      if (delta==0) then
        d(i:i+1) = 0.0_dp
      else
        alpha = d(i)/delta 
        beta = d(i+1)/delta 
        if ((alpha>=0).and.(beta>=0)) then
          if ((alpha<=3).and.(beta<=3)) then
            cycle
          else
            tau = 3/sqrt(alpha*alpha+beta*beta)
            d(i:i+1) = tau*d(i:i+1)
          end if
        else
          d(i:i+1) = 0.0_dp
        end if
      end if
    end do

  end subroutine monotonecubic_init

end module monotonecubic_module
