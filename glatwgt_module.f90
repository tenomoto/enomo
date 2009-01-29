module glatwgt_module

! calculates Gaussian points and weights.

! Source: Based on Swartztrauber (2003)
! Author: T. Enomoto
! Usage:
!   Calculate Gaussian points and weights
!     subroutine f_scgaus(x, w, jMax)
!       x(jMax)  cos(Gaussian colatitudes) = sin(Gaussian latitudes)
!       w(jMax)  Gaussian weights
! History:
!   TE 27 Apr 2003  Fixed a bug in setting SH colatitudes
!   TE 28 Apr 2003  Ported for AFES
!   TE 24 Apr 2003  Implemented Fourier-Legendre formulation.

  use type_module, only : dp, i4b
  private

!  xacc_min  double minimum accuracy available
  real(kind=dp), private, parameter :: xacc_min = 1.0e-15_dp
  integer(kind=i4b), private :: jMid
  real(kind=dp), private, dimension(:), allocatable :: an

  private :: legendre_init, legendre_P, legendre_dp, legendre_clean, newton
  public :: glatwgt_calc, glatwgt_approx

contains

! Public procedures

  subroutine glatwgt_calc(x,w,ny)
    use math_module, only: pi=>math_pi, pih=>math_pih
    implicit none
  
! returns sin(Gaussian latitudes) between 1 and -1
! and Gaussian weights.
! NB. Gaussian colatitudes are used during calculation

    real(kind=dp), dimension(:), intent(out) :: x, w
    integer(kind=i4b), intent(in), optional :: ny

    integer(kind=i4b) :: l, j, jj=1, jMax
    real(kind=dp) :: guess, pn, dpn, pn_max=0.0_dp, s

    if (present(ny)) then
      jMax = ny
    else
      jMax = min(size(x),size(w))
    end if
    jMid = jMax/2

    call legendre_init()

    guess = pih - pih/(jMax + 1)
    call newton(legendre_P, legendre_dP, guess, x(jMid))
    guess = 3.0_dp*x(jMid) - pi
    call newton(legendre_P, legendre_dP, guess, x(jMid-1))
    do l = jMid-2, 1, -1
      guess = 2*x(l+1) - x(l+2) 
      call newton(legendre_P, legendre_dP, guess, x(l))
    end do
    x(jMax:jMid+1:-1) = pi-x(1:jMid)
    do j = 1, jMid
      call legendre_dP(x(j), dpn)
      w(j) = (2.0_dp*jMax + 1.0_dp)/(dpn)**2
    end do

! ***** debug
!  print *, "j, lat, pn, w"
!  do j=1, jMid
!    call legendre_P(x(j), pn)
!    print *, j, (pi/2-x(j))*180/pi, pn, w(j)
!    pn = abs(pn)
!    if (pn>pn_max) then
!      pn_max = pn
!      jj = j
!    end if
!  end do
!  print *, "Largest error:", pn_max, " at", jj
  s = sum(w(1:jMid))
  print *, "sum of weights:", s, " error=", abs(1.0_dp - s)
! *****

    w(jMax:jMid+1:-1) = w(1:jMid)
    x(:) = cos(x(jMax:1:-1))

    call legendre_clean()
  
  end subroutine glatwgt_calc

  subroutine glatwgt_approx(lat1, dlat, ny)
! calculate values of first lat near SP and interval of
! approximate Gaussian latitudes
    use math_module, only: pi2=>math_pi2
    implicit none

    real(kind=dp), intent(out) :: lat1, dlat
    integer(kind=i4b), intent(in) :: ny

    dlat = pi2/(2.0_dp*ny+1)
    lat1 = 0.5_dp*dlat*(1.0_dp-ny)

  end subroutine glatwgt_approx

! private procedures

subroutine legendre_init()
  implicit none

  real(kind=dp), parameter :: a0sq = 1.5_dp
  integer(kind=i4b) :: k, l, n

  n = jMid*2
  allocate(an(0:jMid))
  an(jMid) = a0sq
  do k=2, n
    an(jMid) = (1.0_dp-1.0_dp/(4.0_dp*k*k))*an(jMid)
  end do
  an(jMid) = sqrt(an(jMid))

  do k=1, jMid
    l = 2*k
    an(jMid-k) = (l-1.0_dp)*(2.0_dp*n-l+2.0_dp)/(l*(2.0_dp*n-l+1.0_dp)) * an(jMid-k+1) 
  end do
  an(0) = 0.5_dp*an(0)

end subroutine legendre_init

subroutine legendre_P(theta, pn)
  implicit none

  real(kind=dp), intent(in) :: theta
  real(kind=dp), intent(out) :: pn
  integer(kind=i4b) :: k, l

  pn = 0.0
  do l=0, jMid
    k=l*2 ! k = l*2 + 1 if n odd
    pn = pn + an(l) * cos(k*theta)
  end do
  
end subroutine legendre_P

subroutine legendre_dP(theta, dpn)
  implicit none

  real(kind=dp), intent(in) :: theta
  real(kind=dp), intent(out) :: dpn
  integer(kind=i4b) :: k, l

  dpn = 0.0
  do l=1, jMid
    k=l*2 ! k = l*2 + 1 if n odd
    dpn = dpn - k * an(l) * sin(k*theta)
  end do
  
end subroutine legendre_dP

subroutine legendre_clean()

  deallocate(an)

end subroutine legendre_clean

subroutine newton(f, df, x0, x, tolerance)
  implicit none

! finds the root u

  interface
    subroutine f(x, f_result)
      use type_module, only: dp, i4b
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: f_result
    end subroutine f
    subroutine df(x, f_result)
      use type_module, only: dp, i4b
      real(kind=dp), intent(in) :: x
      real(kind=dp), intent(out) :: f_result
    end subroutine df
  end interface
  real(kind=dp), intent(in) :: x0
  real(kind=dp), intent(out) :: x
  real(kind=dp), optional, intent(in) :: tolerance

  integer(kind=i4b), parameter :: newton_max = 500
  real(kind=dp) :: xacc = xacc_min

  real(kind=dp) :: y, dy
  integer(kind=i4b) :: i

  if (present(tolerance)) then
    if (tolerance < xacc_min) then
      print *, "### Error in newton: tolerance too small"
      stop
    else
      xacc = tolerance
    end if
  else
    xacc = xacc_min
  end if

  x = x0
  do i = 1, newton_max
    call f(x, y)
    call df(x, dy)
    y = y/dy
    if (abs(y) < xacc) then
      return
    end if
    x = x - y
  end do
  print *, "### Error in newton : Too many refinement."
  
end subroutine newton

end module glatwgt_module
