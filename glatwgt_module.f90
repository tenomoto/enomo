module glatwgt_module

! calculates Gaussian points and weights.

! Source: Based on Swartztrauber (2003)
! Author: T. Enomoto
! Usage:
!   Calculate Gaussian points and weights
!     subroutine f_scgaus(x, w, jMax)
!       x(jMax)  Gaussian latitudes
!       w(jMax)  Gaussian weights
! History:
!   TE 27 Apr 2003  Fixed a bug in setting SH colatitudes
!   TE 28 Apr 2003  Ported for AFES
!   TE 24 Apr 2003  Implemented Fourier-Legendre formulation.

  use kind_module, only : dp, i4b
  private

!  epsa: absolute tolerance
!  epsr: relative tolerance
  real(kind=dp), private, parameter :: epsa = 1.0e-70_dp, epsr = 1.0e-15_dp
  integer(kind=i4b), private :: jMid
  real(kind=dp), private, dimension(:), allocatable :: an
  logical, public :: glatwgt_verbose = .false.

  private :: legendre_init, legendre_P, legendre_dp, legendre_clean, newton
  public :: glatwgt_calc, glatwgt_approx, glatwgt_test

contains

! Public procedures

  subroutine glatwgt_calc(x,w,ny)
    use math_module, only: pi=>math_pi, pih=>math_pih, rad2deg=>math_rad2deg
    implicit none
  
! returns Gaussian latitudes and Gaussian weights.
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
    do j = 1, jMid
      call legendre_dP(x(j), dpn)
      w(j) = (2.0_dp*jMax + 1.0_dp)/(dpn)**2
    end do
    w(jMax:jMid+1:-1) = w(1:jMid)

    x(1:jMid) = pih-x(1:jMid)
    x(jMid+1:jMax) = -x(jMid:1:-1)

    if (glatwgt_verbose) then
      do j=1, jMid
        call legendre_P(pih-x(j), pn)
!        print *, j, x(j)*rad2deg, pn, w(j)
        pn = abs(pn)
        if (pn>pn_max) then
          pn_max = pn
          jj = j
        end if
      end do
      print *, "Largest error:", pn_max, " at", jj
      s = sum(w(1:jMid))
      print *, "sum of weights:", s, " error=", abs(1.0_dp - s)
    end if

    call legendre_clean()
  
  end subroutine glatwgt_calc

  subroutine glatwgt_approx(lat1, dlat, ny)
  ! calculate values of first lat near SP and interval of
  ! approximate Gaussian latitudes
  ! Ritchie 1987 (24) and (25)
    use math_module, only: pi2=>math_pi2
    implicit none

    real(kind=dp), intent(out) :: lat1, dlat
    integer(kind=i4b), intent(in) :: ny

    dlat = pi2/(2.0_dp*ny+1)
    lat1 = 0.5_dp*dlat*(1.0_dp-ny)

  end subroutine glatwgt_approx

  subroutine glatwgt_test(nlat,un)
    use math_module, only: rad2deg=>math_rad2deg

    integer(kind=i4b) :: nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp), dimension(:), allocatable :: lat, wgt

    allocate(lat(nlat),wgt(nlat))
    print *, "# ----- glatwgt_test() -----"
    glatwgt_verbose = .true.
    call glatwgt_calc(lat,wgt,nlat)
    if (present(un)) then
      write(unit=un,fmt=*) lat
      write(unit=un,fmt=*) wgt
    end if
    glatwgt_verbose = .false.
    deallocate(lat,wgt)

  end subroutine glatwgt_test

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

    real(kind=dp) :: c = 0.0_dp, y, t
    integer(kind=i4b) :: k, l

    pn = 0.0_dp
    do l=0, jMid
      k=l*2 ! k = l*2 + 1 if n odd
      y = an(l) * cos(k*theta) - c
      t = pn + y
      c = (t - pn) - y
      pn = t
    end do
    
  end subroutine legendre_P

  subroutine legendre_dP(theta, dpn)
    implicit none

    real(kind=dp), intent(in) :: theta
    real(kind=dp), intent(out) :: dpn

    real(kind=dp) :: c = 0.0_dp, y, t
    integer(kind=i4b) :: k, l

    dpn = 0.0_dp
    do l=1, jMid
      k=l*2 ! k = l*2 + 1 if n odd
      y =  -k * an(l) * sin(k*theta) - c
      t = dpn + y
      c = (t - dpn) - y
      dpn = t
    end do
    
  end subroutine legendre_dP

  subroutine legendre_clean()

    deallocate(an)

  end subroutine legendre_clean

  subroutine newton(f, df, x0, x, absolute_tolerance, relative_tolerance)
    implicit none

  ! finds the root u

    interface
      subroutine f(x, f_result)
        use kind_module, only: dp, i4b
        real(kind=dp), intent(in) :: x
        real(kind=dp), intent(out) :: f_result
      end subroutine f
      subroutine df(x, f_result)
        use kind_module, only: dp, i4b
        real(kind=dp), intent(in) :: x
        real(kind=dp), intent(out) :: f_result
      end subroutine df
    end interface
    real(kind=dp), intent(in) :: x0
    real(kind=dp), intent(out) :: x
    real(kind=dp), optional, intent(in) :: absolute_tolerance, relative_tolerance

    integer(kind=i4b), parameter :: newton_max = 500
    real(kind=dp) :: ea = epsa, er = epsr

    real(kind=dp) :: y, dy, xx
    integer(kind=i4b) :: i

    if (present(absolute_tolerance)) then
      if (absolute_tolerance < epsa) then
        print *, "### Error in newton: absolute tolerance too small"
        stop
      else
        ea = absolute_tolerance
      end if
    end if
    if (present(relative_tolerance)) then
      if (relative_tolerance < epsr) then
        print *, "### Error in newton: relative tolerance too small"
        stop
      else
        er = relative_tolerance
      end if
    end if

    x = x0
    do i = 1, newton_max
      call f(x, y)
      call df(x, dy)
      y = y/dy
      xx = x
      x = x - y
      if (abs(x-xx)/(ea+er*(abs(x)+abs(xx)))<1) then
        return
      end if
    end do
    print *, "### Error in newton : Too many refinement."
    
  end subroutine newton

end module glatwgt_module
