module glatwgt_module

! calculates Gaussian points and weights.

! Source: Based on Swartztrauber (2002)
! Author: T. Enomoto
! Usage:
!   Calculate Gaussian points and weights
!     subroutine f_scgaus(x, w, nlat)
!       x(nlat)  Gaussian latitudes
!       w(nlat)  Gaussian weights
! History:
!   TE 27 Apr 2003  Fixed a bug in setting SH colatitudes
!   TE 28 Apr 2003  Ported for AFES
!   TE 24 Apr 2003  Implemented Fourier-Legendre formulation.

  use kind_module, only : dp, i4b
  private

  real(kind=dp), private :: xacc
  integer(kind=i4b), private :: nlat, nlath, nlat2
  real(kind=dp), private, dimension(:), allocatable :: an, dan
  logical, public :: glatwgt_verbose = .false.

  private :: legendre_init, legendre_calc, legendre_clean, newton
  public :: glatwgt_calc, glatwgt_approx, glatwgt_test

contains

! Public procedures

  subroutine glatwgt_calc(x, w, lcolat)
    use math_module, only: pi=>math_pi, pih=>math_pih, rad2deg=>math_rad2deg
    implicit none
  
! returns Gaussian latitudes and Gaussian weights.
! NB. Gaussian colatitudes are used during calculation

    real(kind=dp), dimension(:), intent(out) :: x, w
    logical, optional :: lcolat

    integer(kind=i4b) :: j, jj=1
    real(kind=dp) :: x0, pn, dpn, pn_max=0.0_dp, s

    nlat = min(size(x),size(w))
    nlath = nlat / 2
    nlat = nlath + nlath
    nlat2 = nlat + nlat

    call legendre_init()

!    x0 = pih - pih/(nlat + 1)
    x0 = pih - pih/nlat
    call newton(x0, x(nlath), w(nlath))
    x0 = 3.0_dp*x(nlath) - pi
    call newton(x0, x(nlath-1), w(nlath-1))
    do j = nlath-2, 1, -1
      x0 = x(j+1) + x(j+1) - x(j+2) 
      call newton(x0, x(j), w(j))
    end do
    w(nlat:nlath+1:-1) = w(1:nlath)
    s = sum(w(:))
    w(:) = 2.0_dp*w(:) / s

    if (present(lcolat) .and. lcolat) then
      x(nlath+1:nlat) = pi - x(nlath:1:-1)
    else
      x(1:nlath) = pih-x(1:nlath)
      x(nlath+1:nlat) = -x(nlath:1:-1)
    end if

    if (glatwgt_verbose) then
      do j=1, nlath
        call legendre_calc(pih-x(j), pn, dpn)
!        print *, j, x(j)*rad2deg, pn, w(j)
        pn = abs(pn)
        if (pn>pn_max) then
          pn_max = pn
          jj = j
        end if
      end do
      print *, "Largest error:", pn_max, " at", jj
      s = sum(w(1:nlath))
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
    call glatwgt_calc(lat, wgt)
    if (present(un)) then
      write(unit=un,fmt=*) lat
      write(unit=un,fmt=*) wgt
    end if
    glatwgt_verbose = .false.
    deallocate(lat,wgt)

  end subroutine glatwgt_test

! private procedures

  subroutine legendre_init()
    use machine_module, only : machine_eps
    implicit none

    integer(kind=i4b) :: j
    real(kind=dp) :: t1, t2, b1, b2 

! set tolerance
    xacc = sqrt(machine_eps())
    xacc = xacc * sqrt(xacc) 

    allocate(an(0:nlath), dan(1:nlath))

! Fourier coefficient without scaling factor
! coefficients within [ ] of (8) of Swarztrauber 2002
    an(nlath) = 1.0_dp
    t1 = -1.0_dp
    t2 = nlat + 1.0_dp
    b1 =  0.0_dp
    b2 =  (nlat + nlat) + 1.0_dp
    do j = nlath, 1, -1
      t1 = t1 + 2.0_dp
      t2 = t2 - 1.0_dp
      b1 = b1 + 1.0_dp
      b2 = b2 - 2.0_dp
      an(j - 1) = (t1 * t2) / (b1 * b2) * an(j)
    end do
    do j = 1, nlath
      dan(j) = (j + j) * an(j)
    end do

  end subroutine legendre_init

  subroutine legendre_calc(theta, pn, dpn)
    implicit none

    real(kind=dp), intent(in) :: theta
    real(kind=dp), intent(out) :: pn, dpn

    integer(kind=i4b) :: k, l
    real(kind=dp) :: c0, s0, c, s, c1

    c0 = cos(theta + theta)
    s0 = sin(theta + theta)
    c = c0
    s = s0
    pn = 0.5_dp * an(0)
    dpn = 0.0_dp
    do l=1, nlath
      pn = pn + an(l) * c
      dpn = dpn - dan(l) * s
      c1 = c0 * c - s0 * s
      s  = s0 * c + c0 * s
      c  = c1
    end do
    
  end subroutine legendre_calc

  subroutine legendre_clean()

    deallocate(an, dan)

  end subroutine legendre_clean

  subroutine newton(x0, x, w)
    implicit none

  ! finds the root u

    real(kind=dp), intent(in) :: x0
    real(kind=dp), intent(out) :: x, w

    integer(kind=i4b), parameter :: newton_max = 500

    real(kind=dp) :: y, dy, xx
    integer(kind=i4b) :: i

    x = x0
    do i = 1, newton_max
      xx = x
      call legendre_calc(x, y, dy)
      x = x - y/dy
      if (abs(x - xx) <= xacc * abs(xx)) then
        w = (nlat2 + 1.0_dp)/(dy + y * cos(xx)/sin(xx))**2
        return
      end if
    end do
    print *, "### Error in newton : Too many refinement."
    
  end subroutine newton

end module glatwgt_module
