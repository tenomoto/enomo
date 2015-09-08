module glatwgtq_module

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

  use kind_module, only : qp, i4b
  private

  real(kind=qp), private :: xacc
  integer(kind=i4b), private :: nlat, nlath, nlat2
  real(kind=qp), private, dimension(:), allocatable :: an, dan
  logical, public :: glatwgtq_verbose = .false.

  private :: legendre_init, legendre_calc, legendre_clean, newton
  public :: glatwgtq_calc, glatwgtq_approx, glatwgtq_test

contains

! Public procedures

  subroutine glatwgtq_calc(x, w, lcolat)
    use math_module, only: pi=>math_piq, pih=>math_pihq, rad2deg=>math_rad2degq
    implicit none
  
! returns Gaussian latitudes and Gaussian weights.
! NB. Gaussian colatitudes are used during calculation

    real(kind=qp), dimension(:), intent(out) :: x, w
    logical, optional :: lcolat

    integer(kind=i4b) :: j, jj=1
    real(kind=qp) :: x0, pn, qpn, pn_max=0.0_qp, s

    nlat = min(size(x),size(w))
    nlath = nlat / 2
    nlat = nlath + nlath
    nlat2 = nlat + nlat

    call legendre_init()

    x0 = pih - pih/(nlat + 1)
    call newton(x0, x(nlath), w(nlath))
    x0 = 3.0_qp*x(nlath) - pi
    call newton(x0, x(nlath-1), w(nlath-1))
    do j = nlath-2, 1, -1
      x0 = x(j+1) + x(j+1) - x(j+2) 
      call newton(x0, x(j), w(j))
    end do
    w(nlat:nlath+1:-1) = w(1:nlath)
    s = sum(w(:))
    w(:) = 2.0_qp*w(:) / s

    if (present(lcolat) .and. lcolat) then
      x(nlath+1:nlat) = pi - x(nlath:1:-1)
    else
      x(1:nlath) = pih-x(1:nlath)
      x(nlath+1:nlat) = -x(nlath:1:-1)
    end if

    if (glatwgtq_verbose) then
      do j=1, nlath
        call legendre_calc(pih-x(j), pn, qpn)
!        print *, j, x(j)*rad2deg, pn, w(j)
        pn = abs(pn)
        if (pn>pn_max) then
          pn_max = pn
          jj = j
        end if
      end do
      print *, "Largest error:", pn_max, " at", jj
      s = sum(w(1:nlath))
      print *, "sum of weights:", s, " error=", abs(1.0_qp - s)
    end if

    call legendre_clean()
  
  end subroutine glatwgtq_calc

  subroutine glatwgtq_approx(lat1, dlat, ny)
  ! calculate values of first lat near SP and interval of
  ! approximate Gaussian latitudes
  ! Ritchie 1987 (24) and (25)
    use math_module, only: pi2=>math_pi2q
    implicit none

    real(kind=qp), intent(out) :: lat1, dlat
    integer(kind=i4b), intent(in) :: ny

    dlat = pi2/(2.0_qp*ny+1)
    lat1 = 0.5_qp*dlat*(1.0_qp-ny)

  end subroutine glatwgtq_approx

  subroutine glatwgtq_test(nlat,un)
    use math_module, only: rad2deg=>math_rad2degq

    integer(kind=i4b) :: nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=qp), dimension(:), allocatable :: lat, wgt

    allocate(lat(nlat),wgt(nlat))
    print *, "# ----- glatwgtq_test() -----"
    glatwgtq_verbose = .true.
    call glatwgtq_calc(lat, wgt)
    if (present(un)) then
      write(unit=un,fmt=*) lat
      write(unit=un,fmt=*) wgt
    end if
    glatwgtq_verbose = .false.
    deallocate(lat,wgt)

  end subroutine glatwgtq_test

! private procedures

  subroutine legendre_init()
    use machine_module, only : machine_epsq
    implicit none

    integer(kind=i4b) :: j
    real(kind=qp) :: t1, t2, b1, b2 

! set tolerance
    xacc = sqrt(machine_epsq())
    xacc = xacc * sqrt(xacc) 

    allocate(an(0:nlath), dan(1:nlath))

! Fourier coefficient without scaling factor
! coefficients within [ ] of (8) of Swarztrauber 2002
    an(nlath) = 1.0_qp
    t1 = -1.0_qp
    t2 = nlat + 1.0_qp
    b1 =  0.0_qp
    b2 =  (nlat + nlat) + 1.0_qp
    do j = nlath, 1, -1
      t1 = t1 + 2.0_qp
      t2 = t2 - 1.0_qp
      b1 = b1 + 1.0_qp
      b2 = b2 - 2.0_qp
      an(j - 1) = (t1 * t2) / (b1 * b2) * an(j)
    end do
    do j = 1, nlath
      dan(j) = (j + j) * an(j)
    end do

  end subroutine legendre_init

  subroutine legendre_calc(theta, pn, qpn)
    implicit none

    real(kind=qp), intent(in) :: theta
    real(kind=qp), intent(out) :: pn, qpn

    integer(kind=i4b) :: k, l
    real(kind=qp) :: c0, s0, c, s, c1

    c0 = cos(theta + theta)
    s0 = sin(theta + theta)
    c = c0
    s = s0
    pn = 0.5_qp * an(0)
    qpn = 0.0_qp
    do l=1, nlath
      pn = pn + an(l) * c
      qpn = qpn - dan(l) * s
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

    real(kind=qp), intent(in) :: x0
    real(kind=qp), intent(out) :: x, w

    integer(kind=i4b), parameter :: newton_max = 500

    real(kind=qp) :: y, dy, xx
    integer(kind=i4b) :: i

    x = x0
    do i = 1, newton_max
      xx = x
      call legendre_calc(x, y, dy)
      x = x - y/dy
      if (abs(x - xx) <= xacc * abs(xx)) then
        w = (nlat2 + 1.0_qp)/(dy + y * cos(xx)/sin(xx))**2
        return
      end if
    end do
    print *, "### Error in newton : Too many refinement."
    
  end subroutine newton

end module glatwgtq_module
