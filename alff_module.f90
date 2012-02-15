module alff_module
! calculates normalised associated Legendre functions
! using four-point reccurence
  use kind_module, only: i4b, dp
  implicit none
  private

! Source: Based on Belousov (1962) and Swartztrauber (2003)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  normalised to 1 by default. factor (-1)**m is not included.

  real(kind=dp), public, dimension(:,:,:), allocatable :: alff_pnm

  real(kind=dp), private, dimension(:,:), allocatable :: enm, fnm, gnm
  real(kind=dp), private, dimension(:), allocatable :: sinlat, coslat, sqrtnnr
  real(kind=dp), private, dimension(:,:), allocatable :: ank

  integer(kind=i4b), private :: mmax, jmax, jmaxh
  real(kind=dp), private :: pstart

  public :: alff_init, alff_clean, alff_calc, &
    alff_calcp0, alff_calcp1, alff_calcpn, alff_test
  private :: fouriercoeff

contains

  subroutine alff_init(ntrunc)
    use alf_module, only: alf_init

    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

    call alf_init(ntrunc)
    mmax = ntrunc
    allocate(ank(mmax, 0:mmax/2),sqrtnnr(mmax))
    ank(:,:) = 0.0_dp
    do n=3, mmax
      sqrtnnr(n) = 1.0_dp/sqrt(n*(n+1.0_dp))
    end do
    allocate(enm(mmax,0:mmax),fnm(mmax,0:mmax),gnm(mmax,0:mmax))
    enm(:,:) = 0.0_dp
    fnm(:,:) = 0.0_dp
    gnm(:,:) = 0.0_dp
    do m=2, mmax
      do n=m, mmax
        gnm(n,m) = 1.0_dp / ((n+m-1.0_dp)*(n+m))
        enm(n,m) = gnm(n,m)*(2.0_dp*n+1.0_dp)/(2.0_dp*n-3.0_dp)
        fnm(n,m) = sqrt(enm(n,m)*(n-m)*(n-m-1.0_dp))
        enm(n,m) = sqrt(enm(n,m)*(n+m-2.0_dp)*(n+m-3.0_dp))
        gnm(n,m) = sqrt(gnm(n,m)*(n-m+1.0_dp)*(n-m+2.0_dp))
      end do
    end do

  end subroutine alff_init

  subroutine alff_clean()
    use alf_module, only: alf_clean

    deallocate(ank,alff_pnm)
    if (allocated(enm)) then
      deallocate(enm)
    end if
    if (allocated(fnm)) then
      deallocate(fnm)
    end if
    if (allocated(gnm)) then
      deallocate(gnm)
    end if
    if (allocated(sinlat)) then
      deallocate(sinlat)
    end if
    if (allocated(coslat)) then
      deallocate(coslat)
    end if
    if (allocated(sqrtnnr)) then
      deallocate(sqrtnnr)
    end if

    call alf_clean()

  end subroutine alff_clean
  
  subroutine alff_calc(lat,p00)
    use math_module, only: pih=>math_pih
    use integer_module, only: swap=>integer_swap
    use alf_module, only: cm=>alf_cm, dm=>alf_dm, alf_calcps

    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, k1, k2
    real(kind=dp) :: theta
    real(kind=dp), dimension(:), allocatable :: pmm
    real(kind=dp), dimension(:,:), allocatable :: pn
    
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    call fouriercoeff(pstart)
    jmax = size(lat)
    jmaxh = jmax/2
    if (jmaxh>=1) then
      allocate(alff_pnm(jmaxh, -1:mmax, 0:mmax))
      alff_pnm(:,:,:) = 0.0_dp
      alff_pnm(:,0,0) = pstart
    end if
    allocate(sinlat(jmaxh),coslat(jmaxh))
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:mmax),pn(0:mmax,2))
    pmm(:) = 0.0_dp
    pn(:,:) = 0.0_dp

! calculate Pnm
    do j=1, jmaxh

      theta = pih-lat(j) ! lat => colat

! Pmm and Pm,n=m+1
      pmm(0) = pstart
      call alf_calcps(coslat(j),dm,pmm)

      pn(0,1) = pmm(0) ! m = 0
      pn(1,1) = cm(0)*sinlat(j)*pmm(0)
      call alff_calcp0(theta,pn(:,1))
      alff_pnm(j,0:mmax,0) = pn(0:mmax,1)
      k1 = 1
      k2 = 2
      do m=2, mmax-2, 2 ! m even
        pn(m,k2) = pmm(m)
        pn(m+1,k2) = cm(m)*sinlat(j)*pmm(m)  
        call alff_calcpn(m,pn(:,k1),pn(:,k2))
        alff_pnm(j,m:mmax,m) = pn(m:mmax,k2)
        call swap(k1,k2)
      end do

      pn(1,1) = pmm(1) ! m = 1
      pn(2,1) = cm(1)*sinlat(j)*pmm(1)
      call alff_calcp1(theta,pn(:,1))
      alff_pnm(j,1:mmax,1) = pn(1:mmax,1)
      k1 = 1
      k2 = 2
      do m=3, mmax-2, 2 ! m odd
        pn(m,k2) = pmm(m)
        pn(m+1,k2) = cm(m)*sinlat(j)*pmm(m)  
        call alff_calcpn(m,pn(:,k1),pn(:,k2))
        alff_pnm(j,m:mmax,m) = pn(m:mmax,k2)
        call swap(k1,k2)
      end do
      alff_pnm(j,mmax,mmax-1) = cm(mmax-1)*sinlat(j)*pmm(mmax-1)
      do m=mmax-1, mmax
        alff_pnm(j,m,m) = pmm(m)
      end do

    end do! j

    deallocate(pmm,pn)

  end subroutine alff_calc

  subroutine fouriercoeff(p0)

    real(kind=dp), intent(in), optional :: p0

    integer(kind=i4b) :: n, l, k, n2

! calculate fourier coefficients for Pn
    ank(2,1) = 0.75_dp*sqrt(5.0_dp) * p0
    do n=3, mmax
      ank(n,n/2) = &
        sqrt(1.0_dp-1.0_dp/(4.0_dp*n*n)) * ank(n-1,(n-1)/2)
    end do
    do n=2, mmax
      n2  = n/2
      do k=1, n2
        l = 2*k
        ank(n,n2-k) = (l-1.0_dp)*(2.0_dp*n-l+2.0_dp)/&
          (l*(2.0_dp*n-l+1.0_dp)) * ank(n,n2-k+1)
      end do
      if (n==n2*2) then
        ank(n,0) = 0.5_dp*ank(n,0)
      end if
    end do

  end subroutine fouriercoeff

  subroutine alff_calcp0(theta,p0)
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(0:), intent(inout) :: p0

    integer(kind=i4b) :: n, l, k, n2, nmod, nmax

      nmax = size(p0)-1
      do n=2, nmax
        n2 = n/2
        nmod = n - n2*2
        p0(n) = 0.0_dp
        do l=0, n2
          k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
          p0(n) = p0(n) + ank(n,l)*cos(k*theta)
        end do
      end do

  end subroutine alff_calcp0

  subroutine alff_calcp1(theta,p1)
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(0:), intent(inout) :: p1

    integer(kind=i4b) :: n, l, k, n2, nmod, nmax

      nmax = size(p1)-1
      do n=3, nmax
        n2 = n/2
        nmod = n - n2*2
        p1(n) = 0.0_dp
        do l=0, n2
          k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
          p1(n) = p1(n) + ank(n,l)*k*sqrtnnr(n)*sin(k*theta)
        end do
      end do

  end subroutine alff_calcp1

  subroutine alff_calcpn(m,pn0,pn1)
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(0:), intent(in) :: pn0
    real(kind=dp), dimension(0:), intent(inout) :: pn1

    integer(kind=i4b) :: n, nmax

    if (m<2) then
      return
    end if
    nmax = size(pn1) - 1
!    do n=m, nmax
    do n=m+2, nmax
      pn1(n) = enm(n,m)*pn0(n-2) + fnm(n,m)*pn1(n-2) - gnm(n,m)*pn0(n)
    end do

  end subroutine alff_calcpn

  subroutine alff_test(nt_low,nt_high)
    use math_module, only: rad2deg=>math_rad2deg
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm, alf_checksum

    integer(kind=i4b), intent(in), optional :: nt_low, nt_high

    integer(kind=i4b), parameter :: d = 10
    real(kind=dp), dimension(:), allocatable :: lat, wgt

    integer(kind=i4b) :: ntrunc_low = 159, ntrunc_high = 2159
    integer(kind=i4b) :: nlat, n, m, j, jmaxh, nn, mm, ntrunc
    real(kind=dp) :: x, dx, xx, dd, theta, t1, t2

    print *, "# ----- alff_test() -----" 
    if (present(nt_low)) then
      ntrunc = nt_low
    else
      ntrunc = ntrunc_low
    end if
    nlat = (ntrunc+1)*3/2
    allocate(lat(nlat),wgt(nlat))
    jmaxh = size(alff_pnm,1)
    call glatwgt_calc(lat,wgt,nlat)
    call alff_init(ntrunc)
    call cpu_time(t1)
    call alff_calc(lat)
    call cpu_time(t2)
    print *, "alff_calc cpu time=", t2-t1
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(1)=", lat(1)*rad2deg 
    print *, "m, pm(n=m,j=1), pm(n=ntrunc,j=1)"
    do m=0, ntrunc, ntrunc/d
      print *, m, alff_pnm(1,m,m), alff_pnm(1,ntrunc,m)
    end do
    do m=ntrunc-1, ntrunc
      print *, m, alff_pnm(1,m,m), alff_pnm(1,ntrunc,m)
    end do
    print *, "n, pn(m=0,j=1), pn(m=1,j=1)"
    do n=0, ntrunc, ntrunc/d 
      print *, n, alff_pnm(1,n,0), alff_pnm(1,n,1)
    end do
    do n=ntrunc-1, ntrunc
      print *, n, alff_pnm(1,n,0), alff_pnm(1,n,1)
    end do
    print *, "x=\int pnm pnm dx error"
    xx = 1.0_dp
    dd = 0.0_dp
    nn = 0
    mm = 0
    do m=0, ntrunc, ntrunc
      do n=m, ntrunc, ntrunc/d
        x = alf_checksum(wgt(1:jmaxh),alff_pnm(:,n,m))
        dx = 1.0_dp-x
        if (abs(dx)>dd) then
          xx = x
          dd = dx
          mm = m
          nn = n
        end if
      end do
    end do
    print *, "x=", xx, " with max error= ", dd, " at (n,m)=(", nn, ",", mm, ")"
    deallocate(lat,wgt)
    call alff_clean()

    if (present(nt_high)) then
      ntrunc = nt_high
    else
      ntrunc = ntrunc_high
    end if
    if (ntrunc<0) then
      return
    end if
    nlat = (ntrunc+1)*3/2
    j = nlat/10
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(j=", j, ")=", lat(j)*rad2deg 
    call alff_init(ntrunc)
    m = ntrunc/3
    call alff_calc((/lat(j:j+1)/))
    print *, "n, m, pn(m=", m, ",j=",j,")"
    do n=m, ntrunc, ntrunc/d
      print *, n, m, alff_pnm(1,n,m)
    end do
    do n=ntrunc-1, ntrunc
      print *, n, m, alff_pnm(1,n,m)
    end do
    call alff_clean()
    print *, "# --------------------" 

  end subroutine alff_test

end module alff_module
