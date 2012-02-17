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
  real(kind=dp), private, dimension(:,:), allocatable :: ank

  integer(kind=i4b), private :: mmax, jmax, jmaxh
  real(kind=dp), private :: pstart

  public :: alff_init, alff_clean, alff_calc, &
    alff_calcp0, alff_calcp1, alff_calcpn, alff_test, alff_test_checksum
  private :: fouriercoeff

contains

  subroutine alff_init(ntrunc)
    use alf_module, only: alf_init

    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

    call alf_init(ntrunc)
    mmax = ntrunc
    allocate(ank(mmax, 0:mmax/2))
    ank(:,:) = 0.0_dp
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

    deallocate(ank,enm,fnm,gnm)
    if (allocated(alff_pnm)) then
      deallocate(alff_pnm)
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
    real(kind=dp), dimension(:), allocatable :: pmm, sinlat, coslat
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
    end do! j

    do j=1, jmaxh

      theta = pih-lat(j) ! lat => colat

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

    deallocate(pmm,pn,sinlat,coslat)

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
    real(kind=dp) :: sqrtnnr

      nmax = size(p1)-1
      do n=3, nmax
        n2 = n/2
        nmod = n - n2*2
        p1(n) = 0.0_dp
        sqrtnnr = 1.0_dp/sqrt(n*(n+1.0_dp))
        do l=0, n2
          k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
          p1(n) = p1(n) + ank(n,l)*k*sqrtnnr*sin(k*theta)
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

  subroutine alff_test(ntrunc,nlat,un)
    use math_module, only: rad2deg=>math_rad2deg
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp), dimension(:), allocatable :: lat, wgt

    real(kind=dp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alff_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    call alff_init(ntrunc)
    call cpu_time(t1)
    call alff_calc(lat)
    call cpu_time(t2)
    print *, "alff_calc cpu time=", t2-t1
    if (present(un)) then
      do j=1, nlat
        write(unit=un,rec=j) alff_pnm(j,:,:)
      end do
    end if
    call alff_clean()

  end subroutine alff_test

  subroutine alff_test_checksum(ntrunc,nlat,un)
    use kind_module, only: dp, i4b
    use math_module, only: pih=>math_pih
    use integer_module, only: swap=>integer_swap
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: alf_checksum, alf_calcps, &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp) :: xx, dd, x, dx, p00, theta
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn, k1, k2
    real(kind=dp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pmm, pn
    real(kind=dp), dimension(:,:), allocatable ::  pjm
    real(kind=dp), dimension(:,:,:), allocatable ::  pjn

    print *, "# ----- alff_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    call alff_init(ntrunc)
    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:ntrunc),pn(0:ntrunc), &
      pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc,2))
    pmm(:) = 0.0_dp
    pn(:) = 0.0_dp

    xx = 1.0_dp
    dd = 0.0_dp
    dx = 0.0_dp
    nn = 0
    mm = 0

    p00 = sqrt(0.5_dp)
    call fouriercoeff(p00)
    pjn(:,:,:) = 0.0_dp
    pjn(:,0,:) = p00

    do j=1, jmaxh
      pmm(0) = p00
      call alf_calcps(coslat(j),dm,pmm)
      pjm(j,:) = pmm(:)
    end do
    deallocate(pmm)

! m = 0
    do j=1, jmaxh
      theta = pih-lat(j) ! lat => colat
      pn(0) = pjm(j,0) ! m = 0
      pn(1) = cm(0)*sinlat(j)*pjm(j,0)
      call alff_calcp0(theta,pn)
      pjn(j,0:mmax,1) = pn(:)
    end do
    do n=0, ntrunc
      x = alf_checksum(wgt,pjn(:,n,1))
      dx = 1.0_dp - x
      if (present(un)) then
        write(unit=un,fmt=*) n, m, x, abs(dx)!, dd
      end if
      if (abs(dx)>dd) then
        xx = x
        dd = abs(dx)
        mm = m
        nn = n
      end if
    end do
! m even 
    k1 = 1
    k2 = 2
    do m=2, mmax-2, 2 ! m even
      do j=1, jmaxh
        pjn(j,m,k2) = pjm(j,m)
        pjn(j,m+1,k2) = cm(m)*sinlat(j)*pjm(j,m)
        call alff_calcpn(m,pjn(j,:,k1),pjn(j,:,k2))
      end do
      do n=m, ntrunc
        x = alf_checksum(wgt,pjn(:,n,k2))
        dx = 1.0_dp - x
        if (present(un)) then
          write(unit=un,fmt=*) n, m, x, abs(dx)!, dd
        end if
        if (abs(dx)>dd) then
          xx = x
          dd = abs(dx)
          mm = m
          nn = n
        end if
      end do
      call swap(k1,k2)
    end do

! m = 1
    do j=1, jmaxh
      theta = pih-lat(j) ! lat => colat
      pn(1) = pjm(j,1) ! m = 1
      pn(2) = cm(1)*sinlat(j)*pjm(j,1)
      call alff_calcp1(theta,pn(:))
      pjn(j,1:mmax,1) = pn(1:mmax)
    end do
    do n=1, ntrunc
      x = alf_checksum(wgt,pjn(:,n,1))
      dx = 1.0_dp - x
      if (present(un)) then
        write(unit=un,fmt=*) n, m, x, abs(dx)!, dd
      end if
      if (abs(dx)>dd) then
        xx = x
        dd = abs(dx)
        mm = m
        nn = n
      end if
    end do
! m odd
    k1 = 1
    k2 = 2
    do m=3, mmax-2, 2
      do j=1, jmaxh
        pjn(j,m,k2) = pjm(j,m)
        pjn(j,m+1,k2) = cm(m)*sinlat(j)*pjm(j,m)  
        call alff_calcpn(m,pjn(j,:,k1),pjn(j,:,k2))
      end do
      do n=m, ntrunc
        x = alf_checksum(wgt,pjn(:,n,k2))
        dx = 1.0_dp - x
        if (present(un)) then
          write(unit=un,fmt=*) n, m, x, abs(dx)!, dd
        end if
        if (abs(dx)>dd) then
          xx = x
          dd = abs(dx)
          mm = m
          nn = n
        end if
      end do
      call swap(k1,k2)
    end do
    print *, "x=", xx, " with max error= ", dd, " at (n,m)=(", nn, ",", mm, ")"

    deallocate(lat,sinlat,coslat,wgt,pn,pjm,pjn)
    call alff_clean()
    
  end subroutine alff_test_checksum

end module alff_module
