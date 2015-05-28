module alf_module
! calculates normalised associated Legendre functions
! using three point recurrence
  use kind_module, only: i4b, dp
  implicit none
  private

! Source: Based on Fukushima (2011)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  normalised to 1 by default. factor (-1)**m is not included.


  real(kind=dp), public, dimension(:,:), allocatable :: alf_anm, alf_bnm
  real(kind=dp), public, dimension(:), allocatable :: alf_cm, alf_dm
  
  real(kind=dp), private, parameter :: dbl_min = 2.0_dp**(-1022)

  integer(kind=i4b), private :: alf_ntrunc = 0
  real(kind=dp), private :: pstart
  integer, private :: retain = 0

  public :: alf_init, alf_clean, alf_calc, alf_calc_m, &
    alf_calcps, alf_calcpn, alf_checksum, alf_test, alf_test_checksum

contains

  subroutine alf_init(ntrunc)
    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

!    print *, "alf_init"
    if (alf_ntrunc==ntrunc) then
      retain = retain + 1
!      print *, retain
      return
    end if
    alf_ntrunc = ntrunc
!    print *, "Allocating alf_anm, alf_bnm, alf_cm, alf_dm"

    allocate(alf_cm(0:ntrunc),alf_dm(ntrunc))
    do m=0, ntrunc-1
      alf_cm(m) = sqrt(2.0_dp*m+3.0_dp)
    end do
    do m=1, ntrunc
      alf_dm(m) = sqrt(1.0_dp + 0.5_dp/real(m,kind=dp))
    end do
    allocate(alf_anm(ntrunc,0:ntrunc),alf_bnm(ntrunc,0:ntrunc))
    do m=0, ntrunc
      do n=m+2, ntrunc
        alf_anm(n,m) = sqrt((2.0_dp*n+1.0_dp)/real(n**2-m**2,kind=dp))
        alf_bnm(n,m) = alf_anm(n,m)*sqrt((n+m-1.0_dp)*(n-m-1.0_dp)/(2.0_dp*n-3.0_dp))
        alf_anm(n,m) = alf_anm(n,m)*sqrt(2.0_dp*n-1.0_dp)
      end do
    end do
    retain = retain + 1
!    print *, retain

  end subroutine alf_init

  subroutine alf_clean()

!    print *, "alf_clean"
    retain = retain - 1 
!    print *, retain
    if (retain<1) then
!      print *, "Deallocating alf_anm, alf_bnm, alf_cm, alf_dm"
      deallocate(alf_anm,alf_bnm,alf_cm,alf_dm)
      alf_ntrunc = 0
    end if

  end subroutine alf_clean

  subroutine alf_calc(lat,alf,p00)
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), dimension(0:,0:,:), intent(out) :: alf
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh, mmax

    real(kind=dp), dimension(0:size(alf,1)-1) :: pmm, pnm
    real(kind=dp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alf,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alf(:,:,:) = 0.0_dp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0) = pstart
    do j=1, min(jmax,size(alf,3))
      call alf_calcps(coslat(j),alf_dm,pmm)
      do m=0, mmax-1
        pnm(m) = pmm(m)
        alf(m,m,j) = pmm(m)
        call alf_calcpn(sinlat(j),m,alf_anm(:,m),alf_bnm(:,m),alf_cm(m),pnm)
        alf(m+1:mmax,m,j) = pnm(m+1:mmax)
      end do
      alf(mmax,mmax,j) = pmm(mmax)
    end do

  end subroutine alf_calc

  subroutine alf_calc_m(m,lat,alfm,p00)
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), dimension(0:,:), intent(out) :: alfm
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, n, jmax, jmaxh, mmax

    real(kind=dp), dimension(0:size(alfm,1)-1) :: pmm, pnm
    real(kind=dp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alfm,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alfm(:,:) = 0.0_dp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0) = pstart
    do j=1, min(jmax,size(alfm,2))
      call alf_calcps(coslat(j),alf_dm,pmm)
      if (m/=mmax) then
        pnm(m) = pmm(m)
        alfm(m,j) = pmm(m)
        call alf_calcpn(sinlat(j),m,alf_anm(:,m),alf_bnm(:,m),alf_cm(m),pnm)
        alfm(m+1:mmax,j) = pnm(m+1:mmax)
      else 
        alfm(mmax,j) = pmm(mmax)
      end if
    end do

  end subroutine alf_calc_m

  subroutine alf_calcps(u,d,ps)

    real(kind=dp), intent(in) :: u ! coslat
    real(kind=dp), dimension(:), intent(in) :: d
    real(kind=dp), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
      if (abs(ps(m))<dbl_min) then
        exit
      end if
    end do
    ps(m:) = 0.0_dp

  end subroutine alf_calcps

  subroutine alf_calcpn(t,m,an,bn,c,pn)

    real(kind=dp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(:), intent(in) :: an, bn
    real(kind=dp), intent(in) :: c
    real(kind=dp), dimension(0:), intent(inout) :: pn

    integer(kind=i4b) :: n, nmax

    nmax = size(pn)-1
    if (m+1>nmax) then
      return
    endif

    pn(m+1) = c*t*pn(m)
    do n=m+2, nmax
      pn(n) = an(n)*t*pn(n-1)-bn(n)*pn(n-2)
    end do

  end subroutine alf_calcpn

  function alf_checksum(wgt,pj) result(x)

    real(kind=dp), dimension(:), intent(in) :: wgt
    real(kind=dp), dimension(:), intent(in) :: pj

    real(kind=dp) :: x
    integer(kind=i4b) :: jmaxh

    jmaxh = size(pj)
    x = 2.0_dp*sum(wgt(1:jmaxh)*pj(:)*pj(:))

  end function alf_checksum
  
  subroutine alf_test(ntrunc,nlat,un)
    use math_module, only: rad2deg=>math_rad2deg
    use glatwgt_module, only: glatwgt_calc

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp), dimension(:), allocatable :: lat, wgt
    real(kind=dp), dimension(:,:,:), allocatable :: alf
    real(kind=dp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alf_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    allocate(alf(0:ntrunc,0:ntrunc,nlat/2))
    call glatwgt_calc(lat,wgt,nlat)
    call alf_init(ntrunc)
    call cpu_time(t1)
    call alf_calc(lat(1:nlat/2),alf)
    call cpu_time(t2)
    print *, "alf_calc cpu time=", t2-t1
    if (present(un)) then
      write(unit=un,rec=1) alf
    end if
    deallocate(alf)
    deallocate(lat,wgt)
    call alf_clean()
    
  end subroutine alf_test

  subroutine alf_test_checksum(ntrunc,nlat,un)
    use kind_module, only: dp, i4b
    use glatwgt_module, only: glatwgt_calc

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp) :: xx, dd, x, dx, p00
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn
    real(kind=dp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pmm, pnm
    real(kind=dp), dimension(:,:), allocatable ::  pjm, pjn

    print *, "# ----- alf_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:ntrunc),pnm(0:ntrunc), &
      pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc))

    xx = 1.0_dp
    dd = 0.0_dp
    dx = 0.0_dp
    nn = 0
    mm = 0

    p00 = sqrt(0.5_dp)
    call alf_init(ntrunc)
    do j=1, jmaxh
      pmm(0) = p00
      call alf_calcps(coslat(j),alf_dm,pmm)
      pjm(j,:) = pmm(:)
    end do
    do m=0, ntrunc
      do j=1, jmaxh
        pnm(m) = pjm(j,m)
        pjn(j,m) = pjm(j,m)
        if (m<ntrunc) then
          call alf_calcpn(sinlat(j),m,alf_anm(:,m),alf_bnm(:,m),alf_cm(m),pnm)
          pjn(j,m+1:ntrunc) = pnm(m+1:ntrunc)
        end if
      end do
      do n=m, ntrunc
        x = alf_checksum(wgt,pjn(:,n))
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
    end do
    print *, "x=", xx, " with max error= ", dd, " at (n,m)=(", nn, ",", mm, ")"

    deallocate(lat,sinlat,coslat,wgt,pmm,pnm,pjm,pjn)
    call alf_clean()
    
  end subroutine alf_test_checksum

end module alf_module
