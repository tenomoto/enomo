module alfx_module
! calculate associated Legendre functions
! in extended exponent floating numbers
  use kind_module, only: dp, i4b
  use xreal_module, only: xreal_type
  implicit none
  private

! Source: Based on Fukushima (2011)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  normalised to 1 by default. factor (-1)**m is not included.

  real(kind=dp), public, dimension(:,:,:), allocatable :: alfx_pnm
  real(kind=dp), private :: pstart

  real(kind=dp), private, dimension(:), allocatable :: sinlat, coslat
  integer(kind=i4b), private :: mmax

  interface alfx_calcps
    module procedure alfsp, alfsx
  end interface

  interface alfx_calcpn
    module procedure alfmp, alfmx
  end interface

  public :: alfx_init, alfx_clean, alfx_calc, alfx_calc_inline, &
    alfx_calcps, alfx_calcpn, alfx_test, alfx_test_checksum
  private :: alfsp, alfmp, alfsx, alfmx

contains

  subroutine alfx_init(ntrunc)
    use alf_module, only: alf_init
    integer(kind=i4b), intent(in) :: ntrunc

    mmax = ntrunc

    call alf_init(ntrunc)

  end subroutine alfx_init

  subroutine alfx_clean()
    use alf_module, only: alf_clean

    if (allocated(alfx_pnm)) then
      deallocate(alfx_pnm)
    end if
    if (allocated(sinlat)) then
      deallocate(sinlat)
    end if
    if (allocated(coslat)) then
      deallocate(coslat)
    end if

    call alf_clean()

  end subroutine alfx_clean

  subroutine alfx_calc(lat,p00)
    use xreal_module, only: xreal_type, assignment(=)
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh

    type(xreal_type), dimension(:), allocatable :: pmm, pnm

    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    jmax = size(lat)
    jmaxh = jmax/2
    if (jmaxh<1) then
      return
    end if

    allocate(alfx_pnm(jmaxh, -1:mmax, 0:mmax))
    alfx_pnm(:,:,:) = 0.0_dp
    allocate(sinlat(jmaxh),coslat(jmaxh))
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:mmax),pnm(0:mmax))
    pmm(0) = pstart
    do j=1, jmaxh
      call alfx_calcps(coslat(j),dm,pmm)
      do m=0, mmax
        pnm(m) = pmm(m)
        alfx_pnm(j,m,m) = pmm(m)
        call alfx_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm)
        do n=m+1, mmax
          alfx_pnm(j,n,m) = pnm(n)
        end do
      end do
    end do
    deallocate(pmm,pnm)

  end subroutine alfx_calc

  subroutine alfx_calc_inline(lat,p00)
    use xreal_module, only: big=>xreal_big, bigi=>xreal_bigi
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh

    real(kind=dp), dimension(:), allocatable :: pmm, pnm
    integer(kind=i4b), dimension(:), allocatable :: ipmm, ipnm

    if (present(p00)) then
      pstart = p00
    end if
    jmax = size(lat)
    jmaxh = jmax/2
    if (jmaxh<1) then
      return
    end if

    allocate(alfx_pnm(jmaxh, -1:mmax, 0:mmax))
    alfx_pnm(:,:,:) = 0.0_dp
    allocate(sinlat(jmaxh),coslat(jmaxh))
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:mmax),pnm(0:mmax),ipmm(0:mmax),ipnm(0:mmax))
    pmm(0)  = pstart
    ipmm(0) = 0
    do j=1, jmaxh
      call alfx_calcps(coslat(j),dm,pmm,ipmm)
      do m=0, mmax
        pnm(m)  = pmm(m)
        ipnm(m) = ipmm(m)
        select case (ipmm(m))
          case(0)
            alfx_pnm(j,m,m) = pmm(m)
          case(:-1)
            alfx_pnm(j,m,m) = pmm(m)*bigi
          case(1:)
            alfx_pnm(j,m,m) = pmm(m)*big
        end select
        call alfx_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm,ipnm)
        do n=m+1, mmax
          select case (ipnm(n))
            case(0)
              alfx_pnm(j,n,m) = pnm(n)
            case(:-1)
              alfx_pnm(j,n,m) = pnm(n)*bigi
            case(1:)
              alfx_pnm(j,n,m) = pnm(n)*big
          end select
        end do
      end do
    end do
    deallocate(pmm,pnm,ipmm,ipnm)

  end subroutine alfx_calc_inline

  subroutine alfsp(u,d,ps)
    use xreal_module, only: xreal_type, assignment(=), operator(*)

    real(kind=dp), intent(in) :: u ! coslat
    real(kind=dp), dimension(:), intent(in) :: d
    type(xreal_type), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
    end do

  end subroutine alfsp

  subroutine alfsx(u,d,ps,ips)
    use xreal_module, only: big=>xreal_big, bigi=>xreal_bigi, &
      bigs=>xreal_bigs, bigsi=>xreal_bigsi

    real(kind=dp), intent(in) :: u ! coslat
    real(kind=dp), dimension(:), intent(in) :: d
    real(kind=dp), dimension(0:), intent(inout) :: ps
    integer(kind=i4b), dimension(0:), intent(inout) :: ips

    integer(kind=i4b) :: m, nmax, ix
    real(kind=dp) :: x, y

    nmax = size(ps)-1
    x  = ps(0)
    ix = 0
    do m=1, nmax
      x = (d(m)*u)*x
      y = abs(x)
      if (y>=bigs) then
        x  = x*bigi
        ix = ix + 1
      else if (y<bigsi) then
        x  = x*big
        ix = ix - 1
      end if
      ps(m)  = x
      ips(m) = ix
    end do

  end subroutine alfsx

  subroutine alfmp(t,m,an,bn,c,pn)
    use xreal_module, only: &
      xreal_type, assignment(=), operator(*), xreal_fxpgy

    real(kind=dp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(:), intent(in) :: an, bn
    real(kind=dp), intent(in) :: c
    type(xreal_type), dimension(0:) :: pn

    integer(kind=i4b) :: n, nmax

    nmax = size(pn)-1
    if (m>nmax) then
      return
    endif

    pn(m+1) = c*t*pn(m)
    do n=m+2, mmax
      pn(n) = xreal_fxpgy(an(n)*t,pn(n-1),-bn(n),pn(n-2))
    end do

  end subroutine alfmp

  subroutine alfmx(t,m,an,bn,c,pn,ipn)
    use xreal_module, only: big=>xreal_big, bigi=>xreal_bigi, &
      bigs=>xreal_bigs, bigsi=>xreal_bigsi

    real(kind=dp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(:), intent(in) :: an, bn
    real(kind=dp), intent(in) :: c
    real(kind=dp), dimension(0:) :: pn
    integer(kind=i4b), dimension(0:) :: ipn

    integer(kind=i4b) :: n, nmax, ix, iy, iz, id
    real(kind=dp) :: x, y, z, w

    nmax = size(pn)-1
    if (m>nmax) then
      return
    endif

    y  = pn(m)
    iy = ipn(m)
    x  = c*t*y
    ix = iy
    w  = abs(x)
    if (w>=bigs) then
      x  = x*bigi
      ix = ix + 1
    else if (w<bigsi) then
      x  = x*big
      ix = ix - 1
    end if
    pn(m+1)  = x
    ipn(m+1) = ix
    do n=m+2, mmax
      id = ix - iy
      select case(id)
        case(0)
          z  = (an(n)*t)*x-bn(n)*y
          iz = ix
        case(1)
          z  = (an(n)*t)*x-bn(n)*(y*bigi)
          iz = ix
        case(-1)
          z  = (an(n)*t)*(x*bigi)-bn(n)*y
          iz = iy
        case(2:)
          z  = (an(n)*t)*x
          iz = ix
        case(:-2)
          z  = -bn(n)*y
          iz = iy
      end select
      w = abs(z)
      if (w>bigs) then
        z = z*bigi
        iz = iz + 1
      else if (w<bigsi) then
        z = z*big
        iz = iz - 1
      end if
      pn(n)  = z
      ipn(n) = iz
      y  = x
      iy = ix
      x  = z
      ix = iz
    end do

  end subroutine alfmx

  subroutine alfx_test(ntrunc,nlat,un)
    use math_module, only: rad2deg=>math_rad2deg
    use xreal_module, only: xreal_type, xreal_base10, assignment(=)
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp), dimension(:), allocatable :: lat, wgt
    real(kind=dp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alfx_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    call alfx_init(ntrunc)
    call cpu_time(t1)
    call alfx_calc(lat)
    call cpu_time(t2)
    print *, "alfx_calc cpu time=", t2-t1
    call alfx_clean()
    call alfx_init(ntrunc)
    call cpu_time(t1)
    call alfx_calc_inline(lat)
    call cpu_time(t2)
    print *, "alfx_calc_inline cpu time=", t2-t1
    if (present(un)) then
      do j=1, nlat
        write(unit=un,rec=j) alfx_pnm(j,:,:)
      end do
    end if
    deallocate(lat,wgt)
    call alfx_clean()

  end subroutine alfx_test

  subroutine alfx_test_checksum(ntrunc,nlat,un)
    use kind_module, only: dp, i4b
    use xreal_module, only: xreal_type, big=>xreal_big, bigi=>xreal_bigi
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: alf_checksum, &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp) :: xx, dd, x, dx, p00
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn
    real(kind=dp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pmm, pnm, pj
    integer(kind=i4b), dimension(:), allocatable :: ipmm, ipnm
    real(kind=dp), dimension(:,:), allocatable ::  pjm, pjn
    integer(kind=i4b), dimension(:,:), allocatable :: ipjm, ipjn

    print *, "# ----- alfx_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:ntrunc),pnm(0:ntrunc),ipmm(0:ntrunc),ipnm(0:ntrunc))
    allocate(pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc), &
      ipjm(jmaxh,0:ntrunc),ipjn(jmaxh,0:ntrunc))
    allocate(pj(1:jmaxh))
    pjm(:,:) = 0.0_dp

    xx = 1.0_dp
    dd = 0.0_dp
    dx = 0.0_dp
    nn = 0
    mm = 0

    p00 = sqrt(0.5_dp)
    call alfx_init(ntrunc)
    do j=1, jmaxh
      pmm(0) = p00
      ipmm(0) = 0
      call alfx_calcps(coslat(j),dm,pmm,ipmm)
      pjm(j,:) = pmm(:)
      ipjm(j,:) = ipmm(:)
    end do
    do m=0, ntrunc
      do j=1, jmaxh
        pnm(m) = pjm(j,m)
        ipnm(m) = ipjm(j,m)
        pjn(j,m) = pjm(j,m)
        ipjn(j,m) = ipjm(j,m)
        if (m<ntrunc) then
          call alfx_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm,ipnm)
          pjn(j,m+1:ntrunc) = pnm(m+1:ntrunc)
          ipjn(j,m+1:ntrunc) = ipnm(m+1:ntrunc)
        end if
      end do
      do n=m, ntrunc
        do j=1, jmaxh
          select case(ipjn(j,n))
            case(0)
              pj(j) = pjn(j,n)
            case(:-1)
              pj(j) = pjn(j,n)*bigi
            case(1:)
              pj(j) = pjn(j,n)*big
          end select
        end do
        x = alf_checksum(wgt,pj)
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

    deallocate(lat,sinlat,coslat,wgt,pmm,pnm,ipmm,ipnm,pjm,ipjm,pjn,pj)
    call alfx_clean()
    
  end subroutine alfx_test_checksum

end module alfx_module
