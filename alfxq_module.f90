module alfxq_module
! calculate associated Legendre functions
! in extended exponent floating numbers
  use kind_module, only: qp, i4b
  use xreal_module, only: xreal_quad_type
  implicit none
  private

! Source: Based on Fukushima (2011)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  normalised to 1 by default. factor (-1)**m is not included.

  real(kind=qp), private :: pstart

  integer(kind=i4b), private :: alfxq_ntrunc

  interface alfxq_calcps
    module procedure alfsp, alfsx
  end interface

  interface alfxq_calcpn
    module procedure alfmp, alfmx
  end interface

  public :: alfxq_init, alfxq_clean, alfxq_calc, alfxq_calc_m, alfxq_calc_inline, &
    alfxq_calc_inline_m, alfxq_calcps, alfxq_calcpn, alfxq_test, alfxq_test_checksum
  private :: alfsp, alfmp, alfsx, alfmx

contains

  subroutine alfxq_init(ntrunc)
    use alfq_module, only: alfq_init
    integer(kind=i4b), intent(in) :: ntrunc

!    print *, "alfxq_init"
    alfxq_ntrunc = ntrunc

    call alfq_init(ntrunc)

  end subroutine alfxq_init

  subroutine alfxq_clean()
    use alfq_module, only: alfq_clean

!    print *, "alfxq_clean"
    call alfq_clean()

  end subroutine alfxq_clean

  subroutine alfxq_calc(lat,alf,p00)
    use xreal_module, only: xreal_quad_type, assignment(=)
    use alfq_module, only: &
      anm=>alfq_anm, bnm=>alfq_bnm, cm=>alfq_cm, dm=>alfq_dm
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,0:,:), intent(out) :: alf
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, mmax

    type(xreal_quad_type), dimension(0:size(alf,1)-1) :: pmm, pnm
    real(kind=qp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alf,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_qp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alf(:,:,:) = 0.0_qp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0) = pstart
    do j=1, min(jmax,size(alf,3))
      call alfxq_calcps(coslat(j),dm,pmm)
      do m=0, mmax-1
        pnm(m) = pmm(m)
        alf(m,m,j) = pmm(m)
        call alfxq_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm)
        do n=m+1, mmax
          alf(n,m,j) = pnm(n)
        end do
      end do
      alf(mmax,mmax,j) = pmm(mmax)
    end do

  end subroutine alfxq_calc

  subroutine alfxq_calc_m(m,lat,alfm,p00)
    use xreal_module, only: xreal_quad_type, assignment(=)
    use alfq_module, only: &
      anm=>alfq_anm, bnm=>alfq_bnm, cm=>alfq_cm, dm=>alfq_dm
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,:), intent(out) :: alfm
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, n, jmax, mmax

    type(xreal_quad_type), dimension(0:size(alfm,1)-1) :: pmm, pnm
    real(kind=qp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alfm,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_qp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alfm(:,:) = 0.0_qp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0) = pstart
    do j=1, min(jmax,size(alfm,2))
      call alfxq_calcps(coslat(j),dm,pmm)
      if (m/=mmax) then
        pnm(m) = pmm(m)
        alfm(m,j) = pmm(m)
        call alfxq_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm)
        do n=m+1, mmax
          alfm(n,j) = pnm(n)
        end do
      else
        alfm(mmax,j) = pmm(mmax)
      end if
    end do

  end subroutine alfxq_calc_m

  subroutine alfxq_calc_inline(lat,alf,p00)
    use xreal_module, only: big=>xreal_bigq, bigi=>xreal_bigqi
    use alfq_module, only: &
      anm=>alfq_anm, bnm=>alfq_bnm, cm=>alfq_cm, dm=>alfq_dm
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,0:,:), intent(out) :: alf
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh, mmax

    real(kind=qp), dimension(0:size(alf,1)-1) :: pmm, pnm
    integer(kind=i4b), dimension(0:size(alf,1)-1) :: ipmm, ipnm
    real(kind=qp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alf,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_qp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alf(:,:,:) = 0.0_qp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0)  = pstart
    ipmm(0) = 0
    do j=1, min(jmax,size(alf,3))
      call alfxq_calcps(coslat(j),dm,pmm,ipmm)
      do m=0, mmax-1
        pnm(m)  = pmm(m)
        ipnm(m) = ipmm(m)
        select case (ipmm(m))
          case(0)
            alf(m,m,j) = pmm(m)
          case(:-1)
            alf(m,m,j) = pmm(m)*bigi
          case(1:)
            alf(m,m,j) = pmm(m)*big
        end select
        call alfxq_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm,ipnm)
        do n=m+1, mmax
          select case (ipnm(n))
            case(0)
              alf(n,m,j) = pnm(n)
            case(:-1)
              alf(n,m,j) = pnm(n)*bigi
            case(1:)
              alf(n,m,j) = pnm(n)*big
          end select
        end do
      end do
      select case (ipmm(mmax))
        case(0)
          alf(mmax,mmax,j) = pmm(mmax)
        case(:-1)
          alf(mmax,mmax,j) = pmm(mmax)*bigi
        case(1:)
          alf(mmax,mmax,j) = pmm(mmax)*big
      end select
   end do

  end subroutine alfxq_calc_inline

  subroutine alfxq_calc_inline_m(m,lat,alfm,p00)
    use xreal_module, only: big=>xreal_bigq, bigi=>xreal_bigqi
    use alfq_module, only: &
      anm=>alfq_anm, bnm=>alfq_bnm, cm=>alfq_cm, dm=>alfq_dm
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,:), intent(out) :: alfm
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, n, jmax, jmaxh, mmax

    real(kind=qp), dimension(0:size(alfm,1)-1) :: pmm, pnm
    integer(kind=i4b), dimension(0:size(alfm,1)-1) :: ipmm, ipnm
    real(kind=qp), dimension(size(lat)) :: sinlat, coslat

    mmax =  size(alfm,1) - 1
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_qp)
    end if
    jmax = size(lat)
    if (jmax<1) then
      return
    end if

    alfm(:,:) = 0.0_qp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))
    pmm(0)  = pstart
    ipmm(0) = 0
    do j=1, min(jmax,size(alfm,2))
      call alfxq_calcps(coslat(j),dm,pmm,ipmm)
      if (m/=mmax) then
        pnm(m)  = pmm(m)
        ipnm(m) = ipmm(m)
        select case (ipmm(m))
          case(0)
            alfm(m,j) = pmm(m)
          case(:-1)
            alfm(m,j) = pmm(m)*bigi
          case(1:)
            alfm(m,j) = pmm(m)*big
        end select
        call alfxq_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm,ipnm)
        do n=m+1, mmax
          select case (ipnm(n))
            case(0)
              alfm(n,j) = pnm(n)
            case(:-1)
              alfm(n,j) = pnm(n)*bigi
            case(1:)
              alfm(n,j) = pnm(n)*big
          end select
        end do
      else
        select case (ipmm(mmax))
          case(0)
            alfm(mmax,j) = pmm(mmax)
          case(:-1)
            alfm(mmax,j) = pmm(mmax)*bigi
          case(1:)
            alfm(mmax,j) = pmm(mmax)*big
        end select
     end if
   end do

  end subroutine alfxq_calc_inline_m

  subroutine alfsp(u,d,ps)
    use xreal_module, only: xreal_quad_type, assignment(=), operator(*)

    real(kind=qp), intent(in) :: u ! coslat
    real(kind=qp), dimension(:), intent(in) :: d
    type(xreal_quad_type), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
    end do

  end subroutine alfsp

  subroutine alfsx(u,d,ps,ips)
    use xreal_module, only: big=>xreal_bigq, bigi=>xreal_bigqi, &
      bigs=>xreal_bigqs, bigsi=>xreal_bigqsi

    real(kind=qp), intent(in) :: u ! coslat
    real(kind=qp), dimension(:), intent(in) :: d
    real(kind=qp), dimension(0:), intent(inout) :: ps
    integer(kind=i4b), dimension(0:), intent(inout) :: ips

    integer(kind=i4b) :: m, nmax, ix
    real(kind=qp) :: x, y

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
      xreal_quad_type, assignment(=), operator(*), xreal_fxpgyq

    real(kind=qp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: an, bn
    real(kind=qp), intent(in) :: c
    type(xreal_quad_type), dimension(0:) :: pn

    integer(kind=i4b) :: n, nmax

    nmax = size(pn)-1
    if (m+1>nmax) then
      return
    endif

    pn(m+1) = c*t*pn(m)
    do n=m+2, nmax
      pn(n) = xreal_fxpgyq(an(n)*t,pn(n-1),-bn(n),pn(n-2))
    end do

  end subroutine alfmp

  subroutine alfmx(t,m,an,bn,c,pn,ipn)
    use xreal_module, only: big=>xreal_bigq, bigi=>xreal_bigqi, &
      bigs=>xreal_bigqs, bigsi=>xreal_bigqsi

    real(kind=qp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: an, bn
    real(kind=qp), intent(in) :: c
    real(kind=qp), dimension(0:) :: pn
    integer(kind=i4b), dimension(0:) :: ipn

    integer(kind=i4b) :: n, nmax, ix, iy, iz, id
    real(kind=qp) :: x, y, z, w
    
    nmax = size(pn)-1
    if (m+1>nmax) then
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
    do n=m+2, nmax
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
      if (w>=bigs) then
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

  subroutine alfxq_test(ntrunc,nlat,un)
    use glatwgtq_module, only: glatwgtq_calc

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=qp), dimension(:), allocatable :: lat, wgt
    real(kind=qp), dimension(:,:,:), allocatable :: alf
    real(kind=qp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alfxq_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    allocate(alf(0:ntrunc,0:ntrunc,nlat/2))
    call glatwgtq_calc(lat,wgt)
    call alfxq_init(ntrunc)
    call cpu_time(t1)
    call alfxq_calc(lat(1:nlat/2),alf)
    call cpu_time(t2)
    print *, "alfxq_calc cpu time=", t2-t1
    call alfxq_clean()
    call alfxq_init(ntrunc)
    call cpu_time(t1)
    call alfxq_calc_inline(lat(1:nlat/2),alf)
    call cpu_time(t2)
    print *, "alfxq_calc_inline cpu time=", t2-t1
    if (present(un)) then
      do j=1, nlat
        write(unit=un,rec=1) alf
      end do
    end if
    deallocate(alf)
    deallocate(lat,wgt)
    call alfxq_clean()

  end subroutine alfxq_test

  subroutine alfxq_test_checksum(ntrunc,nlat,un)
    use kind_module, only: qp, i4b
    use xreal_module, only: xreal_quad_type, big=>xreal_bigq, bigi=>xreal_bigqi
    use glatwgtq_module, only: glatwgtq_calc
    use alfq_module, only: alfq_checksum, &
      anm=>alfq_anm, bnm=>alfq_bnm, cm=>alfq_cm, dm=>alfq_dm

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=qp) :: xx, dd, x, dx, p00
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn
    real(kind=qp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pmm, pnm, pj
    integer(kind=i4b), dimension(:), allocatable :: ipmm, ipnm
    real(kind=qp), dimension(:,:), allocatable ::  pjm, pjn
    integer(kind=i4b), dimension(:,:), allocatable :: ipjm, ipjn

    print *, "# ----- alfxq_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgtq_calc(lat,wgt)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:ntrunc),pnm(0:ntrunc),ipmm(0:ntrunc),ipnm(0:ntrunc))
    allocate(pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc), &
      ipjm(jmaxh,0:ntrunc),ipjn(jmaxh,0:ntrunc))
    allocate(pj(1:jmaxh))
    pjm(:,:) = 0.0_qp

    xx = 1.0_qp
    dd = 0.0_qp
    dx = 0.0_qp
    nn = 0
    mm = 0

    p00 = sqrt(0.5_qp)
    call alfxq_init(ntrunc)
    do j=1, jmaxh
      pmm(0) = p00
      ipmm(0) = 0
      call alfxq_calcps(coslat(j),dm,pmm,ipmm)
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
          call alfxq_calcpn(sinlat(j),m,anm(:,m),bnm(:,m),cm(m),pnm,ipnm)
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
        x = alfq_checksum(wgt,pj)
        dx = 1.0_qp - x
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
    call alfxq_clean()
    
  end subroutine alfxq_test_checksum

end module alfxq_module
