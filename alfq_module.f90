module alfq_module
! calculates normalised associated Legendre functions
! using three point recurrence
  use kind_module, only: i4b, dp, qp
  implicit none
  private

! Source: Based on Fukushima (2011)
! Author: T. Enomoto
! Usage:
!   Calculates the values of normalized associated Legendre polynomials
!   at latitudes lat
! NB:
!  normalised to 1 by default. factor (-1)**m is not included.


  real(kind=qp), public, dimension(:,:), allocatable :: alfq_anm, alfq_bnm
  real(kind=qp), public, dimension(:), allocatable :: alfq_cm, alfq_dm
  
  real(kind=qp), private, parameter :: quad_min = 2.0_qp**(-16494)

  integer(kind=i4b), private :: alfq_ntrunc = 0
  real(kind=qp), private :: pstart
  integer, private :: retain = 0

  public :: alfq_init, alfq_clean, alfq_calc, alfq_calc_m, &
    alfq_calcps, alfq_calcpn, alfq_checksum, alfq_test, alfq_test_checksum

contains

  subroutine alfq_init(ntrunc)
    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

!    print *, "alfq_init"
    if (alfq_ntrunc==ntrunc) then
      retain = retain + 1
!      print *, retain
      return
    end if
    alfq_ntrunc = ntrunc
!    print *, "Allocating alfq_anm, alfq_bnm, alfq_cm, alfq_dm"

    allocate(alfq_cm(0:ntrunc),alfq_dm(ntrunc))
    do m=0, ntrunc-1
      alfq_cm(m) = sqrt(2.0_qp*m+3.0_qp)
    end do
    do m=1, ntrunc
      alfq_dm(m) = sqrt(1.0_qp + 0.5_qp/real(m,kind=qp))
    end do
    allocate(alfq_anm(ntrunc,0:ntrunc),alfq_bnm(ntrunc,0:ntrunc))
    do m=0, ntrunc
      do n=m+2, ntrunc
        alfq_anm(n,m) = sqrt((2.0_qp*n+1.0_qp)/real(n**2-m**2,kind=qp))
        alfq_bnm(n,m) = alfq_anm(n,m)*sqrt((n+m-1.0_qp)*(n-m-1.0_qp)/(2.0_qp*n-3.0_qp))
        alfq_anm(n,m) = alfq_anm(n,m)*sqrt(2.0_qp*n-1.0_qp)
      end do
    end do
    retain = retain + 1
!    print *, retain

  end subroutine alfq_init

  subroutine alfq_clean()

!    print *, "alfq_clean"
    retain = retain - 1 
!    print *, retain
    if (retain<1) then
!      print *, "Deallocating alfq_anm, alfq_bnm, alfq_cm, alfq_dm"
      deallocate(alfq_anm,alfq_bnm,alfq_cm,alfq_dm)
      alfq_ntrunc = 0
    end if

  end subroutine alfq_clean

  subroutine alfq_calc(lat,alf,p00)
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,0:,:), intent(out) :: alf
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh, mmax

    real(kind=qp), dimension(0:size(alf,1)-1) :: pmm, pnm
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
      call alfq_calcps(coslat(j),alfq_dm,pmm)
      do m=0, mmax-1
        pnm(m) = pmm(m)
        alf(m,m,j) = pmm(m)
        call alfq_calcpn(sinlat(j),m,alfq_anm(:,m),alfq_bnm(:,m),alfq_cm(m),pnm)
        alf(m+1:mmax,m,j) = pnm(m+1:mmax)
      end do
      alf(mmax,mmax,j) = pmm(mmax)
    end do

  end subroutine alfq_calc

  subroutine alfq_calc_m(m,lat,alfm,p00)
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: lat
    real(kind=qp), dimension(0:,:), intent(out) :: alfm
    real(kind=qp), intent(in), optional :: p00

    integer(kind=i4b) :: j, n, jmax, jmaxh, mmax

    real(kind=qp), dimension(0:size(alfm,1)-1) :: pmm, pnm
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
      call alfq_calcps(coslat(j),alfq_dm,pmm)
      if (m/=mmax) then
        pnm(m) = pmm(m)
        alfm(m,j) = pmm(m)
        call alfq_calcpn(sinlat(j),m,alfq_anm(:,m),alfq_bnm(:,m),alfq_cm(m),pnm)
        alfm(m+1:mmax,j) = pnm(m+1:mmax)
      else 
        alfm(mmax,j) = pmm(mmax)
      end if
    end do

  end subroutine alfq_calc_m

  subroutine alfq_calcps(u,d,ps)

    real(kind=qp), intent(in) :: u ! coslat
    real(kind=qp), dimension(:), intent(in) :: d
    real(kind=qp), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
      if (abs(ps(m))<quad_min) then
        exit
      end if
    end do
    ps(m:) = 0.0_qp

  end subroutine alfq_calcps

  subroutine alfq_calcpn(t,m,an,bn,c,pn)

    real(kind=qp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=qp), dimension(:), intent(in) :: an, bn
    real(kind=qp), intent(in) :: c
    real(kind=qp), dimension(0:), intent(inout) :: pn

    integer(kind=i4b) :: n, nmax

    nmax = size(pn)-1
    if (m+1>nmax) then
      return
    endif

    pn(m+1) = c*t*pn(m)
    do n=m+2, nmax
      pn(n) = an(n)*t*pn(n-1)-bn(n)*pn(n-2)
    end do

  end subroutine alfq_calcpn

  function alfq_checksum(wgt,pj) result(x)

    real(kind=qp), dimension(:), intent(in) :: wgt
    real(kind=qp), dimension(:), intent(in) :: pj

    real(kind=qp) :: x
    integer(kind=i4b) :: jmaxh

    jmaxh = size(pj)
    x = 2.0_qp*sum(wgt(1:jmaxh)*pj(:)*pj(:))

  end function alfq_checksum
  
  subroutine alfq_test(ntrunc,nlat,un)
    use math_module, only: rad2deg=>math_rad2degq
    use glatwgtq_module, only: glatwgtq_calc

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=qp), dimension(:), allocatable :: lat, wgt
    real(kind=qp), dimension(:,:,:), allocatable :: alf
    real(kind=qp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alfq_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    allocate(alf(0:ntrunc,0:ntrunc,nlat/2))
    call glatwgtq_calc(lat,wgt)
    call alfq_init(ntrunc)
    call cpu_time(t1)
    call alfq_calc(lat(1:nlat/2),alf)
    call cpu_time(t2)
    print *, "alfq_calc cpu time=", t2-t1
    if (present(un)) then
      write(unit=un,rec=1) alf
    end if
    deallocate(alf)
    deallocate(lat,wgt)
    call alfq_clean()
    
  end subroutine alfq_test

  subroutine alfq_test_checksum(ntrunc,nlat,un)
    use kind_module, only: qp, i4b
    use glatwgtq_module, only: glatwgtq_calc

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=qp) :: xx, dd, x, dx, p00
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn
    real(kind=qp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pmm, pnm
    real(kind=qp), dimension(:,:), allocatable ::  pjm, pjn

    print *, "# ----- alfq_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgtq_calc(lat,wgt)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:ntrunc),pnm(0:ntrunc), &
      pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc))

    xx = 1.0_qp
    dd = 0.0_qp
    dx = 0.0_qp
    nn = 0
    mm = 0

    p00 = sqrt(0.5_qp)
    call alfq_init(ntrunc)
    do j=1, jmaxh
      pmm(0) = p00
      call alfq_calcps(coslat(j),alfq_dm,pmm)
      pjm(j,:) = pmm(:)
    end do
    do m=0, ntrunc
      do j=1, jmaxh
        pnm(m) = pjm(j,m)
        pjn(j,m) = pjm(j,m)
        if (m<ntrunc) then
          call alfq_calcpn(sinlat(j),m,alfq_anm(:,m),alfq_bnm(:,m),alfq_cm(m),pnm)
          pjn(j,m+1:ntrunc) = pnm(m+1:ntrunc)
        end if
      end do
      do n=m, ntrunc
        x = alfq_checksum(wgt,pjn(:,n))
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

    deallocate(lat,sinlat,coslat,wgt,pmm,pnm,pjm,pjn)
    call alfq_clean()
    
  end subroutine alfq_test_checksum

end module alfq_module
