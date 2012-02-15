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

  real(kind=dp), public, dimension(:,:,:), allocatable :: alf_pnm
  real(kind=dp), public, dimension(:,:), allocatable :: alf_anm, alf_bnm
  real(kind=dp), public, dimension(:), allocatable :: alf_cm, alf_dm

  real(kind=dp), private, dimension(:), allocatable :: sinlat, coslat

  integer(kind=i4b), private :: mmax
  real(kind=dp), private :: pstart

  public :: alf_init, alf_clean, alf_calc, &
    alf_calcps, alf_calcpn, alf_checksum, alf_test

contains

  subroutine alf_init(ntrunc)
    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

    mmax = ntrunc

    allocate(alf_cm(0:mmax),alf_dm(mmax))
    do m=0, mmax-1
      alf_cm(m) = sqrt(2.0_dp*m+3.0_dp)
    end do
    do m=1, mmax
      alf_dm(m) = sqrt(1.0_dp + 0.5_dp/real(m,kind=dp))
    end do
    allocate(alf_anm(mmax,0:mmax),alf_bnm(mmax,0:mmax))
    do m=0, mmax
      do n=m+2, mmax
        alf_anm(n,m) = sqrt((2.0_dp*n+1.0_dp)/real(n**2-m**2,kind=dp))
        alf_bnm(n,m) = alf_anm(n,m)*sqrt((n+m-1.0_dp)*(n-m-1.0_dp)/(2.0_dp*n-3.0_dp))
        alf_anm(n,m) = alf_anm(n,m)*sqrt(2.0_dp*n-1.0_dp)
      end do
    end do

  end subroutine alf_init

  subroutine alf_clean()

    deallocate(alf_anm,alf_bnm,alf_cm,alf_dm)
    if (allocated(alf_pnm)) then
      deallocate(alf_pnm)
    end if
    if (allocated(sinlat)) then
      deallocate(sinlat)
    end if
    if (allocated(coslat)) then
      deallocate(coslat)
    end if

  end subroutine alf_clean

  subroutine alf_calc(lat,p00)
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, jmax, jmaxh

    real(kind=dp), dimension(:), allocatable :: pmm, pnm

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

    allocate(alf_pnm(jmaxh, -1:mmax, 0:mmax))
    alf_pnm(:,:,:) = 0.0_dp
    allocate(sinlat(jmaxh),coslat(jmaxh))
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pmm(0:mmax),pnm(0:mmax))
    pmm(0) = pstart
    do j=1, jmaxh
      call alf_calcps(coslat(j),alf_dm,pmm)
      do m=0, mmax
        pnm(m) = pmm(m)
        alf_pnm(j,m,m) = pmm(m)
        call alf_calcpn(sinlat(j),m,alf_anm(:,m),alf_bnm(:,m),alf_cm(m),pnm)
        alf_pnm(j,m+1:mmax,m) = pnm(m+1:mmax)
      end do
    end do
    deallocate(pmm,pnm)

  end subroutine alf_calc

  subroutine alf_calcps(u,d,ps)

    real(kind=dp), intent(in) :: u ! coslat
    real(kind=dp), dimension(:), intent(in) :: d
    real(kind=dp), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
    end do

  end subroutine alf_calcps

  subroutine alf_calcpn(t,m,an,bn,c,pn)

    real(kind=dp), intent(in) :: t ! sinlat
    integer(kind=i4b), intent(in) :: m
    real(kind=dp), dimension(:), intent(in) :: an, bn
    real(kind=dp), intent(in) :: c
    real(kind=dp), dimension(0:), intent(inout) :: pn

    integer(kind=i4b) :: n, nmax

    nmax = size(pn)-1
    if (m>nmax) then
      return
    endif

    pn(m+1) = c*t*pn(m)
    do n=m+2, mmax
      pn(n) = an(n)*t*pn(n-1)-bn(n)*pn(n-2)
    end do

  end subroutine alf_calcpn

  function alf_checksum(wgt,pj) result(x)

    real(kind=dp), dimension(:), intent(in) :: wgt
    real(kind=dp), dimension(:), intent(in) :: pj

    real(kind=dp) :: x, d, dx
    integer(kind=i4b) :: jmaxh

    jmaxh = size(pj)
    x = 2.0_dp*sum(wgt(1:jmaxh)*pj(:)*pj(:))

  end function alf_checksum
  
  subroutine alf_test(nt_low,nt_high)
    use math_module, only: rad2deg=>math_rad2deg
    use glatwgt_module, only: glatwgt_calc

    integer(kind=i4b), intent(in), optional :: nt_low, nt_high

    integer(kind=i4b), parameter :: d = 10
    real(kind=dp), dimension(:), allocatable :: lat, wgt, pmm, pnm

    integer(kind=i4b) :: ntrunc_low = 159, ntrunc_high = 2159
    integer(kind=i4b) :: nlat, n, m, j, nn, mm, ntrunc
    real(kind=dp) :: x, dx, xx, dd, t1, t2

    print *, "# ----- alf_test() -----" 
    if (present(nt_low)) then
      ntrunc = nt_low
    else
      ntrunc = ntrunc_low
    end if
    nlat = (ntrunc+1)*3/2
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    call alf_init(ntrunc)
    call cpu_time(t1)
    call alf_calc(lat)
    call cpu_time(t2)
    print *, "alf_calc cpu time=", t2-t1
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(1)=", lat(1)*rad2deg 
    print *, "m, pm(n=m,j=1), pm(n=ntrunc,j=1)"
    do m=0, ntrunc, ntrunc/d
      print *, m, alf_pnm(1,m,m), alf_pnm(1,ntrunc,m)
    end do
    do m=ntrunc-1, ntrunc
      print *, m, alf_pnm(1,m,m), alf_pnm(1,ntrunc,m)
    end do
    print *, "n, pn(m=0,j=1), pn(m=1,j=1)"
    do n=0, ntrunc, ntrunc/d 
      print *, n, alf_pnm(1,n,0), alf_pnm(1,n,1)
    end do
    do n=ntrunc-1, ntrunc
      print *, n, alf_pnm(1,n,0), alf_pnm(1,n,1)
    end do
    print *, "x=\int pnm pnm dx error"
    xx = 1.0_dp
    dd = 0.0_dp
    nn = 0
    mm = 0
    do m=0, ntrunc, ntrunc
      do n=m, ntrunc, ntrunc/d
        x = alf_checksum(wgt,alf_pnm(:,n,m))
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
    call alf_clean()

    if (present(nt_high)) then
      ntrunc = nt_high
    else
      ntrunc = ntrunc_high
    end if
    if (ntrunc<0) then
      return
    end if
    allocate(pmm(0:ntrunc),pnm(0:ntrunc))
    pmm(0) = pstart
    nlat = (ntrunc+1)*3/2
    j = nlat/10
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(j=", j, ")=", lat(j)*rad2deg 
    call alf_init(ntrunc)
    call alf_calcps(cos(lat(j)),alf_dm,pmm)
    m = ntrunc/3
    pnm(m) = pmm(m)
    call alf_calcpn(sin(lat(j)),m,alf_anm(:,m),alf_bnm(:,m),alf_cm(m),pnm)
    print *, "n, m, pn(m=", m, ",j=",j,")"
    do n=m, ntrunc, ntrunc/d
      print *, n, m, pnm(n)
    end do
    do n=ntrunc-1, ntrunc
      print *, n, m, pnm(n)
    end do
    deallocate(pmm,pnm,lat,wgt)
    call alf_clean()
    print *, "# --------------------" 

  end subroutine alf_test

end module alf_module
