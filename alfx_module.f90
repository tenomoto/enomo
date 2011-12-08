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
  real(kind=dp), private :: pstart = sqrt(0.5_dp)

  real(kind=dp), private, dimension(:), allocatable :: sinlat, coslat
  integer(kind=i4b), private :: mmax

  public :: alfx_init, alfx_clean, alfx_calc, &
    alfx_calcps, alfx_calcpn, alfx_test

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

  subroutine alfx_calcps(u,d,ps)
    use xreal_module, only: xreal_type, assignment(=), operator(*)

    real(kind=dp), intent(in) :: u ! coslat
    real(kind=dp), dimension(:), intent(in) :: d
    type(xreal_type), dimension(0:), intent(inout) :: ps

    integer(kind=i4b) :: m, nmax

    nmax = size(ps)-1
    do m=1, nmax
      ps(m) = (d(m)*u)*ps(m-1)
    end do

  end subroutine alfx_calcps

  subroutine alfx_calcpn(t,m,an,bn,c,pn)
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

  end subroutine alfx_calcpn

  subroutine alfx_test(nt_low,nt_high)
    use math_module, only: rad2deg=>math_rad2deg
    use xreal_module, only: xreal_type, xreal_base10, assignment(=)
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: &
      anm=>alf_anm, bnm=>alf_bnm, cm=>alf_cm, dm=>alf_dm, alf_checksum

    integer(kind=i4b), intent(in), optional :: nt_low, nt_high

    integer(kind=i4b), parameter :: d = 10
    real(kind=dp), dimension(:), allocatable :: lat, wgt
    type(xreal_type), dimension(:), allocatable :: pmm, pnm

    integer(kind=i4b) :: ntrunc_low = 159, ntrunc_high = 2159
    integer(kind=i4b) :: nlat, n, m, j, jmaxh, nn, mm, ntrunc
    real(kind=dp) :: x, dx, xx, dd

    print *, "# ----- alfx_test() -----" 
    if (present(nt_low)) then
      ntrunc = nt_low
    else
      ntrunc = ntrunc_low
    end if
    nlat = (ntrunc+1)*3/2
    allocate(lat(nlat),wgt(nlat))
    jmaxh = size(alfx_pnm,1)
    call glatwgt_calc(lat,wgt,nlat)
    call alfx_init(ntrunc)
    call alfx_calc(lat)
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(1)=", lat(1)*rad2deg 
    print *, "m, pm(n=m,j=1), pm(n=ntrunc,j=1)"
    do m=0, ntrunc, ntrunc/d
      print *, m, alfx_pnm(1,m,m), alfx_pnm(1,ntrunc,m)
    end do
    do m=ntrunc-1, ntrunc
      print *, m, alfx_pnm(1,m,m), alfx_pnm(1,ntrunc,m)
    end do
    print *, "n, pn(m=0,j=1), pn(m=1,j=1)"
    do n=0, ntrunc, ntrunc/d 
      print *, n, alfx_pnm(1,n,0), alfx_pnm(1,n,1)
    end do
    do n=ntrunc-1, ntrunc
      print *, n, alfx_pnm(1,n,0), alfx_pnm(1,n,1)
    end do
    print *, "x=\int pnm pnm dx error"
    xx = 1.0_dp
    dd = 0.0_dp
    nn = 0
    mm = 0
    do m=0, ntrunc, ntrunc
      do n=m, ntrunc, ntrunc/d
        x = alf_checksum(wgt(1:jmaxh),alfx_pnm(:,n,m))
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
    call alfx_clean()

    if (present(nt_high)) then
      ntrunc = nt_high
    else
      ntrunc = ntrunc_high
    end if
    allocate(pmm(0:ntrunc),pnm(0:ntrunc))
    pmm(0) = pstart
    nlat = (ntrunc+1)*3/2
    j = nlat/10
    allocate(lat(nlat),wgt(nlat))
    call glatwgt_calc(lat,wgt,nlat)
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    print *, "lat(j=", j, ")=", lat(j)*rad2deg 
    call alfx_init(ntrunc)
    call alfx_calcps(cos(lat(j)),dm,pmm)
    m = ntrunc/3
    pnm(m) = pmm(m)
    call alfx_calcpn(sin(lat(j)),m,anm(:,m),bnm(:,m),cm(m),pnm)
    print *, "n, m, pn(m=", m, ",j=",j,")"
    do n=m, ntrunc, ntrunc/d
      print *, n, m, xreal_base10(pnm(n))
    end do
    do n=ntrunc-1, ntrunc
      print *, n, m, xreal_base10(pnm(n))
    end do
    deallocate(pmm,pnm,lat,wgt)
    print *, "# --------------------" 
    call alfx_clean()

  end subroutine alfx_test

end module alfx_module
