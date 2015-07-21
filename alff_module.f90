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

  real(kind=dp), private, dimension(:,:), allocatable :: enm, fnm, gnm
  real(kind=dp), private, dimension(:,:), allocatable :: ank

  integer(kind=i4b), private :: mmax, jmax, jmaxh
  real(kind=dp), private :: pstart

  public :: alff_init, alff_clean, alff_calc, alff_calc_m, &
    alff_calcp0, alff_calcp1, alff_calcpn, alff_test, alff_test_checksum
  private :: fouriercoeff, gamma_poly

contains

  subroutine alff_init(ntrunc)
    use alf_module, only: alf_init

    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m

!    print *, "alff_init"
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

!    print *, "alff_clean"
    deallocate(ank,enm,fnm,gnm)

    call alf_clean()

  end subroutine alff_clean
  
  subroutine alff_calc(lat,alf,p00)
    use math_module, only: pih=>math_pih
    use integer_module, only: swap=>integer_swap

    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), dimension(0:,0:,:), intent(out) :: alf
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, k1, k2
    real(kind=dp) :: theta
    real(kind=dp), dimension(size(lat)) :: sinlat, coslat
    real(kind=dp), dimension(0:mmax,2) :: pn
    
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    call fouriercoeff(pstart)
    jmax = size(lat)
    alf(:,:,:) = 0.0_dp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))

! calculate Pnm
    do j=1, min(jmax,size(alf,3))

      theta = pih-lat(j) ! lat => colat

      pn(:,:) = 0.0_dp
      pn(0,1) = pstart ! m = 0
      pn(1,1) = sqrt(3.0_dp)*sinlat(j)*pn(0,1)
      call alff_calcp0(theta,pn(:,1))
      alf(0:mmax,0,j) = pn(0:mmax,1)
      k1 = 1
      k2 = 2
      do m=2, mmax, 2 ! m even
        call alff_calcpn(m,pn(:,k1),pn(:,k2))
        alf(m:mmax,m,j) = pn(m:mmax,k2)
        call swap(k1,k2)
      end do

      pn(:,:) = 0.0_dp
      pn(1,1) = sqrt(1.5_dp)*coslat(j)*pstart ! m = 1
      pn(2,1) = sqrt(5.0_dp)*sinlat(j)*pn(1,1)
      call alff_calcp1(theta,pn(:,1))
      alf(1:mmax,1,j) = pn(1:mmax,1)
      k1 = 1
      k2 = 2
      do m=3, mmax, 2 ! m odd
        call alff_calcpn(m,pn(:,k1),pn(:,k2))
        alf(m:mmax,m,j) = pn(m:mmax,k2)
        call swap(k1,k2)
      end do

    end do! j

  end subroutine alff_calc

  subroutine alff_calc_m(mout,lat,alfm,p00)
    use math_module, only: pih=>math_pih
    use integer_module, only: swap=>integer_swap

    integer(kind=i4b), intent(in) :: mout 
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), dimension(0:,:), intent(out) :: alfm
    real(kind=dp), intent(in), optional :: p00

    integer(kind=i4b) :: j, m, n, k1, k2
    real(kind=dp) :: theta
    real(kind=dp), dimension(size(lat)) :: sinlat, coslat
    real(kind=dp), dimension(0:mmax) :: pmm
    real(kind=dp), dimension(0:mmax,2) :: pn
    
    if (present(p00)) then
      pstart = p00
    else
      pstart = sqrt(0.5_dp)
    end if
    call fouriercoeff(pstart)
    jmax = size(lat)
    alfm(:,:) = 0.0_dp
    sinlat(:) = sin(lat(:))
    coslat(:) = cos(lat(:))

    if (modulo(mout,2) == 0) then ! even

! calculate Pnm
      do j=1, min(jmax,size(alfm,2))

        theta = pih-lat(j) ! lat => colat

        pn(:,:) = 0.0_dp
        pn(0,1) = pstart ! m = 0
        pn(1,1) = sqrt(3.0_dp)*sinlat(j)*pn(0,1)
        call alff_calcp0(theta,pn(:,1))
        if (mout==0) then
          alfm(0:mmax,j) = pn(0:mmax,1)
        else
          k1 = 1
          k2 = 2
          do m=2, mout, 2 ! m even
            call alff_calcpn(m,pn(:,k1),pn(:,k2))
            alfm(m:mmax,j) = pn(m:mmax,k2)
            call swap(k1,k2)
          end do
        end if
      end do! j

    else ! odd

      do j=1, jmax

        theta = pih-lat(j) ! lat => colat

        pn(:,:) = 0.0_dp
        pn(1,1) = sqrt(1.5_dp)*coslat(j)*pstart ! m = 1
        pn(2,1) = sqrt(5.0_dp)*sinlat(j)*pn(1,1)
        call alff_calcp1(theta,pn(:,1))
        if (m==1) then
          alfm (1:mmax,j) = pn(1:mmax,1)
        else
          k1 = 1
          k2 = 2
          do m=3, mout, 2 ! m odd
            call alff_calcpn(m,pn(:,k1),pn(:,k2))
            alfm(m:mmax,j) = pn(m:mmax,k2)
            call swap(k1,k2)
          end do
        end if

      end do! j

    end if ! even or odd

  end subroutine alff_calc_m

  subroutine fouriercoeff(p0)
    use math_module, only:  pir => math_pir
    real(kind=dp), intent(in) :: p0

    integer(kind=i4b), parameter :: nc = 128

    integer(kind=i4b) :: n, k, n2, nh, l1
    real(kind=dp) :: y

! calculate fourier coefficients for Pn
    ank(2,1) = 0.75_dp*sqrt(5.0_dp) * p0
!    do n=3, mmax
    do n=3, min(nc, mmax)
!      n2 = n * 2
!      ank(n,n/2) = ank(n-1,(n-1)/2) * &
!        sqrt((n2 - 1.0_dp) * (n2 + 1.0_dp)) / n2
        ank(n,n/2) = ank(n-1,(n-1)/2) * &
          sqrt(1.0_dp - 0.25_dp/(n*n))
    end do
    do n=nc+1, mmax
      y = 1.0_dp / n
      ank(n, n/2) = sqrt(2.0_dp * (2.0_dp * n + 1.0_dp) * pir * y) * &
        gamma_poly(0.5_dp * y) / gamma_poly(y) ** 2
    end do
    do n=2, mmax
      nh  = n/2
      do k=1, nh
        l1 = 2 * k - 1
        ank(n,nh-k) = ank(n,nh-k+1) * &
          l1 * (n - k + 1.0_dp) / (k * (2.0_dp*n - l1))
      end do
      if (n==nh*2) then
        ank(n,0) = 0.5_dp*ank(n,0)
      end if
    end do

  end subroutine fouriercoeff

  function gamma_poly(y) result (g)
    real(kind=dp), intent(in) :: y
    real(kind=dp) :: g

    integer(kind=i4b), parameter :: n = 7
    real(kind=dp), dimension(n), parameter :: &
      c = (/1.0_dp/12.0_dp, 1.0_dp/288.0_dp, -139.0_dp/51840.0_dp, -571.0_dp/2488320.0_dp, &
            163879.0_dp/209018880.0_dp, 5246819.0_dp/75246796800.0_dp, -534703531.0_dp/902961561600.0_dp/)
    integer(kind=i4b) :: i

    g = 0.0_dp

    do i=n, 1, -1
      g = (g + c(i)) * y
    end do
    g = g + 1.0_dp

  end function gamma_poly

  subroutine alff_calcp0(theta,p0)
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(0:), intent(inout) :: p0

    integer(kind=i4b) :: n, l, nmax
    real(kind=dp) :: c1, s1, c2, s2, c, s, ct

    nmax = size(p0)-1
    c1 = cos(theta)
    s1 = sin(theta)
    c2 = c1*c1 - s1*s1
    s2 = 2.0_dp * c1 * s1
    do n=2, nmax
      p0(n) = 0.0_dp
      if (mod(n,2) == 0) then
        c = 1.0_dp
        s = 0.0_dp
      else
        c = c1
        s = s1
      end if
      do l=0, n/2 
        p0(n) = p0(n) + ank(n,l) * c
        ct = c2 * c - s2 * s
        s  = s2 * c + c2 * s
        c = ct
      end do
    end do

  end subroutine alff_calcp0

  subroutine alff_calcp1(theta,p1)
    real(kind=dp), intent(in) :: theta
    real(kind=dp), dimension(0:), intent(inout) :: p1

    integer(kind=i4b) :: n, l, k, nmod, nmax
    real(kind=dp) :: sqrtnnr, c1, s1, c2, s2, c, s, ct

    nmax = size(p1)-1
    c1 = cos(theta)
    s1 = sin(theta)
    c2 = c1*c1 - s1*s1
    s2 = 2.0_dp * c1 * s1
    do n=3, nmax
      p1(n) = 0.0_dp
      sqrtnnr = 1.0_dp/sqrt(n*(n+1.0_dp))
      if (mod(n,2) == 0) then
        c = 1.0_dp
        s = 0.0_dp
        nmod = 0
      else
        c = c1
        s = s1
        nmod = 1
      end if
      do l=0, n/2
        k = 2*l + nmod ! n even: k=2*l, n odd: k=2*l+1
        p1(n) = p1(n) + ank(n,l) * (2*l+nmod) * sqrtnnr * s
        ct = c2 * c - s2 * s
        s  = s2 * c + c2 * s
        c = ct
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
    do n=m, nmax
!    do n=m+2, nmax
      pn1(n) = enm(n,m)*pn0(n-2) + fnm(n,m)*pn1(n-2) - gnm(n,m)*pn0(n)
    end do

  end subroutine alff_calcpn

  subroutine alff_test(ntrunc,nlat,un)
    use math_module, only: rad2deg=>math_rad2deg
    use glatwgt_module, only: glatwgt_calc

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp), dimension(:), allocatable :: lat, wgt
    real(kind=dp), dimension(:,:,:), allocatable :: alf

    real(kind=dp) :: t1, t2
    integer(kind=i4b) :: j

    print *, "# ----- alff_test() -----" 
    print *, "ntrunc=", ntrunc, " nlat=", nlat
    allocate(lat(nlat),wgt(nlat))
    allocate(alf(0:ntrunc,0:ntrunc,nlat/2))
    call glatwgt_calc(lat,wgt)
    call alff_init(ntrunc)
    call cpu_time(t1)
    call alff_calc(lat(1:nlat/2),alf)
    call cpu_time(t2)
    print *, "alff_calc cpu time=", t2-t1
    if (present(un)) then
      write(unit=un,rec=1) alf
    end if
    call alff_clean()

  end subroutine alff_test

  subroutine alff_test_checksum(ntrunc,nlat,un)
    use kind_module, only: dp, i4b
    use math_module, only: pih=>math_pih
    use integer_module, only: swap=>integer_swap
    use glatwgt_module, only: glatwgt_calc
    use alf_module, only: alf_checksum

    implicit none

    integer(kind=i4b), intent(in) :: ntrunc, nlat
    integer(kind=i4b), intent(in), optional :: un

    real(kind=dp) :: xx, dd, x, dx, p00, theta
    integer(kind=i4b) :: jmaxh, m, n, j, mm, nn, k1, k2
    real(kind=dp), dimension(:), allocatable :: &
      lat, sinlat, coslat, wgt, pn
    real(kind=dp), dimension(:,:), allocatable ::  pjm
    real(kind=dp), dimension(:,:,:), allocatable ::  pjn

    print *, "# ----- alff_test_checksum() -----" 
    print *, "x=\int pnm pnm dx error"
    print *, "ntrunc=", ntrunc, " nlat=", nlat

    call alff_init(ntrunc)
    jmaxh = nlat/2
    allocate(lat(nlat),sinlat(jmaxh),coslat(jmaxh),wgt(nlat))
    call glatwgt_calc(lat,wgt)
    sinlat(:) = sin(lat(1:jmaxh))
    coslat(:) = cos(lat(1:jmaxh))
    allocate(pn(0:ntrunc), pjm(jmaxh,0:ntrunc),pjn(jmaxh,0:ntrunc,2))
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

! m = 0
    do j=1, jmaxh
      theta = pih-lat(j) ! lat => colat
      pn(0) = p00 ! m = 0
      pn(1) = sqrt(3.0_dp)*sinlat(j)*p00
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
      pn(1) = sqrt(1.5_dp)*coslat(j)*p00 ! m = 1
      pn(2) = sqrt(5.0_dp)*sinlat(j)*pn(1)
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
