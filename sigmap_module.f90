module sigmap_module
! hybrid sigma-p coordinates
  use kind_module, only: i4b, dp
  implicit none
  private

! refences:
! Eckermann (2009) MWR
! Simmons and Burridge (1981)

  logical, public :: sigmap_verbose = .false.

  real(kind=dp), dimension(:), allocatable, private :: &
    pref, heta, a, b, r, dpref
  real(kind=dp), private :: ptop
  integer(kind=i4b), private :: nz

  public :: sigmap_init, sigmap_clean, sigmap_ps2p
  private :: bfunc

contains

  subroutine sigmap_init(h,rp,rs,s,pt,p0,kp,ks)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: h
    real(kind=dp), intent(in) :: rp, rs, s, pt, p0
    integer(kind=i4b), intent(in) :: kp, ks

    integer :: k
    real(kind=dp) :: atansr

    nz = size(h) - 1
    allocate(pref(0:nz),dpref(nz),heta(0:nz),a(0:nz),b(0:nz),r(nz))
    ptop = pt
    heta(:) = h
    pref(:) = ptop + heta(:)*(p0-pt)
    do k=1, nz
      dpref(k) = pref(k) - pref(k-1)
    end do

! boundary conditions
    a(0) = ptop
    b(0) = 0.0_dp
    a(nz) = ptop
    b(nz) = 1.0_dp

! set r
    atansr = 1.0_dp/atan(s)
    do k=1, kp
      r(k) = rp
    end do
    do k=kp+1, nz
      r(k) = rp + (rs-rp)*atansr*atan(s*bfunc(heta(k),heta(kp)))
    end do

! pure p levels
    do k=1, kp
      a(k) = pref(k)
      b(k) = 0.0_dp
    end do

! transition levels
    do k=kp+1, nz-ks
      b(k) = (bfunc(heta(k),heta(kp)))**r(k)
      a(k) = ptop + (heta(k)-b(k))*(pref(nz)-ptop)
    end do

! pure sigma levels
    do k=nz-ks+1, nz
      b(k) = bfunc(heta(k),heta(kp))
      a(k) = ptop + (heta(k)-b(k))*(pref(nz)-ptop)
    end do

    if (sigmap_verbose) then
      print *, "rs=", rs, " rp=", rp, " s=", s
      print *, "ks=", ks, " kp=", kp
      print *, "p0=", p0, " p(nz-ks+1=",nz-ks+1,")=",pref(nz-ks+1), &
               " p(kp=",kp,")=",pref(kp)," ptop=", pt
      print *, "k heta pref dpref a/p0 b r"
      print *, 0, heta(0), pref(0), "N/A", a(0)/p0, b(0), "N/A"
      do k=1, nz
        print *, k,heta(k),pref(k),dpref(k),a(k)/p0,b(k),r(k)
      end do
    end if

  end subroutine sigmap_init

  subroutine sigmap_clean()
    implicit none

    deallocate(pref,dpref,heta,a,b,r)

  end subroutine sigmap_clean

  subroutine sigmap_ps2p(ps,p)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: ps
    real(kind=dp), dimension(:,:,0:), intent(inout) :: p

    integer(kind=i4b) :: k, nz

    nz = size(p,3)-1
    do k=0, nz
      p(:,:,k) = a(k) + b(k)*(ps(:,:)-ptop)
    end do

  end subroutine sigmap_ps2p

  function bfunc(eta,etakp) result(bf)
    implicit none

    real(kind=dp), intent(in) :: eta, etakp
    real(kind=dp) :: bf

    bf = (eta-etakp)/(1.0_dp-etakp)

  end function bfunc

end module sigmap_module
