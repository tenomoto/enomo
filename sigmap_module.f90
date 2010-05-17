module sigmap_module
! hybrid sigma-p coordinates
  use kind_module, only: i4b, dp
  implicit none
  private

! refences:
! Eckermann (2009) MWR
! Ritchie et al. (1994) MWR
! Simmons and Burridge (1981) MWR

  logical, public :: sigmap_verbose = .false.

  real(kind=dp), dimension(:), allocatable, private :: &
    pref, etah, a, b, c, dpref, r, db
  real(kind=dp), private :: ptop
  real(kind=dp), dimension(:,:,:), allocatable, private :: dpf, dpfr, lnpp, alpha

  public :: sigmap_init, sigmap_clean, sigmap_ps2p, sigmap_geopotential, &
    sigmap_integrate, sigmap_advect, sigmap_semiimplicit
  private :: bfunc, setp

contains

  subroutine sigmap_init(h,rp,rs,s,pt,p0,kp,ks)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: h
    real(kind=dp), intent(in) :: rp, rs, s, pt, p0
    integer(kind=i4b), intent(in) :: kp, ks

    integer(kind=i4b) :: k, nz
    real(kind=dp) :: atansr

    nz = size(h) - 1
    allocate(pref(0:nz),etah(0:nz),a(0:nz),b(0:nz))
    allocate(dpref(nz),r(nz),db(nz),c(nz))
    ptop = pt
    etah(:) = h
    pref(:) = ptop + etah(:)*(p0-pt)

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
      r(k) = rp + (rs-rp)*atansr*atan(s*bfunc(etah(k),etah(kp)))
    end do

! pure p levels
    do k=1, kp
      a(k) = pref(k)
      b(k) = 0.0_dp
    end do

! transition levels
    do k=kp+1, nz-ks
      b(k) = (bfunc(etah(k),etah(kp)))**r(k)
      a(k) = ptop + (etah(k)-b(k))*(pref(nz)-ptop)
    end do

! pure sigma levels
    do k=nz-ks+1, nz
      b(k) = bfunc(etah(k),etah(kp))
      a(k) = ptop + (etah(k)-b(k))*(pref(nz)-ptop)
    end do

! set dpref, db, c, f
    dpref(:) = pref(1:nz)-pref(0:nz-1)
    db(:) = b(1:nz) - b(0:nz-1)
    c(:) = a(1:nz)*b(0:nz-1)-a(0:nz-1)*b(1:nz)

    if (sigmap_verbose) then
      print *, "rs=", rs, " rp=", rp, " s=", s
      print *, "ks=", ks, " kp=", kp
      print *, "p0=", p0, " p(nz-ks+1=",nz-ks+1,")=",pref(nz-ks+1), &
               " p(kp=",kp,")=",pref(kp)," ptop=", pt
      print *, "k etah pref a/p0 b"
      do k=0, nz
        print *, k,etah(k),pref(k),a(k)/p0,b(k)
      end do
      print *, "k db r c f"
      do k=1, nz
        print *, k,dpref(k),db(k),r(k),c(k)
      end do
    end if

  end subroutine sigmap_init

  subroutine sigmap_clean()
    implicit none

    deallocate(pref,etah,a,b,dpref,r,db,c)
    if (allocated(dpf)) then
      deallocate(dpf,dpfr,lnpp,alpha)
    end if

  end subroutine sigmap_clean

  subroutine sigmap_ps2p(ps,ph,pf)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: ps
    real(kind=dp), dimension(:,:,0:), intent(inout) :: ph
    real(kind=dp), dimension(:,:,:), intent(inout) :: pf

    integer(kind=i4b) :: k, nz

    nz = size(ph,3)-1
    do k=0, nz
      ph(:,:,k) = a(k) + b(k)*(ps(:,:)-ptop)
    end do
    pf(:,:,:) = 0.5_dp*(ph(:,:,0:nz-1)-ph(:,:,1:nz))

  end subroutine sigmap_ps2p

  subroutine sigmap_geopotential(T,gzs,gz)
    use air_module, only: rd=>air_rd
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: T
    real(kind=dp), dimension(:,:), intent(in) :: gzs
    real(kind=dp), dimension(:,:,:), intent(inout) :: gz 

    integer(kind=i4b) :: nx, ny, nz, i, j, k
    real(kind=dp) :: gzh

    nx = size(T,1)
    ny = size(T,2)
    nz = size(T,3)

    do j=1, ny
      do i=1, nx
        gzh = gzs(i,j)
        gz(i,j,nz) = gzh + alpha(i,j,nz)*rd*T(i,j,nz)
        do k=nz-1, 1, -1
          gzh = gzh + rd*T(i,j,k+1)*lnpp(i,j,k+1)
          gz(i,j,k) = gzh + alpha(i,j,k)*rd*T(i,j,k)
        end do
      end do
    end do

  end subroutine sigmap_geopotential

  subroutine sigmap_integrate(ps,dlnpsdx,dlnpsdy,d,u,v, &
    dlnpsdt,dlnpdx,dlnpdy,omegap,weta,ph,pf)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: ps, dlnpsdx, dlnpsdy
    real(kind=dp), dimension(:,:,:), intent(in) :: d, u, v
    real(kind=dp), dimension(:,:), intent(inout) :: dlnpsdt
    real(kind=dp), dimension(:,:,:), intent(inout) :: dlnpdx, dlnpdy, omegap, pf
! weta: \dot{eta}\frac{\partial p}{\partial\eta}
    real(kind=dp), dimension(0:,:,:), intent(inout) :: weta, ph

    integer(kind=i4b) :: i, j, k, nx, ny, nz
    real(kind=dp) :: f, dvdp

    nx = size(d,1)
    ny = size(d,2)
    nz = size(d,3)

    call sigmap_ps2p(ps,ph,pf)
    call setp(ph)

    dlnpsdt(:,:) = 0.0_dp
    do k=1, nz
      do j=1, ny
        do i=1, nx
          dvdp = d(i,j,k)/ps(i,j)*dpf(i,j,k)+ &
            db(k)*(u(i,j,k)*dlnpsdx(i,j)+v(i,j,k)*dlnpsdy(i,j))
! pressure gradient
          f = ps(i,j)*dpfr(i,j,k)*(db(k)+c(k)*dpfr(i,j,k)*lnpp(i,j,k))
          dlnpdx(i,j,k) = f*dlnpsdx(i,j)
          dlnpdy(i,j,k) = f*dlnpsdy(i,j)
! omega/p
          omegap(i,j,k) = -ps(i,j)*(lnpp(i,j,k)*dlnpsdt(i,j) + alpha(i,j,k)*dvdp &
            + u(i,j,k)*dlnpdx(i,j,k)+v(i,j,k)*dlnpdy(i,j,k))*dpfr(i,j,k)
! sum from 1 to k
          dlnpsdt(i,j) = dlnpsdt(i,j) + dvdp
          weta(i,j,k) = dlnpsdt(i,j)
        end do
      end do
    end do
    dlnpsdt(:,:) = -dlnpsdt(:,:)

! mass weighted vertical velocity
    do k=1, nz
      do j=1, ny
        do i=1, nx
          weta(i,j,k) = -ps(i,j)*(b(k)*dlnpsdt(i,j)+weta(i,j,k))
        end do
      end do
    end do

  end subroutine sigmap_integrate

  subroutine sigmap_advect(weta,x,vertadv)
    implicit none

    real(kind=dp), dimension(:,:,0:), intent(in) :: weta
    real(kind=dp), dimension(:,:,:), intent(in) :: x
    real(kind=dp), dimension(:,:,:), intent(inout) :: vertadv
    
    integer(kind=i4b) :: i, j, k, nx, ny, nz, km, kp
    real(kind=dp) :: f

    nx = size(x,1)
    ny = size(x,2)
    nz = size(x,3)

    do k=1, nz
      kp = min(k+1,nz)
      km = max(k-1,1)
      do j=1, ny
        do i=1, nx
          vertadv(i,j,k) = 0.5_dp*(weta(i,j,k)*(x(i,j,kp)-x(i,j,k)) + &
            weta(i,j,k-1)*(x(i,j,k)-x(i,j,km)))*dpfr(i,j,k)
        end do
      end do
    end do

  end subroutine sigmap_advect

  subroutine sigmap_semiimplicit(Tref,psref,bmatrix)
    use air_module, only: rd=>air_rd
    implicit none

    real(kind=dp), intent(in) :: Tref, psref
    real(kind=dp), dimension(:,:), intent(inout) :: bmatrix

    integer(kind=i4b) :: k, j, nz
    real(kind=dp), dimension(:), allocatable :: lnppref, alpharef, Trefv
    real(kind=dp), dimension(:,:), allocatable :: gmatrix

    nz = size(bmatrix,1)
    allocate(lnppref(nz),alpharef(nz),gmatrix(nz,nz),Trefv(nz))

    lnppref(1) = 2.0_dp
    lnppref(2:nz) = log(pref(2:nz)/pref(1:nz-1))
    alpharef(1) = log(2.0_dp)
    alpharef(2:) = 1.0_dp - pref(1:nz-1)/dpref(2:nz)*lnppref(2:nz)
    Trefv(:) = Tref

! set g matrix of the divergence equation
    gmatrix(:,:) = 0.0_dp
    do j=1, nz
      do k=1, j
        gmatrix(k,j) = rd*lnppref(j) 
      end do
      gmatrix(j,j) = rd*alpharef(j)
    end do

! set matrix T of the energy conversion term
    bmatrix(:,:) = 0.0_dp
    do j=1, nz
      bmatrix(j,j) = alpharef(j)
      do k=j+1, nz
        bmatrix(k,j) = lnppref(k)/dpref(k)*dpref(j)
      end do
    end do
! triangular matrix multiplication: gmatrix tmatrix
    call dtrmm("l","u","n","n",nz,nz,1.0_dp,gmatrix,nz,bmatrix,nz)
! add vector nu of the surface pressure tendency equation
! rank-1 update: bmatrix += Rd Tr dpref/psref
    call dger(nz,nz,1.0_dp,Trefv,1,pref/psref,1,bmatrix,nz)
    deallocate(lnppref,alpharef,gmatrix,trefv)

  end subroutine sigmap_semiimplicit

  function bfunc(eta,etakp) result(bf)
    implicit none

    real(kind=dp), intent(in) :: eta, etakp
    real(kind=dp) :: bf

    bf = (eta-etakp)/(1.0_dp-etakp)

  end function bfunc

  subroutine setp(ph)
    implicit none

    real(kind=dp), dimension(:,:,0:), intent(in) :: ph

    integer(kind=i4b) :: nx, ny, nz

    nx = size(ph,1)
    ny = size(ph,2)
    nz = size(ph,3) - 1

    if (.not.allocated(dpf)) then
      allocate(dpf(nx,ny,nz),dpfr(nx,ny,nz),lnpp(nx,ny,nz),alpha(nx,ny,nz))
    end if
    dpf(:,:,:) = ph(:,:,1:nz)-ph(:,:,0:nz-1)
    dpfr(:,:,:) = 1.0_dp/dpf(:,:,:)
    lnpp(:,:,2:nz) = log(ph(:,:,2:nz)/ph(:,:,1:nz-1))
    lnpp(:,:,1) = 2.0_dp
    alpha(:,:,2:nz) = 1.0_dp-ph(:,:,1:nz-1)*dpfr(:,:,2:nz)*lnpp(:,:,2:nz)
    alpha(:,:,1) = log(2.0_dp)

  end subroutine setp

end module sigmap_module
