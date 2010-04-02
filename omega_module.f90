module omega_module
  use kind_module, only: i4b, dp
  implicit none
  private

  real(kind=dp), dimension(:), allocatable :: hsig, fsig, dsig, ak, bk
  integer(kind=i4b) :: nz

  public :: omega_init, omega_calc, omega_clean

contains

  subroutine omega_init(h, f)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: h, f

    integer(kind=i4b) :: k

    nz = size(f)
    if (size(h)<=nz) then
      print *, "Error in omega_init: size(h) <= size(f)"
      stop
    end if
     
    allocate(hsig(nz+1),fsig(nz),dsig(nz))
    hsig(:) = h(:)
    fsig(:) = f(:)
    do k=1, nz
      dsig(:) = h(k+1) - h(k)
    end do

    ak(nz) = -log(hsig(nz))
    do k=nz-1, 1, -1
      ak(k) = 0.5_dp*log(hsig(k+1)/hsig(k))
    end do
    bk(1) = 0.0_dp
    do k=2, nz
      bk(k) = ak(k)+ak(k-1)
    end do

    if (any(dsig(:)<=0)) then
      print *, "Error in omega_init: negative layer thinckess"
      stop
    end if

  end subroutine omega_init

  subroutine omega_calc(u,v,d,ps,dlnpsx,dlnpsy,omega)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: u, v, d
    real(kind=dp), dimension(:,:), intent(in) :: ps, dlnpsx, dlnpsy
    real(kind=dp), dimension(:,:,:), intent(out) :: omega

    real(kind=dp) :: gk, gksum, udlnps
    integer(kind=i4b) :: nx, ny, i, j, k

    nx = size(u,1)
    ny = size(u,2)
    do j=1, ny
      do i=1, nx 
        gksum = 0.0_dp
        do k=1, nz
          udlnps = u(i,j,k)*dlnpsx(i,j) + v(i,j,k)*dlnpsy(i,j)
          gk = d(i,j,k) + udlnps
          omega(i,j,k) = fsig(k)*ps(i,j) * ( &
            udlnps - ak(k)*gk - bk(k)*gksum )
          gksum = gksum + gk*dsig(k)
        end do
      end do
    end do 

  end subroutine omega_calc

  subroutine omega_clean
    implicit none

    deallocate(hsig,fsig,dsig,ak,bk)

  end subroutine omega_clean

end module omega_module
