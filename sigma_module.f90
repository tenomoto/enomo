module sigma_module
! calculates p velocity (omega) from u, v, ps and sigma levels
  use kind_module, only: i4b, dp
  implicit none
  private

  logical, public :: sigma_verbose = .false.

  real(kind=dp), dimension(:), allocatable, public :: &
    sigma_half, sigma_full, sigma_dsigma, sigma_a, sigma_b
  real(kind=dp), dimension(:), allocatable, private :: dsr
  integer(kind=i4b) :: nz

  public :: sigma_init, sigma_half2full, sigma_dot, sigma_omega, sigma_omega_afes, sigma_clean

contains

  subroutine sigma_init(h)
    implicit none

! Sigma levels are assumed to be increasing from k=1 to k=nz (from the top to the bottom)
    real(kind=dp), dimension(:), intent(in) :: h

    integer(kind=i4b) :: k

    nz = size(h)-1
    allocate(sigma_half(nz+1),sigma_full(nz),sigma_dsigma(nz), &
      sigma_a(nz),sigma_b(nz),dsr(nz))

    sigma_half(:) = h(:)

    call sigma_half2full(sigma_half,sigma_full)
    do k=1, nz
      sigma_dsigma(k) = h(k+1) - h(k)
    end do
    dsr(:) = 1.0_dp/sigma_dsigma(:)

    call sigma_ab(sigma_half,sigma_full,sigma_a,sigma_b)

    if (sigma_verbose) then
      print *, "k sigma_half sigma_full sigma_dsigma sigma_a sigma_b"
      do k=1, nz
        print *,  k, sigma_half(k), sigma_full(k), sigma_dsigma(k), sigma_a(k), sigma_b(k)
      end do
      print *, nz+1, sigma_half(k)
    end if

  end subroutine sigma_init

  subroutine sigma_clean
    implicit none

    deallocate(sigma_half,sigma_full,sigma_dsigma,sigma_a,sigma_b,dsr)

  end subroutine sigma_clean

  subroutine sigma_dot(u,v,d,ps,dlnpsx,dlnpsy,sdot)
! diagnose dsigma/dt NB. sdot is defined at half levels
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: u, v, d
    real(kind=dp), dimension(:,:), intent(in) :: ps, dlnpsx, dlnpsy
    real(kind=dp), dimension(:,:,0:), intent(inout) :: sdot

    real(kind=dp), dimension(nz) :: dsum
    integer(kind=i4b) :: nx, ny, i, j, k

    nx = size(u,1)
    ny = size(u,2)
    do j=1, ny
      do i=1, nx 
        dsum(:) = 0.0_dp ! sdot(0) = sdot(nz) = 0
        do k=1, nz
          dsum(k) = dsum(k) + (d(i,j,k) +  &
            u(i,j,k)*dlnpsx(i,j) + v(i,j,k)*dlnpsy(i,j))*sigma_dsigma(k)
        end do
        do k=1, nz-1
          sdot(i,j,k) = sigma_half(k)*dsum(nz) - dsum(k)
        end do
      end do
    end do 

  end subroutine sigma_dot

  subroutine sigma_omega(u,v,d,ps,dlnpsx,dlnpsy,omega)
! diagnose omega using energy conserving vertical integral (AFES manual)
    use air_module, only: r=>air_rd, cp=>air_cp
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: u, v, d
    real(kind=dp), dimension(:,:), intent(in) :: ps, dlnpsx, dlnpsy
    real(kind=dp), dimension(:,:,:), intent(inout) :: omega

    real(kind=dp), parameter :: cpr = cp/r
    real(kind=dp) :: dsum, dsum1, udlnps
    integer(kind=i4b) :: nx, ny, i, j, k

    nx = size(u,1)
    ny = size(u,2)
    do j=1, ny
      do i=1, nx 
        dsum = 0.0_dp
        dsum1 = 0.0_dp
        do k=1, nz
          udlnps = u(i,j,k)*dlnpsx(i,j) + v(i,j,k)*dlnpsy(i,j)
          dsum = dsum + (d(i,j,k) + udlnps)*sigma_dsigma(k)
          omega(i,j,k) = sigma_full(k)*ps(i,j) * ( &
            udlnps - cpr*dsr(k)*(sigma_a(k)*dsum + sigma_b(k)*dsum1) )
          dsum1 = dsum
        end do
      end do
    end do 

  end subroutine sigma_omega

  subroutine sigma_omega_afes(u,v,d,ps,dlnpsx,dlnpsy,omega)
! diagnose omega as coded in AFES
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: u, v, d
    real(kind=dp), dimension(:,:), intent(in) :: ps, dlnpsx, dlnpsy
    real(kind=dp), dimension(:,:,:), intent(inout) :: omega

    real(kind=dp) :: vivps, vidiv, udlnps
    integer(kind=i4b) :: nx, ny, i, j, k

    nx = size(u,1)
    ny = size(u,2)
    do j=1, ny
      do i=1, nx 
        udlnps = u(i,j,1)*dlnpsx(i,j) + v(i,j,1)*dlnpsy(i,j)
        vidiv = d(i,j,1)*sigma_dsigma(1)
        vivps = udlnps*sigma_dsigma(1)
        omega(i,j,1) = ps(i,j) * (sigma_full(1)*udlnps &
                     -(d(i,j,1)*(sigma_full(1)-sigma_half(1))) )
        do k=2, nz
          udlnps = u(i,j,k)*dlnpsx(i,j) + v(i,j,k)*dlnpsy(i,j)
          omega(i,j,k) = ps(i,j) * (sigma_full(k)*udlnps - &
                        (vivps+udlnps*(sigma_full(k)-sigma_half(k))) &
                       -(vidiv+d(i,j,k)*(sigma_full(k)-sigma_half(k))) )
          vidiv = vidiv + d(i,j,k)*sigma_dsigma(k)
          vivps = vivps + udlnps*sigma_dsigma(k)
        end do
      end do
    end do 

  end subroutine sigma_omega_afes

  subroutine sigma_half2full(half,full)
! Phillips (1974)
    use air_module, only: r=>air_rd, cp=>air_cp

    real(kind=dp), dimension(:), intent(in) :: half 
    real(kind=dp), dimension(:), intent(inout) :: full

    integer(kind=i4b) :: nz, k
    real(kind=dp) :: kr, k1

    nz = size(full)

    kr = cp/r
    k1= (cp+r)/cp ! kappa+1
    do k=1, nz
      full(k) = ((half(k+1)**k1-half(k)**k1)/(k1*(half(k+1)-half(k))))**kr
    end do

  end subroutine sigma_half2full

  subroutine sigma_ab(half,full,a,b)
    use air_module, only: kappa=>air_kappa
    implicit none

    real(kind=dp), dimension(:), intent(in) :: half, full
    real(kind=dp), dimension(:), intent(inout) :: a, b

    integer(kind=i4b) :: k, nz
    real(kind=dp) :: hf

    nz = size(full)

    do k=1, nz
      a(k) = (half(k+1)/full(k))**kappa - 1.0_dp
    end do
    b(1) = 0.0_dp
    do k=2, nz
      b(k) = 1.0_dp - (half(k)/full(k))**kappa
    end do

  end subroutine

end module sigma_module
