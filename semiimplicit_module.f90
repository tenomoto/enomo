module semiimplicit_module
  use kind_module, only: i4b, dp
  use earth_module, only: earth_radius
  implicit none
  private

  public :: semiimplicit_inverse

contains

  subroutine semiimplicit_inverse(bmatrix, beta, dt, invamatrix)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: bmatrix
    real(kind=dp), intent(in) :: beta, dt
    real(kind=dp), dimension(:,:,0:), intent(inout) :: invamatrix

    integer(kind=i4b) :: nz, ntrunc, n, k
    real(kind=dp) :: f

! lapack related variables
    integer(kind=i4b) :: info
    integer(kind=i4b), dimension(size(bmatrix,1)) :: ipiv
    real(kind=dp), dimension(1) :: work_dummy
    real(kind=dp), dimension(:), allocatable :: work

    nz = size(bmatrix,1)
    ntrunc = size(invamatrix,3)-1

    invamatrix(:,:,:) = 0.0_dp
    do k=1, nz
      invamatrix(k,k,:) = 1.0_dp
    end do

    call dgetri(nz,invamatrix,nz,ipiv,work_dummy,-1,info)
    allocate(work(int(work_dummy(1))))

    f = (1.0_dp/earth_radius*beta*dt)**2 ! 1/a^2 (beta dt)^2

    do n=0, ntrunc
      invamatrix(:,:,n) = invamatrix(:,:,n) + n*(n+1)*f*bmatrix(:,:)
      call dgetri(nz,invamatrix,nz,ipiv,work,size(work),info)
    end do

  end subroutine semiimplicit_inverse

end module semiimplicit_module
