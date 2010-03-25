module svd_module
  use kind_module, only: i4b, dp
  implicit none

  integer(kind=i4b), private :: m, n, lda, ldu, ldvt, lwork
  real(kind=dp), dimension(:), allocatable, private :: work

  logical, public :: svd_verbose = .false.

  public :: svd_init, svd_clean, svd_calc, svd_inverse

contains

  subroutine svd_init(a)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: a

    m = size(a, 1)
    n = size(a, 2)
    lda = max(1,m)
    ldu = lda
    ldvt = min(m,n)
    lwork = max(3*ldvt+max(m,n),5*ldvt)
    if (svd_verbose) then
      print *, "m=", m, " n=", n
      print *, "lda=", lda, " ldu=", ldu, " ldvt=", ldvt, " lwork=", lwork
    end if

    allocate(work(lwork))

  end subroutine svd_init

  subroutine svd_calc(x,u,sigma,vt)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: x
    real(kind=dp), dimension(:,:), intent(inout) :: u, vt
    real(kind=dp), dimension(:), intent(inout) :: sigma

    character, parameter :: jobu = "S", jobvt = "S"
  
    integer(kind=i4b) :: info

    call dgesvd(jobu, jobvt, size(x,1), size(x,2), x, lda, &
      sigma, u, ldu, vt, ldvt, work, lwork, info)
    if (svd_verbose) then
      print *, "dgesvd info=", info
      print *, "u max", maxval(u), " min", maxval(u)
      print *, "sigma max", maxval(sigma), " min", maxval(sigma)
      print *, "vt max", maxval(vt), " min", maxval(vt)
    end if

    if (info.ne.0) then
      print *, "### Error in sgesvd."
      stop
    end if

  end subroutine svd_calc

  subroutine svd_inverse(x,v,sigma,ainv)
    use matrix_module, only: matrix_multiply
! calculates inverse of A = XX'
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: x, v
    real(kind=dp), dimension(:), intent(in) :: sigma
    real(kind=dp), dimension(:,:), intent(inout) :: ainv

    real(kind=dp), dimension(size(x,1),size(x,2)) :: xv
    integer(kind=i4b) :: j

    xv = matrix_multiply(x,v)
    do j=1, size(sigma)
      if (sigma(j)/=0) then
        xv(:,j) = xv(:,j)/(sigma(j)**2)
      else
        xv(:,j) = 0
      end if
    end do
    ainv = matrix_multiply(xv,transpose(xv))

  end subroutine svd_inverse

  subroutine svd_clean()
    implicit none

    deallocate(work)

  end subroutine svd_clean

end module svd_module
