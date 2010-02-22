module svd_module
  use kind_module, only: i4b, dp
  implicit none

  integer(kind=i4b), private :: m, n, lda, ldu, ldvt, lwork
  real(kind=dp), dimension(:), allocatable, private :: work

  logical, public :: svd_verbose = .false.

  public :: svd_init, svd_calc, svd_clean

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

  subroutine svd_calc(a,u,sigma,vt)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), intent(inout) :: u, vt
    real(kind=dp), dimension(:), intent(inout) :: sigma

    character, parameter :: jobu = "S", jobvt = "S"
  
    integer(kind=i4b) :: info

    call dgesvd(jobu, jobvt, size(a,1), size(a,2), a, lda, &
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

  subroutine svd_clean()
    implicit none

    deallocate(work)

  end subroutine svd_clean

end module svd_module
