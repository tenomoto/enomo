module matrix_module
  use kind_module, only: i4b, dp
  implicit none
  private

  public :: matrix_multiply

contains

  function matrix_multiply(a,b) result(c)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), intent(in) :: b

    real(kind=dp), dimension(size(a,1),size(b,2)) :: c
    real(kind=dp), parameter :: alpha = 1.0d0, beta = 0.d0

    integer :: m, n, k, ka, kb, lda, ldb, ldc
    character(len=1), parameter :: transa = "N", transb = "N"

    m = size(a,1)
    ka = size(a,2)
    kb = size(b,1)
    n = size(b,2)
    if (ka/=kb) then
      print *, "error in matrix matrix multipy: size(a,2)=", &
        ka, " is not equal to size(b,1)=", kb
      stop
    end if
    k = ka
    lda = m
    ldb = ka
    ldc = m

    call dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

  end function matrix_multiply

end module matrix_module
