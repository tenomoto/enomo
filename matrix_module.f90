module matrix_module
  use kind_module, only: i4b, dp
  implicit none
  private

  public :: matrix_multiply, matrix_vector_multiply

contains

  function matrix_multiply(a,b) result(c)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), intent(in) :: b

    real(kind=dp), dimension(size(a,1),size(b,2)) :: c

    real(kind=dp), parameter :: alpha = 1.0d0, beta = 0.d0
    character(len=1), parameter :: transa = "N", transb = "N"

    integer :: m, n, k, ka, kb, lda, ldb, ldc

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

  function matrix_vector_multiply(a,x) result(y)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:), intent(in) :: x

    real(kind=dp), dimension(size(x)) :: y

    character(len=1), parameter :: trans = "N"
    real(kind=dp), parameter :: alpha = 1.0d0, beta = 0.d0
    integer(kind=i4b), parameter :: incx = 1, incy = 1

    integer :: m, n, lda

    m = size(a,1)
    n = size(a,2)
    lda = m

    call dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
   
  end function matrix_vector_multiply 

end module matrix_module
