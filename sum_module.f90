module sum_module
! calculate a more accurate sum
  use kind_module, only : i4b, dp

  public :: sum_pairwise, sum_kahan, sum_test

contains

  recursive function sum_pairwise(x) result(s)
    implicit none

! http://en.wikipedia.org/wiki/Pairwise_summation

    real(kind=dp), dimension(:), intent(in) :: x
    real(kind=dp) :: s

    integer(kind=i4b), parameter :: n = 1024
    integer(kind=i4b) :: l, i, m

    l = size(x)
    if (l<=n) then
      s = x(1)
      do i=2, l
        s = s + x(i)
      end do
    else
      m = l/2
      s = sum_pairwise(x(1:m)) + sum_pairwise(x(m+1:l))
    end if

  end function sum_pairwise

  function sum_kahan(x) result(s)
    implicit none

! http://en.wikipedia.org/wiki/Kahan_summation_algorithm

    real(kind=dp), dimension(:), intent(in) :: x

    real(kind=dp) :: s, c = 0.0_dp, y, t
    integer(kind=i4b) :: i

    s = 0.0_dp
    do i = 1, size(x)
      y = x(i) - c
      t = s + y
      c = (t - s) - y
      s = t
    end do

  end function sum_kahan

  subroutine sum_test(n)
    implicit none

    integer(kind=i4b), intent(in) :: n

    real(kind=dp), dimension(:), allocatable :: x
    real(kind=dp) :: s, sp, sk
    integer(kind=i4b) :: i

    allocate(x(n))
    call random_number(x)
    x = x * 1.0e16_dp
    s = x(1)
    do i = 2, n
      s = s + x(i)
    end do
    sp = sum_pairwise(x)
    sk = sum_kahan(x)
    print *, "obvious sum = ", s, " pairwise sum = ", sp, " kahan sum = ", sk
    print *, "pairwise - obvious = ", sp - s, "kahan - obvious = ", sk - s
    deallocate(x)
    
  end subroutine sum_test

end module sum_module
