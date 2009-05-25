module search_module
  use kind_module, only: i4b, dp
  implicit none
  private

  public :: search_linear, search_bisection

contains
  function search_linear(y, x, i00) result(i)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: y
    real(kind=dp) :: x
    integer(kind=i4b), intent(in), optional :: i00

    integer(kind=i4b) :: i

    integer(kind=i4b) :: n, i0
    real(kind=dp) :: s

    n = size(y)
! ascending y: s =  1
! decending y: s = -1
    s = sign(1.0_dp, y(n)-y(1))

    if (present(i00)) then
      i0 = i00
    else
      i0 = 1
    end if
! out of range
    if (s*(x-y(1))<0.0_dp) then
      i = 0
      return
    else if (s*(x-y(n))>=0.0_dp) then
      i = n
      return
    end if
    i0 = max(1,min(n-1,i0))

    if (s*(x-y(i0+1))<0.0_dp) then
      if (s*(x-y(i0))>=0.0_dp) then
! ascending y: y(i0) <= x < y(i0+1)
! descending y: y(i0) > x >= y(i0+1)
        i = i0
        return
      else
! search with decreasing i
        do i=i0, 1, -1
          if (s*(x-y(i))>=0.0_dp) then
            exit
          end if
        end do
      end if
    else
! search with increasing i
      do i=i0+1, n-1
        if (s*(x-y(i+1))<0.0_dp) then
          exit
        end if
      end do
    end if

  end function search_linear

  function search_bisection(y, x) result(i)
    implicit none

    real(kind=dp), dimension(:), intent(in) :: y
    real(kind=dp) :: x

    integer(kind=i4b) :: i

    integer(kind=i4b) :: n, il, im, iu
    logical :: lascend

    n = size(y)
    lascend = y(n) > y(1)
    il = 0
    iu = n+1
    do
      if (iu-il==1) then
        exit
      end if
      im = (iu+il)/2
      if ((x > y(im)).eqv.lascend) then
        il = im
      else
        iu = im
      end if
    end do
    i = il

  end function search_bisection

end module search_module
