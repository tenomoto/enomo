module integer_module
  use kind_module, only: i4b
  implicit none
  private

  public :: integer_swap

contains

  subroutine integer_swap(a,b)
    integer(kind=i4b), intent(inout) :: a, b

    a = a + b
    b = a - b
    a = a - b

  end subroutine integer_swap

end module integer_module
