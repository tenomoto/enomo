module machine_module

! Source: Based on zeps() in gaqd.f
! Author: T. Enomoto

  use kind_module, only : dp
  private

  public :: machine_eps

contains

  function machine_eps() result (eps)
! zeps of gaqd.f
    implicit none

    real(kind=dp) :: eps, a, b, c

    a = 4.0_dp/3.0_dp
    b = a - 1.0_dp
    c = b + b + b
    eps = abs(c-1.0_dp)

  end function machine_eps

end module machine_module
