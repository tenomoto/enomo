module machine_module

! Source: Based on zeps() in gaqd.f
! Author: T. Enomoto

  use kind_module, only : dp, qp
  private

  public :: machine_eps, machine_epsq

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

  function machine_epsq() result (eps)
! zeps of gaqd.f
    implicit none

    real(kind=qp) :: eps, a, b, c

    a = 4.0_qp/3.0_qp
    b = a - 1.0_qp
    c = b + b + b
    eps = abs(c-1.0_qp)

  end function machine_epsq

end module machine_module
