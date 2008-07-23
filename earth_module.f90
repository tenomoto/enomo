module earth_module
  use type_module, only: i4b, dp
  use math_module, only: pi=>math_pi
  implicit none

  real(kind=dp), parameter :: &
    earth_radius = 6.371e6_dp, &
    earth_daysec = 86400.0_dp, &
    earth_omega  = 2.0_dp*pi/earth_daysec, &
    earth_gravity = 9.80665_dp

end module earth_module
