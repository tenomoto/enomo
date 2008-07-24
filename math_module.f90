module math_module
  use type_module, only: dp
  implicit none

  real(kind=dp), parameter, public :: &
    math_pi  = 3.141592653589793, &
    math_pi2 = 2.0_dp*math_pi, math_pih = 0.5_dp*math_pi, math_pir = 1.0_dp/math_pi, &
    math_deg2rad = math_pi/180.0_dp, math_rad2deg = 180.0_dp*math_pir

end module math_module
