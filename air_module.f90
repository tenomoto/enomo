module air_module
  use type_module, only : dp
  implicit none

  real(kind=dp), parameter, public :: &
    air_rd = 287.0_dp,  & ! gas constant of dry air, J/deg/kg
    air_cp = 1004.0_dp, & ! specific heat at constant pressure, J/deg/kg
    air_cv =  717.0_dp, & ! specific heat at constant volume, J/deg/kig
    air_kappa = air_rd/air_cp, &
    air_gamma = air_cp/air_cv

end module air_module
