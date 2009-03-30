module water_module
  use kind_module, only: dp
  use air_module, only: air_rd
  implicit none

  real(kind=dp), parameter, public ::  &
    water_t0  = 273.15_dp, & ! freezing temperature
    water_rv  = 461.0_dp, & ! gas constant for water vapour, J/deg/kg
    water_lv0 = 2.5e6_dp, & ! latent heat of vapouriztion at 0C , J/kg
    water_eps = air_rd/water_rv

end module water_module
