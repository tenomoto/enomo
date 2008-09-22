module math_module
  use type_module, only: dp, dpc
  implicit none

  real(kind=dp), parameter, public :: &
    math_pi  = 3.141592653589793, &
    math_pi2 = 2.0_dp*math_pi, math_pih = 0.5_dp*math_pi, math_pir = 1.0_dp/math_pi, &
    math_deg2rad = math_pi/180.0_dp, math_rad2deg = 180.0_dp*math_pir
  complex(kind=dpc), parameter, public :: &
    math_i = (0.0_dp,1.0_dp)

  public :: math_arg

contains

  function math_arg(z) result(theta)
    implicit none

    complex(kind=dpc), intent(in) :: z
    real(kind=dp) :: theta

    real(kind=dp) :: x, y

    x = real(z,kind=dp)
    y = aimag(z)

    if ((x/=0.0_dp).and.(y/=0.0_dp)) then
      theta = atan2(y,x)
      if (theta<0) then
        theta = theta + math_pi2
      end if
    else if (x==0.0_dp) then
			if (y>0.0_dp) then
				theta = math_pih
			else if (y<0.0_dp) then
				theta = math_pi+math_pih
      else
        theta = 0.0_dp
			end if
    else if (y==0.0_dp) then
			if (x>0.0_dp) then
				theta = 0.0_dp
			else
				theta = math_pi
			end if 
    end if

  end function math_arg

end module math_module
