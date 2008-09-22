module confmap_module
! generates grid and map factors using conformal mapping
  use type_module, only: i4b, dp, dpc
  implicit none
  private

! References
! Bentsen et al. (1999), Mon. Wea. Rev.
! Murray (1996), J. Comp. Phys.
! Couriter and Geleyn (1988) QJRMS

! History
! 2008-9-18 Takeshi Enomoto <eno@jamstec.go.jp>
! * first version

  public :: confmap_stereo, confmap_invstereo, confmap_linear, &
            confmap_grid2earth
contains

  function confmap_stereo(theta, phi) result(z)
    use math_module, only: math_i
    implicit none

    real(kind=dp), intent(in) :: theta, phi
    complex(kind=dpc) :: z

! 0<=phi<pi
    z = tan(0.5_dp*phi)*exp(math_i*theta)
     
  end function confmap_stereo

  subroutine confmap_invstereo(z,theta,phi)
	  use math_module, only: math_arg
    implicit none

    complex(kind=dpc), intent(in) :: z
    real(kind=dp), intent(out) :: theta, phi

    real(kind=dp) :: x, y

    theta = math_arg(z)
    phi = 2.0_dp*atan(abs(z))

  end subroutine confmap_invstereo 

  function confmap_linear(z,a,b,c,d) result (w)
    implicit none

    complex(kind=dpc), intent(in) :: z, a, b, c, d
    complex(kind=dpc) :: w

    w = (a*z+b)/(c*z+d)
 
  end function confmap_linear

  subroutine confmap_grid2earth(lon, lat, a, b, c, lone, late)
    use math_module, only: pih=>math_pih, pi2=>math_pi2
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    complex(kind=dp), intent(in) :: a, b, c

    integer(kind=i4b) :: nx, ny, i, j
    real(kind=dp) :: colat
    real(kind=dp), dimension(size(lon),size(lat)), intent(inout) :: lone, late
    complex(kind=dp) :: w, z, ai, bi, ci, di

    nx = size(lon)
    ny = size(lat)

    ai = -b*(c-a)
    bi =  a*(c-b)
    ci =   -(c-a)
    di =     c-b
    do j=1, ny
      do i=1, nx
        w = confmap_stereo(lon(i),pih-lat(j))
        z = confmap_linear(z,ai,bi,ci,di)
        call confmap_invstereo(z,lone(i,j),colat) 
        if (lone(i,j)<0) then
          lone(i,j) = lone(i,j) + pi2
        end if
        late(i,j) = pih-colat
      end do
    end do

  end subroutine confmap_grid2earth

end module confmap_module
