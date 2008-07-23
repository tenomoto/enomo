module sphere_module
! utility for a spherical topology
	use type_module, only: i4b, dp
	use math_module, only: pi=>math_pi, pi2=>math_pi2, pih=>math_pih, pir=>math_pir
	private

	public :: sphere_xy2lon, sphere_lonlat2xyz, sphere_uv2xyz, sphere_xyz2uv, &
            sphere_lon2i, sphere_lat2j, sphere_lat2jg, &
            sphere_orthodrome, sphere_cosine, sphere_trcosine, sphere_sine 

contains

	subroutine sphere_lonlat2xyz(lon, lat, x, y, z)
		implicit none

		real(kind=dp), intent(in) :: lon, lat
		real(kind=dp), intent(out) :: x, y, z

    x = cos(lon)*cos(lat)
		y = sin(lon)*cos(lat)
		z = sin(lat)

	end subroutine sphere_lonlat2xyz

	subroutine sphere_uv2xyz(u, v, lon, lat, xd, yd, zd)
		implicit none

		real(kind=dp), intent(in) :: u, v, lon, lat
		real(kind=dp), intent(out) :: xd, yd, zd

		xd = -u*sin(lon) - v*cos(lon)*sin(lat)
		yd =  u*cos(lon) - v*sin(lon)*sin(lat)
		zd =  v*cos(lat)

	end subroutine sphere_uv2xyz

	subroutine sphere_xyz2uv(xd, yd, zd, lon, lat, u, v)
		implicit none

		real(kind=dp), intent(in) :: xd, yd, zd, lon, lat
		real(kind=dp), intent(out) :: u, v

		u = cos(lon)*yd - sin(lon)*xd
		if (abs(lat) == pih) then
			v = -cos(lon)*xd-sin(lon)*yd ! omitted division by sin(pi/2) = 1
		else
			v = zd/cos(lat)
		end if

	end subroutine sphere_xyz2uv

	function sphere_xy2lon(x,y) result(lon)
		implicit none

		real(kind=dp), intent(in) :: x, y
		real(kind=dp) :: lon

		if (x==0.0_dp) then
			if (y>=0.0_dp) then
				lon = pih
			else
				lon = pi+pih
			end if
		else
			lon = atan(y/x)
			if (x<0.0_dp) then
				lon = lon + pi
			else if (y<0.0_dp) then ! x1 > 0.0
				lon = lon + pi2
			end if
		end if

	end function sphere_xy2lon

	function sphere_lon2i(lon,nx) result(i)
	! returns the closest longitudinal point i, not exceeding lon
  ! 1 <= return value <= nx
		implicit none

		real(kind=dp), intent(in) :: lon ! radians
		integer(kind=i4b), intent(in) :: nx

		integer(kind=i4b) :: i
		real(kind=dp) :: dlonr

		dlonr = 0.5_dp*nx*pir
		i = floor(lon*dlonr+1.0_dp) ! lon = 2pi/nx*(i-1)=dlon*(i-1)

	end function sphere_lon2i

  function sphere_lat2j(lat,ny) result(j)
	! returns the closest latitudinal point j, not exceeding lat
  ! 1 <= return value <= ny
    implicit none

    real(kind=dp), intent(in) :: lat ! radians
    integer(kind=i4b), intent(in) :: ny

    integer(kind=i4b) :: j
    real(kind=dp) :: dlatr

    dlatr = (ny-1.0_dp)*pir
    j = floor((lat+pih)*dlatr+1.0_dp) ! lat = -pi/2 + pi/(J-1)*(j-1)

  end function sphere_lat2j

	function sphere_lat2jg(lat,ny) result(j)
	! returns the closest Gaussian point j using approximation
	! lat varies from SP to NP
  ! 0 <= return value <= ny
		implicit none

		real(kind=dp), intent(in) :: lat ! radians
		integer(kind=i4b) :: ny

		integer(kind=i4b) :: j

    j = floor(0.5_dp*(ny+1.0_dp+(2.0_dp*ny+1.0_dp)*lat*pir)) ! lat = (-J-1+2j)pi/(2J+1)

	end function sphere_lat2jg

	function sphere_orthodrome(lon1, lat1, lon2, lat2) result(l)
	! returns the shortest distance between two points on the unit sphere
		implicit none

		real(kind=dp), intent(in) :: lon1, lat1, lon2, lat2
		real(kind=dp) :: l

		l = acos(cos(lon1-lon2)*cos(lat1)*cos(lat2)+sin(lat1)*sin(lat2))

	end function sphere_orthodrome

! spherical triagle composed of points pa, pb, pc and arcs a, b, c
  function sphere_cosine(cosb, cosc, sinb, sinc, pa) result(cosa)
    implicit none

    real(kind=dp), intent(in) :: cosb, cosc, sinb, sinc, pa
    real(kind=dp) :: cosa

    cosa = cosb*cosc+sinb*sinc*cos(pa)

  end function sphere_cosine

  function sphere_sine(sina, sinb, pa) result(sinpb)
    implicit none

    real(kind=dp), intent(in) :: sina, sinb, pa
    real(kind=dp) :: sinpb

    sinpb = sinb/sina*sin(pa)

  end function sphere_sine

  function sphere_trcosine(sina, sinb, sinc, cosb, cosc, pa) result(cospb)
    implicit none

    real(kind=dp), intent(in) :: sina, sinb, sinc, cosb, cosc, pa
    real(kind=dp) :: cospb

    cospb = (cosb*sinc-sinb*cosc*cos(pa))/sina

   end function sphere_trcosine

end module sphere_module
