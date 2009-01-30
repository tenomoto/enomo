module sphere_module
! utility for a spherical topology
	use type_module, only: i4b, dp
	use math_module, only: pi=>math_pi, pi2=>math_pi2, pih=>math_pih, pir=>math_pir
	private

	public :: sphere_xy2lon, sphere_lonlat2xyz, sphere_uv2xyz, sphere_xyz2uv,     &
            sphere_j2j, sphere_ij2i, &
            sphere_orthodrome, sphere_cosine, sphere_trcosine, sphere_sine,     &
            sphere_lat2y, sphere_y2lat

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

  function sphere_j2j(j,ny) result(jout)
  ! returns 1<=jout<=ny from j that could be beyond poles (j<1, j>ny).
    implicit none

    integer(kind=i4b), intent(in) :: j, ny
    integer(kind=i4b) :: jout

    if (j<1) then
      jout = 1 - j
    else if (j>ny) then
      jout = ny + (ny-j) + 1
    else
      jout = j
    end if

  end function sphere_j2j

  function sphere_ij2i(i, j, nx, ny) result(iout)
  ! returns 1<=iout<=ny from i that could be i>nx and 
  !                          j that could be beyond poles (j<1, j>ny).
    implicit none

    integer(kind=i4b), intent(in) :: i, j, nx, ny
    integer(kind=i4b) :: iout

    if ((j<1).or.(j>ny)) then
      iout = i+sign(nx/2,nx/2-i) ! i+nx/2 (i<=nx/2), i-nx/2 (i>nx/2)
    else
      iout = mod(i-1,nx)+1
    end if
      
  end function sphere_ij2i

	function sphere_orthodrome(lon1, lat1, lon2, lat2) result(l)
	! returns the shortest distance between two points on the unit sphere
		implicit none

		real(kind=dp), intent(in) :: lon1, lat1, lon2, lat2
		real(kind=dp) :: l

    l = acos(sphere_cosine(sin(lat1),sin(lat2),cos(lat1),cos(lat2),lon1-lon2))

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

! mercator projection
   function sphere_lat2y(lat) result(y)
     implicit none

     real(kind=dp), intent(in) :: lat
     real(kind=dp) :: y

     y = log(tan(0.5_dp*(pih+lat)))

   end function sphere_lat2y

   function sphere_y2lat(y) result(lat)
     implicit none

     real(kind=dp), intent(in) :: y
     real(kind=dp) :: lat

     lat = 2.0_dp*atan(exp(y)) - pih

   end function sphere_y2lat

end module sphere_module

