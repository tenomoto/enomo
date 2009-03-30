module upstream_module

! finds departure and mid-points

! Reference: Ritchie (1987), Ritchie and Beaudoin (1994)
! Method:
!   use the Cartesian coordinates with the origin at the centre of the sphere (Rithchie 1987)
!   approximations far from the poles (Ritchie and Beaudoin 1994)

! Author: T. Enomoto
! History: 
!   2009-01-27 Modified for enomo
!   2004-02-26 First version

  use kind_module, only: i4b, dp
  use math_module, only: pi2=>math_pi2
  use earth_module, only: a=>earth_radius
  private

  integer(kind=i4b), private :: kmax = 2, jc, nx, ny, n = 2, m = 2
  real(kind=dp), private :: &
    small = 1.0_dp-10, &
    latc = 0.0_dp ! use find_lonlat at |lat| < latc

  real(kind=dp), dimension(:), allocatable, private :: seclat, seclat2, tanlat
  real(kind=dp), dimension(:,:), allocatable, private :: ubuf, vbuf, buf

  public :: upstream_init, upstream_clean, upstream_find, upstream_calcxy
  private :: find_xyz, find_lonlat

contains

  subroutine upstream_init(lon, lat, nn, mm, itermax, delta, latpol)
    use glatwgt_module, only: glatwgt_within
    use regrid_module, only: regrid_init
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    integer(kind=i4b), optional :: nn, mm, itermax
    real(kind=dp), optional :: delta, latpol

    real(kind=dp) :: dlatr
    integer(kind=i4b) :: j

    nx = size(lon)
    ny = size(lat)
    if (present(nn)) then
      n = nn
    end if
    if (present(mm)) then
      m = mm
    end if
    allocate(ubuf(nx+2*n,ny+2*m),vbuf(nx+2*n,ny+2*m),buf(nx+2*n,ny+2*m))
    call regrid_init(lon,lat,n,m)
    if (present(itermax)) then
      kmax = itermax
    end if
    if (present(delta)) then
      small = delta
    end if
    if (present(latpol)) then
      latc = abs(latpol)
    end if
    latc = sign(latc,lat(1))
    if ((lat(2)-lat(1))==(lat(3)-lat(2))) then ! equispaced
      dlatr = 1.0_dp/(lat(2)-lat(1))
      jc = floor((latc-lat(1))*dlatr+1)
    else ! assume Gaussian latitudes
      jc = glatwgt_within(lat,latc)
    end if
    jc = max(1,jc)
    if (jc<ny/2) then ! use find_lonlat
      print *, "jc =", jc, " approximation between", lat(jc+1), " and", lat(ny-jc)
      allocate(seclat(ny), seclat2(ny), tanlat(ny))
      do j=jc+1, ny-jc
        seclat(j) = 1.0_dp/cos(lat(j))
        seclat2(j) = seclat(j)*seclat(j)
        tanlat(j) = tan(lat(j))
      end do
    end if

  end subroutine upstream_init

  subroutine upstream_clean()
    use regrid_module, only: regrid_clean
    implicit none

    deallocate(ubuf, vbuf, buf)
    call regrid_clean()

  end subroutine upstream_clean

  subroutine upstream_find(u, v, dt, deplon, deplat, midlon, midlat, umid, vmid)
    use regrid_module, only: regrid_extend
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: u, v
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(:,:), intent(inout) :: deplon, deplat
    real(kind=dp), dimension(:,:), optional, intent(inout) :: &
      midlon, midlat, umid, vmid

    integer(kind=i4b) :: nx, ny, i, j, jr
    real(kind=dp) :: mlon, mlat, um, vm
    logical :: lmid

    nx = size(u,1)
    ny = size(u,2)

    lmid = present(midlon).and.present(midlat).and. &
           present(umid)  .and.present(vmid)
    call regrid_extend(u, ubuf, .true.)
    call regrid_extend(v, vbuf, .true.)
    ubuf = ubuf/a
    vbuf = vbuf/a
    ! use RB94 approximation far from the poles
    do j=jc+1, ny-jc
      do i=1, nx
        call find_lonlat(i,j,dt,deplon(i,j),deplat(i,j),mlon,mlat,um,vm)
        if (lmid) then
          midlon(i,j) = mlon 
          midlat(i,j) = mlat
          umid(i,j) = um*a
          vmid(i,j) = vm*a
        end if   
      end do
    end do
    do j=1, jc
      jr = ny-j+1
      do i=1, nx
        call find_xyz(i,j,dt,deplon(i,j),deplat(i,j),mlon,mlat,um,vm)
        if (lmid) then
          midlon(i,j) = mlon 
          midlat(i,j) = mlat
          umid(i,j) = um*a
          vmid(i,j) = vm*a
        end if   
        call find_xyz(i,jr,dt,deplon(i,jr),deplat(i,jr),mlon,mlat,um,vm)
        if (lmid) then
          midlon(i,jr) = mlon 
          midlat(i,jr) = mlat
          umid(i,jr) = um*a
          vmid(i,jr) = vm*a
        end if   
      end do
    end do
    ! use find_xyz() first and last jc latitudes
! debug start
!    print *, "upstream_find"
!    print *, maxval(u), minval(u)
!    print *, maxval(v), minval(v)
!    print *, maxval(deplon), "at", maxloc(deplon)
!    print *, minval(deplon), "at", minloc(deplon)
!    print *, deplat(::4,61)
! debug end

  end subroutine upstream_find

  subroutine find_xyz(i, j, dt, dlon, dlat, mlon, mlat, un, vn)
    use sphere_module, only: &
      xy2lon=>sphere_xy2lon, lonlat2xyz=>sphere_lonlat2xyz, uv2xyz=>sphere_uv2xyz
    use regrid_module, only: regrid_lon, regrid_lat, &
      regrid_bilinear, regrid_linpol, regrid_spcher
    implicit none

    integer(kind=i4b), intent(in) :: i, j
    real(kind=dp), intent(in) :: dt
    real(kind=dp), intent(out) :: dlon, dlat, mlon, mlat, un, vn

    integer(kind=i4b) :: k
    real(kind=dp) :: & 
      xd, yd, zd, & ! Cartesian velocity
      xg, yg, zg, & ! arrival point in Cartesian coordinates
      x0, y0, z0, & ! present point in Cartesian coordinates
      x1, y1, z1, & ! updated point in Cartesian coordinates
      b,          & ! correction factor
      alon, alat

    alon = regrid_lon(i+n)
    alat = regrid_lat(j+m)
    ! backward initial guess (use wind at arrival point)
    un = ubuf(i+n,j+m)
    vn = vbuf(i+n,j+m)
    mlon = alon
    mlat = alat
    ! transform into Cartesian coordinates
    call lonlat2xyz(alon, alat, xg, yg, zg) 
    do k=1, kmax
      if (k>1) then
        un = regrid_linpol(ubuf, mlon, mlat)
        vn = regrid_linpol(vbuf, mlon, mlat)
      end if
      ! normalised Cartesian velocity
      call uv2xyz(un,vn,mlon,mlat,xd,yd,zd)
      ! correction factor
      b = 1.0_dp/sqrt(1.0_dp+dt*dt*(xd*xd+yd*yd+zd*zd)-2.0_dp*dt*(xd*xg+yd*yg+zd*zg))
      ! calculate new points
      x1 =  b*(xg - dt*xd)
      y1 =  b*(yg - dt*yd)
      z1 =  b*(zg - dt*zd)
      ! calculate (lon,lat) from (x,y,z)
      mlat = asin(z1)
      mlon = xy2lon(x1,y1)
    end do
    b = 2.0_dp*(x1*xg+y1*yg+z1*zg) ! calculate the departure point
    x1 = b*x1 - xg
    y1 = b*y1 - yg
    z1 = b*z1 - zg
    dlon = xy2lon(x1,y1)
    dlat = asin(z1)
    
  end subroutine find_xyz

  subroutine find_lonlat(i,j,dt,dlon,dlat,mlon,mlat,un,vn)
    use regrid_module, only: regrid_lon, regrid_lat, &
      regrid_bilinear, regrid_linpol, regrid_spcher
    implicit none

    integer(kind=i4b), intent(in) :: i, j
    real(kind=dp), intent(in) :: dt
    real(kind=dp), intent(out) :: dlon, dlat, mlon, mlat, un, vn

    real(kind=dp), parameter :: f6 = 1.0_dp/6.0_dp, f23 = 2.0_dp/3.0_dp

    real(kind=dp) :: alon, alat, undt, vndt, undt2, vmagdt2
    integer(kind=i4b) :: k

    alon = regrid_lon(i+n)
    alat = regrid_lat(j+m)
    ! backward initial guess (use wind at arrival point)
    un = ubuf(i+n,j+m)
    vn = vbuf(i+n,j+m)
    do k=1, kmax
      if (k>1) then
        un = regrid_linpol(ubuf, mlon, mlat)
        vn = regrid_linpol(vbuf, mlon, mlat)
      end if
      undt = un*dt
      undt2 = undt*undt
      vndt = vn*dt
      vmagdt2 = (un*un+vn*vn)*dt
      mlon = alon - undt*seclat(j)*(1.0_dp+f6*(undt2*seclat2(j)-vmagdt2))
      if (mlon<0.0_dp) then
        mlon = mlon + pi2
      else if (mlon>=pi2) then
        mlon = mlon - pi2
      end if
      mlat = alat - vndt+0.5_dp*tanlat(j)*undt2
    end do
    dlon = alon - 2.0_dp*seclat(j)*undt*(1.0_dp-tanlat(j)*vndt)
    if (dlon<0.0_dp) then
      dlon = dlon + pi2
    else if (dlon>=pi2) then
      dlon = dlon - pi2
    end if
    dlat = alat - 2.0_dp*vndt+(seclat2(j)-f23)*undt2*vndt

  end subroutine find_lonlat

  subroutine upstream_calcxy(u, v, mlat, dta, xd, yd, xa, ya, lapprox)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: u, v, mlat ! u,v at midpoint
    real(kind=dp), intent(in) :: dta ! dta = dt/a
    real(kind=dp), dimension(:,:), intent(inout) :: xd, yd, xa, ya
    logical, intent(in), optional :: lapprox ! use RB94 approximation

    real(kind=dp) :: vm, sml, cml, al, sal, cal, cal1, sgm, cgm, x1, x2, y1, y2
    integer(kind=i4b) :: i, j, nx, ny

    nx = size(u,1)
    ny = size(u,2)
    do j=1, ny
      do i=1, nx
        vm = sqrt(u(i,j)**2+v(i,j)**2)
        if (vm.ne.0) then
           sml = sin(mlat(i,j))
           cml = cos(mlat(i,j))
           if (present(lapprox).and.lapprox) then
!   approximation by Ritchie and Beaudoin (1994)
             x2 = dta*sml
             x1 = u(i,j)*x2
             y1 = v(i,j)*x2
             xd(i,j) = -x1
             yd(i,j) = cml + y1
             xa(i,j) = x1
             ya(i,j) = cml - y1
           else
!   Ritchie (1988) formulation
            al = vm*dta
            sal = sin(al)
            cal = cos(al)
            sgm = v(i,j)/vm
            cgm = u(i,j)/vm
            cal1 = 1.0d0-cal
            x1 = sal*cgm*sml
            x2 = cal1*sgm*cgm*cml
            y1 = sal*sgm*sml
            y2 = cal1*sgm*sgm*cml
            xd(i,j) = -x1 + x2
            yd(i,j) = cml + y1 - y2
            xa(i,j) = x1 + x2
            ya(i,j) = cml - y1 - y2
          end if
        else
          xd(i,j) = 0.0_dp
          yd(i,j) = cml
          xa(i,j) = 0.0_dp
          ya(i,j) = cml
        end if
      end do
    end do

  end subroutine upstream_calcxy

end module upstream_module

