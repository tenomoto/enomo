module vis5d_module
  use kind_module, only : i4b, sp, dp
  implicit none
  private

  integer(kind=i4b), parameter, private :: &
    compress = 1, vertical = 3, len_varname = 15
  integer(kind=i4b), private :: nr, nc, numtimes, projection
  real(kind=sp), dimension(5) :: proj_args
  integer(kind=i4b), dimension(:), alloctable :: nl
  character(len=len_varname), dimension(:), allocatable :: varname
  integer(kind=i4b), dimension(:), allocatable :: datestamp, timestamp
  real(kind=dp), dimension(:,:), allocatable :: lon, lat, t, u
  integer(kind=i4b), dimension(:,:), allocatable :: i0, i1, j0, j1

  public :: vis5d_init, vis5d_clean, vis5d_interpolate
  private :: set_grid_lonlat, set_grid_xy, create_timestamp

contains

  subroutine vis5d_init(outfname, nrow, ncol, lpolar, loc, plev, vars, nlev, &
    t0, dt, nt, lon_in, lat_in, dlat)
    use time_module, only : time_type
    implicit none

    character(len=*), intent(in) :: outfname
    integer(kind=i4b), intent(in) :: nrow, ncol
    logical, intent(in) :: lpolar                            ! T for polar stereo
    real(kind=sp), dimension(:), intent(in) :: &
      loc,            & ! xy: centre lon and lat, lon-lat: lon0, lon1, lat0, lat1 (rad)
    real(kind=sp), dimension(:), intent(in) :: &
      plev              ! pressure levels (hPa)
    character(len_varname), dimension(:), intent(in) :: vars ! names of variables
    integer(kind=i4b), dimension(:), intent(in) :: nlev      ! # of levels of each var
    type(time_type), intent(in) :: t0, dt                    ! init time and increment
    real(kind=sp), dimension(:), intent(in) :: &
      lon_in, lat_in, & ! longitudes and latitudes of input (rad)
    real(kind=sp), optional, intent(in) :: dlat ! approximate dlat for Gaussian grid

    integer(kind=i4b) :: ierr

    nr = nrow
    nc = ncol
    numtimes = nt

    if (lpolar) then
      projection = 3
      call set_grid_xy(loc(1),loc(2))
    else
      projection = 1
      call set_grid_lonlat(loc(1),loc(2),loc(3),loc(4))
    end if

    allocate(lon(nr,nc),lat(nr,nc),t(nr,nc),u(nr,nc),ii(nr,nc),jj(nr,nc))
    call set_tu(lon_in, lat_in)

    allocate(vert_args(size(plev)), varname(size(vars)), nl(size(vars)))
    vert_args(:) = plev(:)
    varname(:) = vars(:)
    nl(:) = nlev(:)

    allocate(timestamp(numtimes), datestamp(numtimes))
    call create_timestamp(t0, dt)

    ierr = v5dcreate(trim(adjustl(outfname)), numtimes, numvars, nr, nc, nl, varname, &
      timestamp, datestamp, compress, projection, proj_args, vertical, vert_args)

  end subroutine vis5d_init

  subroutine vis5d_clean()
    implicit none

    deallocate(vert_args, varname, nl, timestamp, datestamp, lon, lat, t, u, ii, jj)

  end subroutine vis5d_clean

  subroutine vis5d_write(data_in, time, var)
    use interpolate_module, only : interpolate_bilinear
    implicit none

    real(kind=sp), dimension(:,:,:), intent(in) :: data_in
    integer(kind=i4b), intent(in) :: time, var

    real(kind=dp) :: f(4)
    integer(kind=i4b) :: i, j, i0, j0, i1, j1, k, nz

    nz = nl(var)
    do k=1, nz
      do i=1, nc
        do j=1, nr
          j0 = jj(j,i)
          i0 = ii(j,i)
          i1 = modulo(i0, nc) + 1
          f(1) = data_in(ii(j,i), jj(j,i), )
          data_grid(j, i, k) = interpolate_bilinear(f, t(j,i), u(j,i))
        end do
      end do  
    end do
    ierr = v5dwrite(time, var, data_grid)

  end subroutine vis5d_write

! private routines

  subroutine set_grid_lonlat(lon0, lon1, lat0, lat1)
    use math_module, only : pi2=>math_pi2, rad2deg=>math_rad2deg
    implicit none

    real(kind=dp), dimension, intent(in) :: lon0, lon1, lat0, lat1

    integer(kind=i4b) :: i, j
    real(kind=dp) :: dlon, dlat

    dlon = (lon1-lon0)/(nc-1)
    do i=1, nc
      lon(1, i) = lon0 + dlon*(i-1)
      if (lon(1,i)>pi2) then
        lon(1,i) = lon(1,i) - pi2
      end if
    end do
    lon(2:nc,:) = lon(1,:)

    dlat = (lat1-lat0)/(nc-1)
    do j=1, nr
      lat(j, 1) = lat1 - dlat*(j-1)
    end do
    lat(:,2:nr) = lat(:,1)

    proj_args(1) = lat1*rad2deg
    proj_args(2) = lon0*rad2deg
    proj_args(3) = dlat*rad2deg
    proj_args(4) = dlon*rad2deg

  end subroutine set_grid_lonlat

  subroutine set_grid_xy(lat0, lat1)
    use math_module, only : pih=>math_pih, rad2deg=>math_rad2deg
    use earth_module, only : earth_radius
    implicit none

    real(kind=dp), intent(in) :: lat0, lat1

    integer(kind=i4b) :: i, j
    real(kind=dp) :: rmax, sinlat1, dx, dy, f

    sinlat1 = sin(lat1) 
    rmax = sqrt((1.0_dp-sinlat1)/(1.0_dp+sinlat1))
    dx = 2*rmax/(nc-1) 
    dy = 2*rmax/(nr-1)
    f = sign(1, lat0)
    proj_args(1) = lat0*rad2deg ! centre lat
    proj_args(2) = lat0*rad2deg ! centre lon
    proj_args(3) = 0.5*nr
    proj_args(4) = 0.5*nc
    proj_args(5) = 2.0*atan(sqrt(dx*dx))*earth_radius*0.001 ! km

    x = -rmax
    do i=1, nc
      y = rmax
      do j=1, nr
        lon(j,i) = math_atan2(y,x)
        lat(j,i) = f*(pih-abs(2.0_dp*atan(sqrt(x*x+y*y))))
        y = y - dy
      end do
      x = x + dx
    end do

  end subroutine set_grid_xy

  subroutine set_tu(lon_in, lat_in)
    implicit none
    real(kind=dp), dimension, intent(in) :: lon_in, lat_in

    integer(kind=i4b) :: i, j, ii, jj
    real(kind=sp) :: lon1, lat1, dlonr, dlatr
    
    nlon = size(lon_in)
    nlat = size(lat_in)
    lon1 = minval(lon_in)
    lat1 = minval(lat_in)
    dlonr = 1.0_dp/abs(lon_in(2)-lon_in(1))
    dlatr = 2.0_dp/(abs(lat_in(2)-lat_in(1))+abs(lat_in(nlat)-lat_in(nlat-1)))
    do i=1, nc
      do j=1, nr
        ii(j,i) = nint((lon(j,i)-lon1)*dlonr+1)
        jj(j,i) = nint((lat(j,i)-lat1)*dlatr+1)
      end do
    end do

  end subroutine set_tu

  subroutine create_timestamp(t0, dt)
    use time_module, only : time_yyyymmdd, time_hhmmss, operator(+)
    implicit none

    type(time_type), intent(in) :: t0, dt

    integer(kind=i4b) :: nt, i
    type(time_type) :: t

    nt = size(timestamp)

    t = t0
    do i=1, nt
      t = t + dt
      datestamp(i) = time_yyyymmdd(t)
      timestamp(i) = time_hhmmss(t)
    end do

  end subroutine create_timestamp

end module vis5d_module

