module regrid_module

! interpolate values in a given lon-lat grid

  use type_module, only: i4b, dp
  implicit none
  private

  integer(kind=i4b), private :: nx, nxh, ny, n, m
  real(kind=dp), private :: dlon, dlonr
  real(kind=dp), dimension(:), allocatable, public :: regrid_lon, regrid_lat

  public :: regrid_init, regrid_clean, regrid_extend, &
    regrid_bilinear, regrid_linpol, regrid_spcher, regrid_bicubic

contains

  subroutine regrid_init(lon,lat,nn,mm)
    use math_module, only: pi=>math_pi, pi2=>math_pi2
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    integer(kind=i4b), intent(in) :: nn, mm

    integer(kind=i4b) :: i, j

    nx  = size(lon)
    ny  = size(lat)
    nxh = nx/2
    n = nn
    m = mm
    allocate(regrid_lon(nx+2*n),regrid_lat(ny+2*m))

    dlon = lon(2)-lon(1)
    dlonr = 1.0_dp/dlon
    do i=1, nx
      regrid_lon(i+n) = lon(i)
    end do
    do i=1, n
      regrid_lon(i) = -(n-i+1)*dlon
      regrid_lon(nx+n+i) = pi2+(i-1)*dlon
    end do
    do j=1, ny
      regrid_lat(j+m) = lat(j)
    end do
    do j=1, m
      regrid_lat(j)        = sign(pi,lat(1)) - lat(m-j+1)
      regrid_lat(ny+2*m-j+1) = -regrid_lat(j)
    end do

  end subroutine regrid_init

  subroutine regrid_clean()
    implicit none

    deallocate(regrid_lon,regrid_lat)

  end subroutine regrid_clean

  subroutine regrid_extend(a, b, lreverse)
    implicit none

  ! copy a to b and sets halo regions
  ! assumes cyclic boundary condition in x
  !         180 degree shift beyond poles in y
  ! values beyond poles are multiplied by -1 if lreverse

    real(kind=dp), dimension(:,:), intent(in) :: a
    real(kind=dp), dimension(:,:), intent(inout) :: b
    logical, optional, intent(in) :: lreverse

    integer(kind=i4b) :: i, j
    real(kind=dp) :: d

    if (present(lreverse).and.lreverse) then
      d = -1.0_dp
    else
      d =  1.0_dp
    end if

    ! copy input into buffer
    do j=1, ny
      do i=1, nx
        b(i+n,j+m) = a(i,j)
      end do
    end do
    ! fill each m lat beyond poles with 180 degree shift
    do j=1, m
      do i=1, nxh
        b(i+nxh+n,j)      = d*a(i,    -j+m+1)
        b(i+n,    j)      = d*a(i+nxh,-j+m+1)
        b(i+nxh+n,j+ny+m) = d*a(i,    -j+ny+1)
        b(i+n,    j+ny+m) = d*a(i+nxh,-j+ny+1)
      end do
    end do
    ! cyclic condition in x
    do j=1, ny+2*m
      do i=1, n
        b(i,     j) = b(i+nx,j)
        b(i+nx+n,j) = b(i+n, j)
      end do
    end do

  end subroutine regrid_extend

  function regrid_bilinear(b,lon,lat) result(ai)
    use glatwgt_module, only: glatwgt_within
    use interpolate_module, only: interpolate_bilinear
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: b
    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j
    real(kind=dp) :: t, u
    real(kind=dp), dimension(4) :: f

    i = floor(lon*dlonr+1) + n
    j = glatwgt_within(regrid_lat(1+m:ny+m),lat) + m
    f(1) = b(i,  j  )
    f(2) = b(i+1,j  )
    f(3) = b(i+1,j+1)
    f(4) = b(i,  j+1)
    t = (lon-regrid_lon(i))*dlonr
    u = (lat-regrid_lat(j))/(regrid_lat(j+1)-regrid_lat(j))
    ai = interpolate_bilinear(f,t,u)

  end function regrid_bilinear

  function regrid_linpol(b,lon,lat,wolin) result(ai)
    use glatwgt_module, only: glatwgt_within
    use interpolate_module, only: interpolate_linpol
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: b
    real(kind=dp), intent(in) :: lon, lat
    logical, intent(in), optional :: wolin
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j, ii, jj
    real(kind=dp) :: t
    real(kind=dp), dimension(4) :: lon4, lat4
    real(kind=dp), dimension(4,4) :: f

    i = floor(lon*dlonr+1) + n
    j = glatwgt_within(regrid_lat(1+m:ny+m),lat) + m
    lon4(:) = regrid_lon(i-1:i+2)
    lat4(:) = regrid_lat(j-1:j+2)
    f(:,:) = b(i-1:i+2,j-1:j+2)
    if (present(wolin).and.wolin) then
      t = -1
    else
      t = (lon-regrid_lon(i))*dlonr
    end if
    ai = interpolate_linpol(f,lon4,lat4,lon,lat,t)

  end function regrid_linpol

  function regrid_spcher(b,bx,lon,lat,usecubic) result(ai)
    use glatwgt_module, only: glatwgt_within
    use interpolate_module, only: interpolate_spcher
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: b,bx
    real(kind=dp), intent(in) :: lon, lat
    logical, intent(in), optional :: usecubic
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j, ii, jj
    real(kind=dp) :: t
    real(kind=dp), dimension(6) :: lat6
    real(kind=dp), dimension(2,6) :: f, fx

    i = floor(lon*dlonr+1) + n
    t = (lon-regrid_lon(i))*dlonr
    j = glatwgt_within(regrid_lat(1+m:ny+m),lat) + m
    lat6(:) = regrid_lat(j-2:j+3)
    f(:,:)  =  b(i:i+1,j-2:j+3)
    fx(:,:) = bx(i:i+1,j-2:j+3)
    if (present(usecubic).and.usecubic) then
      ai = interpolate_spcher(f(:,2:5),fx(:,2:5),lat6(2:5),dlon,lat,t)
    else
      ai = interpolate_spcher(f,fx,lat6,dlon,lat,t)
    end if

  end function regrid_spcher

  function regrid_bicubic(b,bx,by,bxy,lon,lat) result(ai)
    use glatwgt_module, only: glatwgt_within
    use interpolate_module, only: interpolate_bicubic
    implicit none
 
    real(kind=dp), dimension(:,:), intent(in) :: b, bx, by, bxy
    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j
    real(kind=dp) :: t, u, dlat
    real(kind=dp), dimension(4) :: f, fx, fy, fxy

    i = floor(lon*dlonr+1) + n
    j = glatwgt_within(regrid_lat(1+m:ny+m),lat) + m
    dlat = regrid_lat(j+1)-regrid_lat(j)
    f(1)   =   b(i,  j  )
    f(2)   =   b(i+1,j  )
    f(3)   =   b(i+1,j+1)
    f(4)   =   b(i,  j+1)
    fx(1)  =  bx(i,  j  )
    fx(2)  =  bx(i+1,j  )
    fx(3)  =  bx(i+1,j+1)
    fx(4)  =  bx(i,  j+1)
    fy(1)  =  by(i,  j  )
    fy(2)  =  by(i+1,j  )
    fy(3)  =  by(i+1,j+1)
    fy(4)  =  by(i,  j+1)
    fxy(1) = bxy(i,  j  )
    fxy(2) = bxy(i+1,j  )
    fxy(3) = bxy(i+1,j+1)
    fxy(4) = bxy(i,  j+1)
    t = (lon-regrid_lon(i))*dlonr
    u = (lat-regrid_lat(j))/dlat
    ai = interpolate_bicubic(f,fx,fy,fxy,dlon,dlat,t,u)

  end function regrid_bicubic

end module regrid_module
  
