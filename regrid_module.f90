module regrid_module

! interpolate values in a given lon-lat grid

  use kind_module, only: i4b, dp
  use math_module, only: pi=>math_pi, pi2=>math_pi2, pih=>math_pih
  use search_module, only: search_bisection
  implicit none
  private

  integer(kind=i4b), private :: nx, nxh, ny, n, m
  real(kind=dp), private :: dlon, dlonr
  real(kind=dp), dimension(:), allocatable, public :: regrid_lon, regrid_lat

  public :: regrid_init, regrid_clean, regrid_extend, &
    regrid_bilinear, regrid_linpol, regrid_spcher, regrid_bicubic
  private :: handle_beyondpole

contains

  subroutine regrid_init(lon,lat,nn,mm)
    use glatwgt_module, only: glatwgt_approx
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, lat
    integer(kind=i4b), intent(in) :: nn, mm

    integer(kind=i4b) :: i, j
    real(kind=dp) :: dlat

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
    call glatwgt_approx(lat1, dlat, ny)
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

  function regrid_bilinear(b,lon,lat,bpole) result(ai)
    use interpolate_module, only: interpolate_bilinear
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: b
    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp), intent(in), optional :: bpole
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j
    real(kind=dp) :: t, u, lat1
    real(kind=dp), dimension(4) :: f
    logical :: lpole

    lpole = .false.
    lat1=lat
    i = floor(lon*dlonr+1) + n
    j = search_bisection(regrid_lat,lat)
    if (present(bpole))
      call handle_beyondpole(lat1,i,j)
    end if
    lpole = present(bpole).and.((j==m).or.(j==ny+m))
    t = (lon-regrid_lon(i))*dlonr
    u = (lat1-regrid_lat(j))/(regrid_lat(j+1)-regrid_lat(j))
    f(1) = b(i,  j  )
    f(2) = b(i+1,j  )
    if (lpole)
      f(3) = bpole
      f(4) = bpole
      u = u*2.0_dp
    else
      f(3) = b(i+1,j+1)
      f(4) = b(i,  j+1)
    end if
    ai = interpolate_bilinear(f,t,u)

  end function regrid_bilinear

  function regrid_linpol(b,lon,lat,usequasicubic,bpole) result(ai)
    use interpolate_module, only: interpolate_linpol
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: b
    real(kind=dp), intent(in) :: lon, lat
    logical, intent(in), optional :: usequasicubic
    real(kind=dp), intent(in), optional :: bpole
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j, ii, jj
    real(kind=dp) :: t, lat1
    real(kind=dp), dimension(4) :: lon4, lat4
    real(kind=dp), dimension(4,4) :: f
    logical :: lpole

    lpole = .false.
    lat1 = lat
    i = floor(lon*dlonr+1) + n
    j = search_bisection(regrid_lat,lat)
    if (present(bpole))
      call handle_beyondpole(lat1,i,j)
    end if
    lpole = present(bpole).and. &
      (((m<=j).and.(j<=m+1)).or.((ny+m-1<=j).and.(j<=ny+m))) then
    lon4(:) = regrid_lon(i-1:i+2)
    if (.not.lpole) then
      lat4(:) = regrid_lat(j-1:j+2)
      f(:,:) = b(i-1:i+2,j-1:j+2)
    else
      if (j==m) then
        lat4(1) = regrid_lat(j)
        lat4(2) = 0.5_dp*(regrid_lat(j)+regrid_lat(j+1)) ! -pih or pih
        lat4(3:4) = regrid_lat(j+1:j+2)
        f(:,1) = -b(i-1:i+2,j)
        f(:,2) = bpole
        f(:,3:4) = b(i-1:i+2,j+1:j+2)
      else if (j==m+1) then
        lat4(1) = 0.5_dp*(regrid_lat(j-1)+regrid_lat(j+1)) ! -pih or pih
        lat4(2:4) = regrid_lat(j:j+2)
        f(:,1) = bpole
        f(:,2:4) = b(i-1:i+2,j:j+2)
      else if (j==ny+m-1) then
        lat4(1:3) = regrid_lat(j-1:j+1)
        lat4(4) = 0.5_dp*(regrid_lat(j+1)+regrid_lat(j+2)) ! -pih or pih
        f(:,1:3) = b(i-1:i+2,j-1:j+1)
        f(:,4) = bpole
      else if (j==ny+m) then
        lat4(1:2) = regrid_lat(j-1:j)
        lat4(3) = 0.5_dp*(regrid_lat(j)+regrid_lat(j+1)) ! -pih or pih
        lat4(4) = regrid_lat(j+1)
        f(:,1:2) = b(i-1:i+2,j-1:j)
        f(:,3) = bpole
        f(:,4) = -b(i-1:i+2,j+1)
      end if
    end if
    t = -1
    if (present(usequasicubic).and.usequasicubic) then
      t = (lon-regrid_lon(i))*dlonr
    end if
    ai = interpolate_linpol(f,lon4,lat4,lon,lat,t)

  end function regrid_linpol

  function regrid_spcher(b,bx,lon,lat,usecubic,bpole) result(ai)
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
    j = search_bisection(regrid_lat,lat)
    if (present(bpole))
      call handle_beyondpole(lat1,i,j)
    end if
    lpole = present(bpole).and. &
      (((m<=j).and.(j<=m+2)).or.((ny+m-2<=j).and.(j<=ny+m))) then
    if (.not.lpole) then
      lat6(:) = regrid_lat(j-2:j+3)
      f(:,:)  =  b(i:i+1,j-2:j+3)
      fx(:,:) = bx(i:i+1,j-2:j+3)
    else
      if (j==m) then
        lat6(1:2) = regrid_lat(j-1:j)
        lat6(3) = 0.5_dp*(regrid_lat(j)+regrid_lat(j+1))
        lat6(4:6) = regrid_lat(j+1:j+3)
        f(:,1:2) = -b(i:i+1,j-1,j)
        f(:,3) = bpole
        f(:4:6) = regrid_lat(j+1:j+3)
        fx(:,1:2) = bx(i:i+1,j-1,j)
        fx(:,3) = 0.0_dp
        fx(:4:6) = regrid_lat(j+1:j+3)
      else (j==m+1) then
        lat6(1) = regrid_lat(j-1) 
        lat6(2) = 0.5_dp*(regrid_lat(j-1)+regrid_lat(j))
        lat6(3:6) = regrid_lat(j:j+3)
        f(:,1) = -b(i:i+1,j-1)
        f(:,2) = bpole
        f(:,3:6) = b(i:i+1,j:j+3)
        fx(:,1) = -bx(i:i+1,j-1)
        fx(:,2) = 0.0_dp
        fx(:,3:6) = bx(i:i+1,j:j+3)
      else (j==m+2) then
        lat6(1) = 0.5_dp*(regrid_lat(j-2)+regrid_lat(j-1))
        lat6(2:6) = regrid_lat(j-1:j+3)
        f(:,1) = bpole
        f(:,2:6) = bx(i:i+1,j-1:j+3)
        fx(:,1) = 0.0_dp
        fx(:,2:6) = bx(i:i+1,j-1:j+3)
      else (j=ny+m-2) then
        lat6(1:5) = regrid_lat(j-2:j+2)
        lat6(6) = 0.5_dp*(regrid_lat(j+2)+regrid_lat(j+3))
        f(:,1:5) = b(i:i+1,j-2:j+2)
        f(:,6) = bpole
        fx(:,1:5) = b(i:i+1,j-2:j+2)
        fx(:,6) = 0.0_dp
      else (j=ny+m-1) then
        lat6(1:4) = regrid_lat(j-2:j+1)
        lat6(5) = 0.5_dp*(regrid_lat(j+1)+regrid_lat(j+2))
        lat6(6) = regrid_lat(j+2)
        f(:,1:4) = f(i:i+1,j-2:j+1)
        f(:,5) = bpole
        f(:,6) = -b(i:i+1,j+2)
        fx(:,1:4) = bx(i:i+1,j-2:j+1)
        fx(:,5) = 0.0_d0
        fx(:,6) = bx(i:i+1,j+2)
      else (j=ny+m) then
        lat(1:3) = regrid_lat(j-2:j)
        lat(4) = 0.5_dp*(regrid_lat(j)+regrid_lat(j+1))
        lat(5:6) = regrid_lat(j+1:j+2)
        f(:,1:3) = b(i:i+1,j-2:j)
        f(:,4) = bpole
        f(:,5:6) = -b(i:i+1,j+1:j+2)
        fx(:,1:3) = bx(i:i+1,j-2:j)
        fx(:,4) = 0.0_dp
        fx(:,5:6) = bx(i:i+1,j+1:j+2)
      end if
    end if
    if (present(usecubic).and.usecubic) then
      ai = interpolate_spcher(f(:,2:5),fx(:,2:5),lat6(2:5),dlon,lat,t)
    else
      ai = interpolate_spcher(f,fx,lat6,dlon,lat,t)
    end if

  end function regrid_spcher

  function regrid_bicubic(b,bx,by,bxy,lon,lat) result(ai)
    use interpolate_module, only: interpolate_bicubic
    implicit none
 
    real(kind=dp), dimension(:,:), intent(in) :: b, bx, by, bxy
    real(kind=dp), intent(in) :: lon, lat
    real(kind=dp) :: ai

    integer(kind=i4b) :: i, j
    real(kind=dp) :: t, u, dlat
    real(kind=dp), dimension(4) :: f, fx, fy, fxy

    i = floor(lon*dlonr+1) + n
    j = search_bisection(regrid_lat,lat)
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

! private routines

  subroutine handle_beyondpole(lat,i,j)
    implicit none

    real(kind=dp), intent(inout) :: lat
    integer(kind=i4b), intent(inout) :: i
    integer(kind=i4b), intent(in) :: j

    if ((j==0).or.(j==ny+2*m)) then
      print *, "Error in regrid. Out of bounds."
      stop
    end if

    if (abs(lat)>pih) then
      lat = sign(1.0_dp,lat)*pi-lat
      i = modulo(i-n+nx/2-1, nx)+n+1
    end if

  end subroutine handle_beyondpole(lat,i,j)

end module regrid_module
  
