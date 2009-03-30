module besttrack_module
  use kind_module, only : i4b, dp
  use time_module, only : time_type
  implicit none

  type, public :: besttrack_type
    type(time_type) :: time
    integer(kind=i4b) :: indicator, grade, direction_r25, direction_r15
    real(kind=dp) :: lon, lat, pressure, maxwind, &
      r25, r15, r25_long, r25_short, r15_long, r15_short
  end type besttrack_type

  integer(kind=i4b), parameter :: &
    besttrack_td = 2, besttrack_ts = 3, besttrack_sts = 4, besttrack_ty = 5, &
    bestttack_l = 6

  character(len=20), dimension(:), allocatable, public :: besttrack_storm_name
  integer, dimension(:), allocatable, public :: besttrack_ndata

  private :: parse_header, parse_data
  public :: besttrack_open, besttrack_read, besttrack_close

contains

  subroutine besttrack_open(u, fname)
    implicit none

    integer(kind=i4b), intent(in) :: u
    character(len=*), intent(in) :: fname

    integer(kind=i4b) :: a, n

! count number of typhoons
    n = 0
    open(unit=u, file=fname, status="old", action="read")
    do
      read(unit=u,fmt=*,end=100) a
      if (a==66666) then
        n = n + 1
      end if
    end do
100 continue

    allocate(besttrack_storm_name(n),besttrack_ndata(n)) 

    rewind(unit=u)
    n = 1
    do
      read(unit=u,fmt=*,end=200) a
      if (a==66666) then
        backspace(unit=u)
        call parse_header(u, besttrack_storm_name(n), besttrack_ndata(n))
        n = n + 1
      end if
    end do
200 continue

  end subroutine besttrack_open

  subroutine besttrack_close(u)
    implicit none

    integer(kind=i4b), intent(in) :: u

    deallocate(besttrack_storm_name,besttrack_ndata)
    close(unit=u)

  end subroutine besttrack_close

  subroutine besttrack_read(u, t, bt)
    implicit none

    integer(kind=i4b), intent(in) :: u, t ! t typhoon id
    type(besttrack_type), dimension(:), intent(inout) :: bt

    integer(kind=i4b) :: i, j, a

    rewind(unit=u)
    do i=1, t-1
      read(unit=u,fmt=*) a
      do j=1, besttrack_ndata(i)
        read(unit=u,fmt=*) a
      end do
    end do
    read(unit=u,fmt=*) a
    do j=1, besttrack_ndata(t)
      call parse_data(u, bt(j))
    end do

  end subroutine besttrack_read

! private routines

  subroutine parse_header(u, nm, nd)
    implicit none

    integer(kind=i4b), intent(in) :: u
    character(len=20), intent(out) :: nm ! name of the storm
    integer(kind=i4b), intent(out) :: nd
    
    integer(kind=i4b) :: a, b, c, d, e, f, g, i
    character(len=20) :: h

    read(unit=u,fmt=*) a, b, c, d, e, f, g, h, i
    nm = adjustl(h)
    nd = c

  end subroutine parse_header

  subroutine parse_data(u, bt)
    use math_module, only : math_nm2m, math_knot2ms
    use string_module, only : string_digits
    implicit none

    integer(kind=i4b), intent(in) :: u
    type(besttrack_type), intent(inout) :: bt

    integer(kind=i4b), dimension(11) :: x

    x = 0
    read(unit=u,fmt=*) x(1:7)
    x(1) = 2000000000+x(1)
    bt % time = time_type( &
      string_digits(x(1), 1, 4), &
      string_digits(x(1), 5, 6), &
      string_digits(x(1), 7, 8), &
      string_digits(x(1), 9, 10), 0, 0.0)
    bt % indicator = x(2)
    bt % grade = x(3)
    bt % lat = 0.1_dp*x(4)
    bt % lon = 0.1_dp*x(5)
    bt % pressure = 100.0_dp*x(6)
    bt % maxwind = x(7)*math_knot2ms

    if ((x(3)==besttrack_ts).or. &
        (x(3)==besttrack_sts).or.&
        (x(3)==besttrack_ty)) then
      backspace(unit=u)
      read(unit=u,fmt=*) x
      if (x(8)/=0) then
        bt % direction_r25 = string_digits(x(8),1,1)
        bt % r25_long = string_digits(x(8),2,5)*math_nm2m
        bt % r25_short = x(9)*math_nm2m
        bt % r25 = 0.5_dp*(bt % r25_long+bt % r25_short)
      end if
      bt % direction_r15 = string_digits(x(10),1,1)
      bt % r15_long = string_digits(x(10),2,5)*math_nm2m
      bt % r15_short = x(11)*math_nm2m
      bt % r15 = 0.5_dp*(bt % r15_long+bt % r15_short)
    end if

  end subroutine parse_data

end module besttrack_module
