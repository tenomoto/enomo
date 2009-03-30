module time_module
  use kind_module, only : i4b, sp
  implicit none
  private

  type, public :: time_type
    integer(kind=i4b) :: year, month, day, hour, minute
    real(kind=sp) :: second
  end type time_type

  private :: add_time
  public :: time_yyyymm, time_yyyymmdd, time_yyyymmddhh, time_hhmm, time_hhmmss, operator(+)

  interface operator(+)
    module procedure add_time
  end interface

contains

  function time_yyyymm(time) result(t)
    implicit none

    type(time_type), intent(in) :: time
    integer(kind=i4b) :: t

    t = time % year * 100 + time % month

  end function time_yyyymm 

  function time_yyyymmdd(time) result(t)
    implicit none

    type(time_type), intent(in) :: time
    integer(kind=i4b) :: t

    t = time_yyyymm(time)
    t = t * 100 + time % day

  end function time_yyyymmdd 

  function time_yyyymmddhh(time) result(t)
    implicit none

    type(time_type), intent(in) :: time
    integer(kind=i4b) :: t

    t = time_yyyymmdd(time)
    t = t * 100 + time % hour

  end function time_yyyymmddhh

  function time_hhmm(time) result(t)
    implicit none

    type(time_type), intent(in) :: time
    integer(kind=i4b) :: t

    t = time % hour * 100 + time % minute

  end function time_hhmm

  function time_hhmmss(time) result(t)
    implicit none

    type(time_type), intent(in) :: time
    integer(kind=i4b) :: t

    t = time % hour * 10000 + time % minute * 100 + nint(time % second)

  end function time_hhmmss

  function add_time(t1, t2) result (t)
    use calendar_module, only : calendar_dayofyear, calendar_date
    implicit none

    type(time_type), intent(in) :: t1, t2
    type(time_type) :: t

    integer(kind=i4b) :: inc, doy
    integer(kind=i4b), dimension(3) :: ymd

    t % second = t1 % second + t2 %second
    inc = floor(t % second / 60.0)
    t % second = t % second - inc*60.0

    t % minute = t1 % minute + t2 % minute + inc
    inc = t % minute / 60
    t % minute = modulo(t % minute, 60)

    t % hour = t1 % hour + t2 % hour + inc
    inc = t % hour / 24
    t % hour = modulo(t % hour, 24)

    t % day = t1 % day + t2 % day + inc

    t % month = t1 % month + t2 % month
    inc = (t % month-1) / 12
    t % month = modulo((t % month-1), 12) + 1

    t % year = t1 % year + t2 % year + inc

    ymd = calendar_date(t % year, &
      calendar_dayofyear(t % year, t % month, t % day))

    t % year  = ymd(1) 
    t % month = ymd(2) 
    t % day   = ymd(3) 

  end function add_time

end module time_module
