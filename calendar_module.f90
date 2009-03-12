module calendar_module
! provides date utility fuctions
  use type_module, only: i4b
  implicit none

  character(len=*), dimension(12), parameter, public :: calendar_monthnames = &
    (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  character(len=*), dimension(0:6), parameter, public :: calendar_weekdaynames = &
    (/"Sun","Mon","Tue","Wed","Thr","Fri","Sat"/)

  public :: calendar_isleapyear, calendar_daysinmonth, &
            calendar_weekday, calendar_dayofyear, calendar_date, &
            calendar_pentad, calendar_pentaddates, calendar_ndays

contains

  function calendar_isleapyear(y) result(d)
!   returns d=1 if argument y is a leap year otherwise d=0
    implicit none

    integer(kind=i4b), intent(in) :: y
    integer(kind=i4b) :: d

    d = 0
    if (modulo(y,4)==0) then
      if (modulo(y,100)/=0) then
        d = 1
      else
        if (modulo(y,400)==0) then
          d = 1
        end if
      end if
    end if

  end function calendar_isleapyear

  function calendar_daysinmonth(y,m) result(d)
!   returns the number of days in the month
    implicit none

    integer(kind=i4b), intent(in) :: y, m
    integer(kind=i4b) :: d

    select case(m)
      case(4,6,9,11)
        d = 30
      case(2)
        d = 28 + calendar_isleapyear(y)
      case default
        d = 31
    end select

  end function calendar_daysinmonth

  function calendar_weekday(y,m,d) result(w)
!   returns the day of the week using Zeller's algorithm
!   0-6 correspond to Sun-Sat
    implicit none

    integer(kind=i4b), intent(in) :: y, m, d
    integer(kind=i4b) :: x,n,w

    x = y
    n = m
! treat Jan and Feb as 13th and 14th months of the previous year
    if (m<=2) then
      n = n + 12
      x = x - 1
    end if

    w = modulo((d-1+(13*(n+1)/5)+x+x/4-x/100+x/400),7)

  end function calendar_weekday

  function calendar_dayofyear(y,m,d) result(n)
! returns the day of the year
    implicit none

    integer(kind=i4b), intent(in) :: y, m, d
    integer(kind=i4b) :: i, n

    n = d
    do i=1, m-1
      n = n + calendar_daysinmonth(y,i)
    end do

  end function calendar_dayofyear

  function calendar_date(y,n) result(ymd)
!   returns year, month, date from the year and day of the year
    implicit none

    integer(kind=i4b), intent(in) :: y, n
    integer(kind=i4b) :: x, m, d, k
    integer(kind=i4b), dimension(3) :: ymd

    x = y
    m = 1
    d = n
    substloop: do
      k = calendar_daysinmonth(x,m)
      if (d <= k) then
        exit substloop
      else
        if (m<12) then
          m = m + 1
        else
          x = x + 1
          m = 1
        end if
        d = d - k
      end if
    end do substloop
    ymd(1) = x
    ymd(2) = m
    ymd(3) = d

  end function calendar_date

  function calendar_pentad(y,m,d) result(p)
!   returns pentad the date(y,m,d) belogs to
    implicit none

    integer(kind=i4b), intent(in) :: y, m, d
    integer(kind=i4b) :: l,x,p
    integer(kind=i4b), parameter :: z = 60 ! 29 Feb

    l = calendar_isleapyear(y)
    x = calendar_dayofyear(y,m,d)
    if (x>z) then
      x = x - l
    end if
    p = (x-1)/5 + 1

  end function calendar_pentad

  function calendar_pentaddates(p) result(n)
!   returns the beginning and ending dates of the pentad
    implicit none

    integer(kind=i4b), intent(in) :: p
    integer(kind=i4b), dimension(4) :: n
    integer(kind=i4b), dimension(3) :: b, e
    integer(kind=i4b) :: d
    integer(kind=i4b), parameter :: y = 1999 ! non-leap dummy year

    d = (p-1)*5+1

    b = calendar_date(y,d) 
    e = calendar_date(y,d+4) 
    
    n(1:2) = b(2:3)
    n(3:4) = e(2:3)

  end function calendar_pentaddates

  function calendar_ndays(y1,m1,d1,y2,m2,d2) result(n)
! returns number of days
    implicit none

    integer(kind=i4b), intent(in) :: y1, m1, d1, y2, m2, d2
    integer(kind=i4b) :: n

    integer(kind=i4b) :: n1, n2, m, ny, y

    n1 = calendar_dayofyear(y1,m1,d1)
    n2 = calendar_dayofyear(y2,m2,d2)
    ny = y2-y1
    m = 0
    do y=y1, y2-1
      m = m + calendar_isleapyear(y)
    end do
    n = max(0,n2 + ny*365 + m - n1 + 1)
 
  end function calendar_ndays

end module calendar_module
