module string_module
  use kind_module, only : i2b, i4b, i8b, sp, dp
  implicit none

  public :: string_tolower, string_toint, string_tolong, string_tofloat, string_todouble, &
    string_digits

contains

  function string_tolower(s) result(l)
    implicit none

    character(len=*), intent(in) :: s
    character(len=len(s)) :: l

    integer(kind=i4b), parameter :: cap_a = 65, cap_z = 90, off = 32
    integer(kind=i4b) :: i, a

    l = adjustl(s)
    do i=1, len_trim(l)
      a = iachar(l(i:i))
      if ((a>=cap_a).and.(a<=cap_z)) then
        a = a + off
        l(i:i) = achar(a)
      end if
    end do

  end function string_tolower

  function string_toshort(s) result(i)
    implicit none

    character(len=*), intent(in) :: s
    integer(kind=i2b) :: i

    read(s,fmt=*) i

  end function string_toshort

  function string_toint(s) result(i)
    implicit none

    character(len=*), intent(in) :: s
    integer(kind=i4b) :: i

    read(s,fmt=*) i

  end function string_toint

  function string_tolong(s) result(i)
    implicit none

    character(len=*), intent(in) :: s
    integer(kind=i8b) :: i

    read(s,fmt=*) i

  end function string_tolong

  function string_tofloat(s) result(x)
    implicit none

    character(len=*), intent(in) :: s
    real(kind=sp) :: x

    read(s,fmt=*) x

  end function string_tofloat

  function string_todouble(s) result(x)
    implicit none

    character(len=*), intent(in) :: s
    real(kind=dp) :: x

    read(s,fmt=*) x

  end function string_todouble

  function string_digits(a, i, j) result(b)
    implicit none

    integer(kind=i4b), intent(in) :: a, i, j
    integer(kind=i4b) :: b

    character(len=16) astr, bstr

    write(astr, *) a
    astr = adjustl(astr)
    bstr = astr(i:j)
    read(bstr, *) b

  end function string_digits

end module string_module
