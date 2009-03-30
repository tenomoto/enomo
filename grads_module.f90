module grads_module
  use kind_module, only : i4b, sp
  use string_module, only : string_tolower
  implicit none
  private

  type, public :: grads_ctl
    character(len=256) :: dset, title
    real(kind=sp) :: undef
    integer(kind=i4b) :: xdef, ydef, zdef, tdef
    integer(kind=i4b), dimension(6) :: t
    character(len=2) :: tunit
    real(kind=sp), dimension(:), allocatable :: xlevs, ylevs, zlevs
    integer(kind=i4b) :: vars
    integer(kind=i4b), dimension(:), allocatable :: levs
    character(len=15), dimension(:), allocatable :: varnames
  end type grads_ctl

  public :: grads_parse, grads_clean
  private :: read_dset, read_title, read_undef, read_options, &
    read_axis, read_tdef, read_vars, read_var
  
contains

  function grads_parse(u, ctl) result(f)
    implicit none

    character(len=*), intent(in) :: ctl
    integer(kind=i4b), intent(in) :: u
    type(grads_ctl) :: f

    integer(kind=i4b) :: i, n
    character(len=16) :: a

    open(unit=u, file=ctl, status="old", action="read")
    do
      read(unit=u, fmt=*, end=100) a
      a = string_tolower(adjustl(a))
      select case(trim(a))
        case("dset")
          f % dset = read_dset(u)
        case("title")
          f % title = read_title(u)
        case("undef")
          f % undef = read_undef(u)
        case("xdef")
          f % xdef = read_axis(u, f % xlevs)
        case("ydef")
          f % ydef = read_axis(u, f % ylevs)
        case("zdef")
          f % zdef = read_axis(u, f % zlevs)
        case("vars")
          n = read_vars(u) 
          f % vars = n
          allocate(f % levs(n), f % varnames(n))
          do i=1, n
            f % levs(i) = read_var(u, f % varnames(i))
          end do
        case default
      end select
    end do
100 continue
    close(unit=u)

  end function grads_parse

  subroutine grads_clean(f)
    implicit none

    type(grads_ctl), intent(inout) :: f

    deallocate(f % xlevs, f % ylevs, f % zlevs, f % levs, f % varnames)

  end subroutine grads_clean

  function read_dset(u) result(dset)
    implicit none

    integer(kind=i4b), intent(in) :: u
    character(len=256) :: dset

    character(len=4) :: a

    backspace(unit=u)
    read(unit=u, fmt=*) a, dset
    if (dset(1:1)=="^") then
      dset(1:1) = " "
    end if
    dset = adjustl(dset)

  end function read_dset

  function read_title(u) result(title)
    implicit none

    integer(kind=i4b), intent(in) :: u
    character(len=256) :: title

    character(len=5) :: a

    backspace(unit=u)
    read(unit=u, fmt=*) a, title
    title = adjustl(title)

  end function read_title

  function read_undef(u) result(undef)
    implicit none

    integer(kind=i4b), intent(in) :: u
    real(kind=sp) :: undef

    character(len=5) :: a

    backspace(unit=u)
    read(unit=u, fmt=*) a, undef

  end function read_undef

  function read_options(u,n) result(s)
    implicit none

    integer(kind=i4b), intent(in) :: u, n
    character(len=16), dimension(:), allocatable :: s

    character(len=6) :: a

    backspace(unit=u)
    allocate(s(n))
    read(unit=u, fmt=*) a, s

  end function read_options

  function read_axis(u,axis) result(n)
    implicit none

    integer(kind=i4b), intent(in) :: u
    real(kind=sp), dimension(:), allocatable, intent(out) :: axis
    integer(kind=i4b) :: n

    character(len=4) :: a
    character(len=6) :: b
    integer(kind=i4b) :: i
    real(kind=sp) :: x0, dx

    backspace(unit=u)
    read(unit=u, fmt=*) a, n, b
    b=string_tolower(b)
    allocate(axis(n))
    backspace(unit=u)
    if (b=="linear") then
      read(unit=u, fmt=*) a, n, b, x0, dx
      do i=1, n
        axis(i) = x0 + dx*(i-1)
      end do
    else
      read(unit=u, fmt=*) a, n, b, axis
    end if

  end function read_axis

  function read_tdef(u,t,tu) result(n)
    implicit none

    integer(kind=i4b), intent(in) :: u
    integer(kind=i4b), dimension(6), intent(out) :: t
    character(len=2), intent(out) :: tu
    integer(kind=i4b) :: n

    character(len=*), parameter :: &
      monthname = "jan feb mar apr may jun jul aug sep oct nov dec", &
      tunits = "mnhrdymoyr"

    character(len=4) :: a
    character(len=6) :: b
    character(len=15) :: c
    character(len=4) :: d
    integer(kind=i4b) :: ic, iz, im, ik

    backspace(unit=u)
    read(unit=u, fmt=*) a, n, b, c, d
    c = adjustl(string_tolower(c))
    ic = index(c,":")
    iz = index(c,"z")
    im = iz+scan(c(iz+1:),monthname)
    read(c(im+3:),fmt=*) t(1)
    if (t(1)<100) then
      if (t(1)>=50) then
        t(1) = t(1) + 1900
      else
        t(1) = t(1) + 2000
      end if
    end if
    t(2) = (index(monthname,c(im:im+2))+3)/4
    read(c(iz+1:im-1),fmt=*) t(3)
    if (ic/= 0) then
      read(c(1:ic-1),fmt=*) t(4)
      read(c(ic+1:iz-1),fmt=*) t(5)
    else
      t(4:5) = 0
    end if
    d = adjustl(d)
    ik = scan(d,tunits)
    read(d(1:ik-1),fmt=*) t(6)
    tu = d(ik:ik+1)
    
  end function read_tdef

  function read_vars(u) result(n)
    implicit none

    integer(kind=i4b), intent(in) :: u
    integer(kind=i4b) :: n

    character(len=4) :: a

    backspace(unit=u)
    read(unit=u, fmt=*) a, n

  end function read_vars

  function read_var(u,s) result(n)
    implicit none

    integer(kind=i4b), intent(in) :: u
    character(len=15), intent(out) :: s
    integer(kind=i4b) :: n

    read(unit=u, fmt=*) s, n

  end function read_var

end module grads_module
