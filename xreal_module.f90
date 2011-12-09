module xreal_module
! Extended exponent of floating-point numbers

! Reference:
!   Fukushima, Toshio, 2011: Numerical computation of spherical
!     harmonics of arbitrary degree and order by extending
!     exponent of floating point numbers. J. Geodesy,
!     doi:10.1007//s00190-011-0519-2

  use kind_module, only: dp, i4b
  implicit none
  private

  integer(kind=i4b), parameter :: ind = 960, indh = ind/2
  real(kind=dp), parameter, public :: &
    xreal_big  = 2.0_dp**ind,  xreal_bigi  = 2.0_dp**(-ind), &
    xreal_bigs = 2.0_dp**indh, xreal_bigsi = 2.0_dp**(-indh)

  type xreal_type
    real(kind=dp) :: p     ! principal component
    integer(kind=i4b) :: i ! auxiliary index
  end type xreal_type

  interface assignment(=)
    module procedure x_assign_f, f_assign_x
  end interface

  interface operator(+)
    module procedure xadd
  end interface

  interface operator(-)
    module procedure xsub
  end interface

  interface operator(*)
    module procedure xmul, fx
  end interface

  interface operator(/)
    module procedure xdiv
  end interface

  public :: xreal_type, xreal_norm, xreal_fxpgy, xreal_base10, xreal_test, &
    assignment(=), operator(+), operator(-), operator(*), operator(/)
  private :: x_assign_f, f_assign_x ,xmul, xdiv, xadd, xsub

contains

  subroutine x_assign_f(x,f)
    type(xreal_type), intent(out) :: x
    real(kind=dp), intent(in) :: f

    x%p = f
    x%i = 0

  end subroutine x_assign_f

  subroutine f_assign_x(g,x)
    real(kind=dp), intent(out) :: g
    type(xreal_type), intent(in) :: x

    select case(x%i)
      case(0)
        g = x%p
      case(:-1) ! underflow
        g = x%p*xreal_bigi
      case(1:)  ! overflow
        g = x%p*xreal_big
    end select

  end subroutine f_assign_x

  function xreal_norm(x) result(y)
    type(xreal_type), intent(in) :: x
    type(xreal_type) :: y

    real(kind=dp) :: w

    w = abs(x%p)
    if (w>=xreal_bigs) then
      y%p = x%p*xreal_bigi
      y%i = x%i + 1
    else if (w<xreal_bigsi) then
      y%p = x%p*xreal_big
      y%i = x%i - 1
    else
      y = x
    end if

  end function xreal_norm

  function xmul(x,y) result(z)
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: s, z

    z = xreal_norm(x)
    s = xreal_norm(y)
    z%p = z%p * s%p
    z%i = z%i + s%i
    z = xreal_norm(z)

  end function xmul

  function fx(f,x) result(z)
    real(kind=dp), intent(in) :: f
    type(xreal_type), intent(in) :: x
    type(xreal_type) :: z

    z = xreal_norm(x)
    z%p = f * z%p
    z%i = z%i
    z = xreal_norm(z)

  end function fx

  function xdiv(x,y) result(z)
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: s, z

    z = xreal_norm(x)
    s = xreal_norm(y)
    z%p = z%p / s%p
    z%i = z%i - s%i
    z = xreal_norm(z)

  end function xdiv

  function xreal_fxpgy(f,x,g,y) result(z)
    real(kind=dp), intent(in) :: f, g
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: z

    integer(kind=i4b) :: id

    id = x%i - y%i
    select case(id)
      case(0)
        z%p = f*x%p + g*y%p
        z%i = x%i
      case(1)
        z%p = f*x%p + g*xreal_bigi*y%p
        z%i = x%i
      case(-1)
        z%p = f*xreal_bigi*x%p + g*y%p
        z%i = y%i
      case(2:)
        z%p = f*x%p
        z%i = x%i
      case(:-2)
        z%p = g*y%p
        z%i = y%i
    end select
    z = xreal_norm(z)

  end function xreal_fxpgy

  function xadd(x,y) result(z)
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: z

    z = xreal_fxpgy(1.0_dp,x,1.0_dp,y)

  end function xadd

  function xsub(x,y) result(z)
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: z

    z = xreal_fxpgy(1.0_dp,x,-1.0_dp,y)

  end function xsub

  function xreal_base10(x) result(y)
    type(xreal_type), intent(in) :: x
    type(xreal_type) :: y

    integer(kind=i4b), parameter :: i10 = nint(log10(xreal_big))
    real(kind=dp), parameter :: p10 = xreal_big*10.0_dp**(-i10)

    y = x
    if (x%i/=0) then
      y%i = nint(log10(abs(x%p)))
      y%p = (x%p*10.0_dp**(-y%i))*(p10**(x%i))
      y%i = y%i+i10*x%i
    end if

  end function xreal_base10

  subroutine xreal_test()

    real(kind=dp) :: f, g
    integer(kind=i4b) :: i, j
    type(xreal_type) :: x, y, z

    print *, "----- xreal_test() -----"
    f = 3.0d100
    g = 5.0d99
    x = f
    y = g
    print *, "f=", f, " x=", x, " g=", g, " y=", y
    z = x * y
    print *, "x*y=", z, xreal_base10(z)
    z = x + y
    print *, "x+y=", z, xreal_base10(z)
    z = x - y
    print *, "x-y=", z, xreal_base10(z)
    z = x / y
    print *, "x/y=", z, xreal_base10(z)
    z = xreal_fxpgy(f,x,g,y)  
    print *, "fx+gy=", z, xreal_base10(z)
    print *, "--------------------"

  end subroutine xreal_test

end module xreal_module
