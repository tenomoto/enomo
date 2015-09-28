module xreal_module
! Extended exponent of floating-point numbers

! Reference:
!   Fukushima, Toshio, 2011: Numerical computation of spherical
!     harmonics of arbitrary degree and order by extending
!     exponent of floating point numbers. J. Geodesy,
!     doi:10.1007//s00190-011-0519-2

  use kind_module, only: dp, qp, i4b
  implicit none
  private

  integer(kind=i4b), parameter :: &
    ind = 960, indh = ind/2, indq = 1600, indqh = indq/2
  real(kind=dp), parameter, public :: &
    xreal_big  = 2.0_dp**ind,  xreal_bigi  = 2.0_dp**(-ind), &
    xreal_bigs = 2.0_dp**indh, xreal_bigsi = 2.0_dp**(-indh)
  real(kind=qp), parameter, public :: &
    xreal_bigq  = 2.0_qp**indq,  xreal_bigqi  = 2.0_qp**(-indq), &
    xreal_bigqs = 2.0_qp**indqh, xreal_bigqsi = 2.0_qp**(-indqh)

  type xreal_type
    real(kind=dp) :: p     ! principal component
    integer(kind=i4b) :: i ! auxiliary index
  end type xreal_type

  type xreal_quad_type
    real(kind=qp) :: p     ! principal component
    integer(kind=i4b) :: i ! auxiliary index
  end type xreal_quad_type

  interface assignment(=)
    module procedure x_assign_f, f_assign_x, xq_assign_fq, fq_assign_xq
  end interface

  interface operator(+)
    module procedure xadd
  end interface

  interface operator(-)
    module procedure xsub
  end interface

  interface operator(*)
    module procedure xmul, fx, xmulq, fxq
  end interface

  interface operator(/)
    module procedure xdiv
  end interface

  interface operator(**)
    module procedure xpowi
  end interface

  interface operator(==)
    module procedure xeq
  end interface

  interface operator(/=)
    module procedure xne
  end interface

  interface operator(>)
    module procedure xgt
  end interface

  interface operator(>=)
    module procedure xge
  end interface

  interface operator(<)
    module procedure xlt
  end interface

  interface operator(<=)
    module procedure xle
  end interface

  public :: xreal_type, xreal_quad_type, xreal_norm, xreal_normq, xreal_fxpgy, xreal_fxpgyq, xreal_base10, xreal_test, &
    assignment(=), operator(+), operator(-), operator(*), operator(/), &
    operator(**), operator(==), operator(/=), operator(>), operator(>=), &
    operator(<), operator(<=)
  private :: x_assign_f, f_assign_x, xq_assign_fq, fq_assign_xq, xmul, xmulq, xdiv, xadd, xsub, &
    xpowi, xeq, xne, xgt, xge, xlt, xle, fx, fxr, fxq

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

  subroutine xq_assign_fq(x,f)
    type(xreal_quad_type), intent(out) :: x
    real(kind=qp), intent(in) :: f

    x%p = f
    x%i = 0

  end subroutine xq_assign_fq

  subroutine fq_assign_xq(g,x)
    real(kind=qp), intent(out) :: g
    type(xreal_quad_type), intent(in) :: x

    select case(x%i)
      case(0)
        g = x%p
      case(:-1) ! underflow
        g = x%p*xreal_bigqi
      case(1:)  ! overflow
        g = x%p*xreal_bigq
    end select

  end subroutine fq_assign_xq

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

  function xreal_normq(x) result(y)
    type(xreal_quad_type), intent(in) :: x
    type(xreal_quad_type) :: y

    real(kind=qp) :: w

    w = abs(x%p)
    if (w>=xreal_bigqs) then
      y%p = x%p*xreal_bigqi
      y%i = x%i + 1
    else if (w<xreal_bigqsi) then
      y%p = x%p*xreal_bigq
      y%i = x%i - 1
    else
      y = x
    end if

  end function xreal_normq

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


  function xmulq(x,y) result(z)
    type(xreal_quad_type), intent(in) :: x, y
    type(xreal_quad_type) :: s, z

    z = xreal_normq(x)
    s = xreal_normq(y)
    z%p = z%p * s%p
    z%i = z%i + s%i
    z = xreal_normq(z)

  end function xmulq

  function fxq(f,x) result(z)
    real(kind=qp), intent(in) :: f
    type(xreal_quad_type), intent(in) :: x
    type(xreal_quad_type) :: z

    z = xreal_normq(x)
    z%p = f * z%p
    z%i = z%i
    z = xreal_normq(z)

  end function fxq

  function xdiv(x,y) result(z)
    type(xreal_type), intent(in) :: x, y
    type(xreal_type) :: s, z

    z = xreal_norm(x)
    s = xreal_norm(y)
    z%p = z%p / s%p
    z%i = z%i - s%i
    z = xreal_norm(z)

  end function xdiv

  function fxr(f,x) result(z)
    real(kind=dp), intent(in) :: f
    type(xreal_type), intent(in) :: x
    type(xreal_type) :: z

    z = xreal_norm(x)
    z%p = f / z%p
    z%i = z%i
    z = xreal_norm(z)

  end function fxr

  function xpowi(x,n) result(y)
    type(xreal_type), intent(in) :: x
    integer(kind=i4b), intent(in) :: n

    integer :: m
    type(xreal_type) :: y, z

    y = 1.0_dp;
    z = xreal_norm(x)
    m = n
    if (n < 0) then
      z = fxr(1.0_dp, z)
      m = -n
    else
    end if
    do
      if (iand(m, 1)==1) then
        y = y * z
      end if
      m = ishft(m, -1)
      z = z * z
      if (m <= 0) then
        exit
      end if
    end do
    y = xreal_norm(y)

  end function xpowi

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

  function xreal_fxpgyq(f,x,g,y) result(z)
    real(kind=qp), intent(in) :: f, g
    type(xreal_quad_type), intent(in) :: x, y
    type(xreal_quad_type) :: z

    integer(kind=i4b) :: id

    id = x%i - y%i
    select case(id)
      case(0)
        z%p = f*x%p + g*y%p
        z%i = x%i
      case(1)
        z%p = f*x%p + g*xreal_bigqi*y%p
        z%i = x%i
      case(-1)
        z%p = f*xreal_bigqi*x%p + g*y%p
        z%i = y%i
      case(2:)
        z%p = f*x%p
        z%i = x%i
      case(:-2)
        z%p = g*y%p
        z%i = y%i
    end select
    z = xreal_normq(z)

  end function xreal_fxpgyq

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

  function xeq(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    l = (xx%p == yy%p) .and. (xx%i == yy%i)
    
  end function xeq

  function xne(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    l = (xx%p /= yy%p) .or. (xx%i /= yy%i)
    
  end function xne

  function xgt(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    if (xx%i==yy%i) then
      l = (xx%p > yy%p)
    else
      l = xx%i > yy %i
    end if
    if (xx%p*yy%p<0) then
      l = .not. l
    end if
    
  end function xgt

  function xge(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    if (xx%i==yy%i) then
      l = (xx%p >= yy%p)
    else
      l = xx%i >= yy %i
    end if
    if (xx%p*yy%p<0) then
      l = .not. l
    end if
    
  end function xge

  function xlt(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    if (xx%i==yy%i) then
      l = (xx%p < yy%p)
    else
      l = xx%i < yy %i
    end if
    if (xx%p*yy%p<0) then
      l = .not. l
    end if
    
  end function xlt

  function xle(x,y) result(l)
    type(xreal_type), intent(in) :: x, y
    logical :: l 
    type(xreal_type) :: xx, yy
    
    xx = xreal_norm(x)
    yy = xreal_norm(y)
    if (xx%i==yy%i) then
      l = (xx%p <= yy%p)
    else
      l = xx%i <= yy %i
    end if
    if (xx%p*yy%p<0) then
      l = .not. l
    end if
    
  end function xle

  function xreal_base10(x) result(y)
    type(xreal_type), intent(in) :: x
    type(xreal_type) :: y

    integer(kind=i4b) :: i10
    real(kind=dp) :: p10

    i10 = nint(log10(xreal_big))
    p10 = xreal_big*10.0_dp**(-i10)

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
    z = x * x * x * x * x
    print *, "x*x*x*x*x=", z, xreal_base10(z)
    z = x**5
    print *, "x**5=", z, xreal_base10(z)
    z = x**(-2)
    print *, "x**-2=", z, xreal_base10(z)
    print *, "x**3==y:", (x**3)==y
    print *, "x**3/=y:", (x**3)/=y
    print *, "x**3>y:", (x**3)>y
    print *, "x**3>=y:", (x**3)>=y
    print *, "x**3>=x**3:", (x**3)>=(x**3)
    print *, "x**3<y:", (x**3)<y
    print *, "x**3<=x**3:", (x**3)<=(x**3)

  end subroutine xreal_test

end module xreal_module
