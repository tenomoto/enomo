module interpolate_module

! interpolate in a stencil

! References:
!   - Ritchie (1987) describes selection of points across the pole
!   - Interpolation schemes from Numerical Recipes
! Author: T. Enomoto
! History: 
! 2004-09-10 some simplification
! 2004-03    first version

  use type_module, only: i4b, dp
  implicit none
  private

  private :: cubiclagrange, quinticlagrange
  public :: interpolate_bilinear, interpolate_bicubic, interpolate_linpol, interpolate_spcher

contains

  function interpolate_bilinear(f, t, u) result(fi)
    implicit none

    real(kind=dp), dimension(4), intent(in) :: f
    real(kind=dp), intent(in) :: t, u
    real(kind=dp) :: fi

    real(kind=dp) :: t1, u1

    t1 = 1.0_dp-t
    u1 = 1.0_dp-u
    fi = u1*(t1*f(1)+t*f(2)) + u*(t*f(3)+t1*f(4))

  end function interpolate_bilinear

  function interpolate_bicubic(f, fx, fy, fxy, dlon, dlat, t, u) result(fi)
    implicit none

    real(kind=dp), dimension(4), intent(in) :: f, fx, fy, fxy
    real(kind=dp), intent(in) :: dlon, dlat, t, u
    real(kind=dp) :: fi

    integer(kind=i4b) :: i, j
    real(kind=dp), dimension(16) :: c
    real(kind=dp), dimension(16) :: x
    real(kind=dp) :: dxdy, x12, x14, x58, x67, x910, &
      x1112, x1316, x1314, x1516, x1234, x5678, x9101112

! 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
! 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
!-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
! 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
! 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
! 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
! 0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
! 0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
!-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
! 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
! 9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
!-6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
! 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
! 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
!-6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
! 4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1

    dxdy = dlon*dlat
    do i=1, 4
      x(i) = f(i)
      x(i+4) = fx(i)*dlon
      x(i+8) = fy(i)*dlat
      x(i+12) = fxy(i)*dxdy
    end do
    x12 = x(1)-x(2)
    x14 = x(1)-x(4)
    x58 = x(5)-x(8)
    x67 = x(6)-x(7)
    x910 = x(9)-x(10)
    x1112 = x(11)-x(12)
    x1314 = x(13)+x(14)
    x1316 = x(13)+x(16)
    x1516 = x(15)+x(16)
    x1234 = x12+x(3)-x(4)
    x5678 = x58+x67
    x9101112 = x910-x1112
    c( 1) = x( 1)
    c( 2) = x( 9)
    c( 3) = -3.d0*x14-2.d0*x( 9)-x(12)
    c( 4) =  2.d0*x14+     x( 9)+x(12)
    c( 5) = x( 5)
    c( 6) = x(13)
    c( 7) = -3.d0*x58-2.d0*x(13)-x(16)
    c( 8) =  2.d0*x58+x1316
    c( 9) = -3.d0*x12 -2.d0*x( 5)-x( 6)
    c(10) = -3.d0*x910-2.d0*x(13)-x(14)
    c(11) =  9.d0*x1234+6.d0*(x58+x910)+4.d0*x(13) &
              +3.d0*(x67-x1112)+2.d0*(x(14)+x(16))+x(15)
    c(12) = -6.d0*x1234-4.d0*x58-3.d0*x9101112 &
              -2.d0*(x67+x1316)-x(14)-x(15)
    c(13) =  2.d0*x12 +x( 5)+x( 6)
    c(14) =  2.d0*x910+x1314
    c(15) = -6.d0*x1234-4.d0*x910-3.d0*x5678 &
              +2.d0*(x1112-x1314)-x1516
    c(16) =  4.d0*x1234+2.d0*(x5678+x9101112)+x1314+x1516

    fi = 0.0d0
    do j=4,1,-1
      i = 4*(j-1)
      fi = t*fi+((c(i+4)*u+c(i+3))*u+c(i+2))*u+c(i+1)
    end do

  end function interpolate_bicubic

  function interpolate_linpol(f,lon,lat,x,y,t) result(fi)
    implicit none

    integer(kind=i4b), parameter :: nx = 4, ny = 4
    real(kind=dp), dimension(nx,ny), intent(in) :: f
    real(kind=dp), dimension(nx), intent(in) :: lon
    real(kind=dp), dimension(ny), intent(in) :: lat
    real(kind=dp), intent(in) :: x,y,t
    real(kind=dp) :: fi

! full cubic lagrange interpolation if t<0
    real(kind=dp), dimension(ny) :: ytmp
    integer(kind=i4b) :: j
    if (t>=0) then
      ytmp(1) = (1.0_dp-t)*f(2,1)+t*f(3,1)
      ytmp(4) = (1.0_dp-t)*f(2,4)+t*f(3,4)
    else
      ytmp(1) = cubiclagrange(lon,f(1:nx,1),x)
      ytmp(4) = cubiclagrange(lon,f(1:nx,4),x)
    end if
    ytmp(2) = cubiclagrange(lon,f(1:nx,2),x)
    ytmp(3) = cubiclagrange(lon,f(1:nx,3),x)
    fi = cubiclagrange(lat,ytmp,y)

  end function interpolate_linpol

  function interpolate_spcher(f,fx,lat,dlon,y,t) result(fi)
    implicit none

    integer(kind=i4b), parameter :: nx = 2
    real(kind=dp), dimension(:), intent(in) :: lat
    real(kind=dp), dimension(nx,size(lat)), intent(in) :: f, fx
    real(kind=dp), intent(in) :: dlon, y, t
    real(kind=dp) :: fi

    real(kind=dp), dimension(size(lat)) :: ytmp
    real(kind=dp) :: tt, ttt, c1, c2, c3, c4
    integer(kind=i4b) :: j, ny

! Cubic Hermie spline interpolation
! Source: Wikipedia (Cubic Hermite spline), www.cubic.org/docs/hermite.htm
    tt  = t*t
    ttt = t*tt
    ny = size(lat)
    do j=1, ny
!   coeff | 2, -3,  0,  1|  & ! hermite
!         |-2,  3,  0,  0|  &
!         | 1, -2,  1,  0|  &
!         | 1, -1,  0,  0|
      c1 =  2.d0*ttt-3.d0*tt  +1.0
      c2 = -2.d0*ttt+3.d0*tt
      c3 =       ttt-2.d0*tt+t
      c4 =       ttt-     tt
      ytmp(j) = c1*f(1,j)+c2*f(2,j)+(c3*fx(1,j)+c4*fx(2,j))*dlon
    end do
 
    if (ny==6) then
      fi = quinticlagrange(lat,ytmp,y)
    else if (ny==4) then
      fi = cubiclagrange(lat,ytmp,y)
    else
      print *, "*** error in interpolate_spcher ny=", ny
      stop
    end if

  end function interpolate_spcher

! private routines

  function cubiclagrange(xa, ya, x) result(y)
    implicit none

    real(kind=dp), dimension(4), intent(in) :: xa, ya
    real(kind=dp), intent(in) :: x 
    real(kind=dp) :: y

    real(kind=dp) :: y12, y23, y34, y123, y234

    y12  = ((x-xa(2))*ya(1)+(xa(1)-x)*ya(2))/(xa(1)-xa(2))
    y23  = ((x-xa(3))*ya(2)+(xa(2)-x)*ya(3))/(xa(2)-xa(3))
    y34  = ((x-xa(4))*ya(3)+(xa(3)-x)*ya(4))/(xa(3)-xa(4))

    y123 = ((x-xa(3))*y12+(xa(1)-x)*y23)/(xa(1)-xa(3))
    y234 = ((x-xa(4))*y23+(xa(2)-x)*y34)/(xa(2)-xa(4))

    y    = ((x-xa(4))*y123+(xa(1)-x)*y234)/(xa(1)-xa(4))

  end function cubiclagrange
      
  function quinticlagrange(xa, ya, x) result(y)
    implicit none

    real(kind=dp), dimension(6), intent(in) :: xa, ya
    real(kind=dp), intent(in) :: x 
    real(kind=dp) :: y

    real(kind=dp) :: y12, y23, y34, y45, y56, &
      y123, y234, y345, y456, y1234, y2345, y3456, y12345, y23456

    y12    = ((x-xa(2))*ya(1)+(xa(1)-x)*ya(2))/(xa(1)-xa(2))
    y23    = ((x-xa(3))*ya(2)+(xa(2)-x)*ya(3))/(xa(2)-xa(3))
    y34    = ((x-xa(4))*ya(3)+(xa(3)-x)*ya(4))/(xa(3)-xa(4))
    y45    = ((x-xa(5))*ya(4)+(xa(4)-x)*ya(5))/(xa(4)-xa(5))
    y56    = ((x-xa(6))*ya(5)+(xa(5)-x)*ya(6))/(xa(5)-xa(6))
 
    y123   = ((x-xa(3))*y12+(xa(1)-x)*y23)/(xa(1)-xa(3))
    y234   = ((x-xa(4))*y23+(xa(2)-x)*y34)/(xa(2)-xa(4))
    y345   = ((x-xa(5))*y34+(xa(3)-x)*y45)/(xa(3)-xa(5))
    y456   = ((x-xa(6))*y45+(xa(4)-x)*y56)/(xa(4)-xa(6))

    y1234  = ((x-xa(4))*y123+(xa(1)-x)*y234)/(xa(1)-xa(4))
    y2345  = ((x-xa(5))*y234+(xa(2)-x)*y345)/(xa(2)-xa(5))
    y3456  = ((x-xa(6))*y345+(xa(3)-x)*y456)/(xa(3)-xa(6))

    y12345 = ((x-xa(5))*y1234+(xa(1)-x)*y2345)/(xa(1)-xa(5))
    y23456 = ((x-xa(6))*y2345+(xa(2)-x)*y3456)/(xa(2)-xa(6))

    y      = ((x-xa(6))*y12345+(xa(1)-x)*y23456)/(xa(1)-xa(6))

  end function quinticlagrange

end module interpolate_module
