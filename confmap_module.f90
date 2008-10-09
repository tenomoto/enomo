module confmap_module
! generates grid and map factors using conformal mapping
  use type_module, only: i4b, dp, dpc
  implicit none
  private

! References
! Bentsen et al. (1999), Mon. Wea. Rev.
! Murray (1996), J. Comp. Phys.
! Couriter and Geleyn (1988) QJRMS

! History
! 2008-10-08 Takeshi Enomoto <eno@jamstec.go.jp>
! * allows the South Pole at either a, b, c
! 2008-09-18 Takeshi Enomoto <eno@jamstec.go.jp>
! * first version

  public :: confmap_stereo, confmap_invstereo, confmap_linear, &
            confmap_grid2earth, confmap_antipode, confmap_eqcolat, &
            confmap_midpoint

contains

  function confmap_stereo(lon, colat) result(z)
    use math_module, only: math_i, math_inf, math_pi
    implicit none

    real(kind=dp), intent(in) :: lon, colat
    complex(kind=dpc) :: z

! 0<=colat<pi
    if (colat < math_pi) then
      z = tan(0.5_dp*colat)*exp(math_i*lon)
    else
      z = math_inf
    end if
     
  end function confmap_stereo

  subroutine confmap_invstereo(z,lon,colat)
	  use math_module, only: math_arg, math_inf, math_pi
    implicit none

    complex(kind=dpc), intent(in) :: z
    real(kind=dp), intent(out) :: lon, colat

    if (z/=math_inf) then
      lon = math_arg(z)
      colat = 2.0_dp*atan(abs(z))
    else
      lon = 0.0_dp
      colat = math_pi
    end if

  end subroutine confmap_invstereo 

  function confmap_linear(z,a,b,c,d) result (w)
    use math_module, only: math_inf
    implicit none

    complex(kind=dpc), intent(in) :: z, a, b, c, d
    complex(kind=dpc) :: w

    complex(kind=dpc) :: denom

    denom = c*z + d
    if (abs(denom)==0.0_dp) then
      w = math_inf
    else
      w = (a*z+b)/denom
    end if
 
  end function confmap_linear

  function confmap_antipode(lona, colata) result(b)
    use math_module, only: pi=>math_pi, pih=>math_pih
    implicit none

    real(kind=dp), intent(in) :: lona, colata
    complex(kind=dpc) :: b

    real(kind=dp) :: lonb, colatb

    if (lonb>=pi) then
      lonb = lona - pi
    else
      lonb = lona + pi
    end if
    colatb = pi-colata
    b = confmap_stereo(lonb,colatb)

  end function confmap_antipode

  function confmap_eqcolat(lona,colata,alpha) result(c)
    use math_module, only: pi=>math_pi
    implicit none

    real(kind=dp), intent(in) :: lona,colata, alpha 
    complex(kind=dpc) :: c

    real(kind=dp) :: lonc, colatc

    lonc = lona
    colatc = colata + 2.0_dp*atan(1.0_dp/alpha)
    if (colatc>pi) then
      if (lonc>=pi) then
        lonc = lonc - pi
      else
        lonc = lonc + pi
      end if
      colatc = pi - colatc
    end if
    c = confmap_stereo(lonc,colatc)
 
  end function confmap_eqcolat

  function confmap_midpoint(lona,colata,lonb,colatb) result(c)
    use math_module, only: math_atan2
    implicit none

    real(kind=dp), intent(in) :: lona, colata, lonb, colatb
    complex(kind=dpc) :: c

    real(kind=dp) :: lonc, colatc, cx, cy, cz

    cx = cos(lona)*sin(colata)+cos(lonb)*sin(colatb)
    cy = sin(lona)*sin(colata)+sin(lonb)*sin(colatb)
    cz = cos(colata)+cos(colatb)

    lonc = math_atan2(cy,cx)
    colatc = math_atan2(sqrt(cx*cx+cy*cy),cz)
    c = confmap_stereo(lonc,colatc)

  end function confmap_midpoint

  subroutine confmap_grid2earth(lon, colat, a, b, c, lone, colate, m)
    use math_module, only: pih=>math_pih, pi=>math_pi, pi2=>math_pi2, math_inf
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, colat
    complex(kind=dp), intent(in) :: a, b, c

    integer(kind=i4b) :: nx, ny, i, j
    real(kind=dp), dimension(size(lon),size(colat)), intent(inout) :: lone, colate, m

    real(kind=dp) :: lonb, colatb
    complex(kind=dp) :: w, z, ai, bi, ci, di, det

    nx = size(lon)
    ny = size(colat)

    if ((a/=math_inf).and.(b/=math_inf).and.(c/=math_inf)) then
      ai =  b*(c-a) ! -d
      bi = -a*(c-b) !  b
      ci =    (c-a) !  c
      di =   -(c-b) ! -a
      det = (c-b)*(c-a)*(a-b) ! ad-bc
    else if (a==math_inf) then
      ai = b
      bi = c-b
      ci = 1.0_dp
      di = 0.0_dp
      det = -(c-b)
    else if (b==math_inf) then
      ai = -(c-a)
      bi = -a
      ci =  0.0_dp
      di = -1.0_dp
      det = c-a
    else if (c==math_inf) then
      ai =  b
      bi = -a
      ci =  1.0_dp
      di = -1.0_dp
      det = -b+a
    else
      print *, "Invalid a, b, or c"
      stop
    end if
    do j=1, ny
      if (colat(j)==pi) then
        call confmap_invstereo(b,lonb,colatb) 
        lone(:,j) = lonb
        colate(:,j) = colatb
        if (b/=math_inf) then
          m(:,j) = (1.0_dp+abs(b)**2)*abs(ci/det)
        else
          m(:,j) = 1.0_dp
        end if
      else
        do i=1, nx
          w = confmap_stereo(lon(i),colat(j))
          z = confmap_linear(w,ai,bi,ci,di)
          if ((a==math_inf).and.(abs(w)==0.0_dp)) then
            m(i,j) = abs(det)
          else if ((c==math_inf).and.(abs(w-1.0_dp)==0.0_dp)) then
            m(i,j) = 0.5_dp*abs(det)
          else
            m(i,j) = (1.0_dp+abs(z)**2)/(1.0_dp+abs(w)**2)*abs((ci*w+di)**2/det)
          end if
          call confmap_invstereo(z,lone(i,j),colate(i,j)) 
          if (lone(i,j)<0) then
            lone(i,j) = lone(i,j) + pi2
          end if
        end do
      end if
    end do

  end subroutine confmap_grid2earth

end module confmap_module
