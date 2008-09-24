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
! 2008-9-18 Takeshi Enomoto <eno@jamstec.go.jp>
! * first version

  public :: confmap_stereo, confmap_invstereo, confmap_linear, &
            confmap_grid2earth, confmap_antipode, confmap_eqcolat, &
            confmap_midpoint

contains

  function confmap_stereo(lon, colat) result(z)
    use math_module, only: math_i
    implicit none

    real(kind=dp), intent(in) :: lon, colat
    complex(kind=dpc) :: z

! 0<=colat<pi
    z = tan(0.5_dp*colat)*exp(math_i*lon)
     
  end function confmap_stereo

  subroutine confmap_invstereo(z,lon,colat)
	  use math_module, only: math_arg
    implicit none

    complex(kind=dpc), intent(in) :: z
    real(kind=dp), intent(out) :: lon, colat

    lon = math_arg(z)
    colat = 2.0_dp*atan(abs(z))

  end subroutine confmap_invstereo 

  function confmap_linear(z,a,b,c,d) result (w)
    implicit none

    complex(kind=dpc), intent(in) :: z, a, b, c, d
    complex(kind=dpc) :: w

    w = (a*z+b)/(c*z+d)
 
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

  subroutine confmap_grid2earth(lon, colat, a, b, c, lone, colate)
    use math_module, only: pih=>math_pih, pi=>math_pi, pi2=>math_pi2
    implicit none

    real(kind=dp), dimension(:), intent(in) :: lon, colat
    complex(kind=dp), intent(in) :: a, b, c

    integer(kind=i4b) :: nx, ny, i, j
    real(kind=dp), dimension(size(lon),size(colat)), intent(inout) :: lone, colate

    real(kind=dp) :: lonb, colatb
    complex(kind=dp) :: w, z, ai, bi, ci, di

    nx = size(lon)
    ny = size(colat)

    ai = -b*(c-a)
    bi =  a*(c-b)
    ci =   -(c-a)
    di =     c-b
    do j=1, ny
      if (colat(j)==pi) then
        call confmap_invstereo(b,lonb,colatb) 
        lone(:,j) = lonb
        colate(:,j) = colatb
      else
        do i=1, nx
          w = confmap_stereo(lon(i),colat(j))
          z = confmap_linear(w,ai,bi,ci,di)
          call confmap_invstereo(z,lone(i,j),colate(i,j)) 
          if (lone(i,j)<0) then
            lone(i,j) = lone(i,j) + pi2
          end if
        end do
      end if
    end do

  end subroutine confmap_grid2earth

end module confmap_module
