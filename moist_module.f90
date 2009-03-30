module moist_module
! collection of functions to calculate moist variables
  use kind_module, only: i4b, dp
  use math_module, only: undef=>math_undef
  use air_module, only: Rd=>air_rd, pr=>air_pr
  use water_module, only: T0=>water_t0, eps=>water_eps
  implicit none
  private

! Source: Bolton (1980)

! Author: Takeshi Enomoto <eno@jamstec.go.jp>
! History:
!   2007-11-23 first version

  logical, private :: debug = .false.

  public :: moist_es, moist_es2T, moist_q2e, moist_q2r, moist_e2r, moist_e2q, &
            moist_RH, moist_Tv, moist_Td, moist_Tl, moist_pl, &
            moist_theta, moist_theta2T, moist_thetae, moist_Tp, moist_cape
  private :: kappa, fthetae

contains

  function moist_es(T) result(es)
    implicit none

    real(kind=dp), intent(in) :: T ! K 
    real(kind=dp) :: es ! Pa

    real(kind=dp) :: Tc

    Tc = T - T0
    es = 100.0_dp*exp(19.482_dp - 4303.4_dp/(Tc+243.5_dp))

  end function moist_es

  function moist_es2T(es) result(T)
    implicit none

    real(kind=dp), intent(in) :: es ! Pa
    real(kind=dp) :: T ! K

    T = 4303.4_dp/(19.482_dp-log(0.01_dp*es)) - 243.5_dp + T0 

  end function moist_es2T

  function moist_q2e(q,p) result(e)
    implicit none

    real(kind=dp), intent(in) :: q, p
    real(kind=dp) ::e

    e = p*q/(eps+(1.0_dp-eps)*q)

  end function moist_q2e

  function moist_q2r(q) result(r)
    implicit none

    real(kind=dp), intent(in) :: q
    real(kind=dp) :: r

    r = q / (1.0_dp - q)

  end function moist_q2r

  function moist_e2r(e,p) result(r)
    implicit none

    real(kind=dp), intent(in) :: e, p
    real(kind=dp) :: r

    r = eps*e/(p-e)

  end function moist_e2r

  function moist_e2q(e,p) result(q)
    implicit none

    real(kind=dp), intent(in) :: e, p
    real(kind=dp) :: q

    q = eps*e/(p-(1.0_dp-eps)*e)

  end function moist_e2q

  function moist_RH(T,q,p) result(rh)
    implicit none

    real(kind=dp), intent(in) :: T, q, p
    real(kind=dp) :: rh

    rh = moist_q2e(q,p)/moist_es(T)

  end function moist_RH

  function moist_Tv(T,q) result(Tv)
    implicit none

    real(kind=dp), intent(in) :: T, q
    real(kind=dp) :: Tv

    real(kind=dp), parameter :: epsd = (1.-eps)/eps

    Tv = T*(1.+epsd*q)

  end function moist_Tv

  function moist_Td(q,p) result(Td)
    implicit none

    real(kind=dp), intent(in) :: q, p
    real(kind=dp) :: Td

    real(kind=dp) :: e

    e = moist_q2e(q,p)
    Td = moist_es2T(e)

  end function moist_Td

  function moist_Tl(T,q,p) result(Tl)
! Bolton 1980
    implicit none

    real(kind=dp), intent(in) :: T, q, p
    real(kind=dp) :: Tl

    real(kind=dp) :: e

    e = moist_q2e(q,p)
! NB. log(e=0) = -Inf
    Tl = 2840.0_dp/(3.5_dp*log(T)-log(e*0.01_dp)-4.805_dp) + 55.0_dp

  end function moist_Tl

  function moist_pl(theta,q,Tl) result(pl)
    implicit none

    real(kind=dp), intent(in) :: theta, q, Tl
    real(kind=dp) :: pl

    real(kind=dp) :: r

    r = moist_q2r(q)
    pl = pr*(Tl/theta)**(1./kappa(r))

  end function moist_pl

  function moist_thetae(T,q,p,Tl) result(thetae)
    implicit none

    real(kind=dp), intent(in) :: T, q, p, Tl
    real(kind=dp) :: thetae

    real(kind=dp) :: r

    r = moist_q2r(q)
    thetae = T * ((pr/p)**kappa(r)) * exp(fthetae(Tl,r))

  end function moist_thetae

  function moist_Tp(thetae,Tl,T,p) result(Tp)
    implicit none

    real(kind=dp), parameter :: small = 1.e-3
    integer(kind=i4b), parameter :: n = 2
    real(kind=dp), intent(in) :: thetae, Tl, T, p
    real(kind=dp) :: Tp

    integer(kind=i4b) :: i
    real(kind=dp) :: te, dte, rs, Tc, s

    Tp = T ! environmental temperature as initial guess
    do i=1, n
      rs = moist_e2r(moist_es(Tp),p)
      te = Tp * ((pr/p)**kappa(rs)) * exp(fthetae(Tl,rs))
      Tc = Tp - T0
! s = dthetae/dT
! Assumptions:
! - kappa is constant
! - thetae formula by Betts and Dugan 1973
! - r = eps e/p
      s = thetae*(1./Tp+(4304.4*2675.*rs)/(Tl*(Tc+243.5)**2))
      dte =  te-thetae
      if (abs(dte)<small) then
        exit
      end if
      Tp = Tp - (te-thetae)/s
    end do

  end function moist_Tp

  subroutine moist_cape(Tsrc,qsrc,psrc, T,q,p, Tp,cin,cape,pfc,pnb, &
                        use_virtual,verbose)
! input: T, q, p at source level
! output: parcel T, CIN, CAPE, pressure at LFC anc LNB
    implicit none

    real(kind=dp), intent(in) :: Tsrc, qsrc, psrc
    real(kind=dp), dimension(:), intent(in) :: T, q, p
    real(kind=dp), dimension(:), intent(inout) :: Tp
    real(kind=dp), intent(out) :: cin, cape, pfc, pnb
    logical, optional :: use_virtual, verbose

    real(kind=dp) :: pmin = 100.d2

    integer(kind=i4b) :: nz, k
    real(kind=dp) :: theta, thetae, Tl, pl, dlnp, b, bb, qp, Tvp, Tve
    logical :: lbuoy, lTv

    nz = size(p)
    lbuoy = .false.
    lTv = present(use_virtual).and.use_virtual
    debug = present(verbose).and.verbose

    if (debug) then
      print *, "Tsrc(C) qsrc(g/kg) psrc(hPa)"
      print *, Tsrc-273.15_dp, qsrc*1000.0_dp, psrc*0.01_dp
      print *, "e=", moist_q2e(qsrc,psrc)*0.01_dp, "hPa"
    end if

    Tl = moist_Tl(Tsrc,qsrc,psrc) ! T at LCL
    theta = moist_theta(Tsrc,qsrc,psrc)
    pl = moist_pl(theta,qsrc,Tl)
    thetae = moist_thetae(Tsrc,qsrc,psrc,Tl)

    if (debug) then
      print *, "Tl(C) pl(hPa) q(g/kg) theta thetae"
      print *, Tl-273.15_dp, pl*0.01_dp, q(1)*1000, theta, thetae
      print *, "p Tp T T-Tp lbuoy dlnp CIN CAPE"
    end if

    cape = 0.0_dp
    cin = 0.0_dp
    pfc = undef
    pnb = undef
    do k=1, nz
      if (p(k)<pmin) then
        exit
      end if
      if (p(k)>pl) then ! dry adiabat below LCL; theta conserved
        Tp(k) = moist_theta2T(theta,qsrc,p(k))
        qp = qsrc
      else              ! moist adiabat above LCL; thetae conserved
        Tp(k) = moist_Tp(thetae,Tl,T(k),p(k))
        qp = moist_e2q(moist_es(Tp(k)),p(k))
      end if
      dlnp = 0.5_dp*(log(p(max(1,k-1)))-log(p(min(nz,k+1))))
      if (lTv) then
        Tvp = moist_Tv(Tp(k),qp)
        Tve = moist_Tv(T(k),q(k))
        b = Tvp - Tve
      else
        b = Tp(k) - T(k)
      end if
      if (b>0.0_dp) then
        if (.not.lbuoy) then
          lbuoy = .true.
          if (k/=1) then
            pfc = p(k-1)*((p(k)/p(k-1))**(bb/(bb-b)))
          else
            pfc = psrc
          end if
        end if
        cape = cape + Rd*b*dlnp
      else
        if (lbuoy) then ! T => F
! Initially lbuoy = F thus k > 1
          pnb = p(k-1)*((p(k)/p(k-1))**(bb/(bb-b)))
          exit
        else
          cin = cin + Rd*b*dlnp
        end if
      end if
      bb = b

      if (debug) then
        print *, 0.01_dp*p(k), Tp(k)-273.15_dp, T(k)-273.15_dp, b, lbuoy, cin, cape
      end if

    end do
    if (cape==0.) then
      cin = 0.
    end if

    if (debug) then
      print *, "pfc=", 0.01_dp*pfc, " hPa pnb=", 0.01_dp*pnb, " hPa", " cape=", cape, " cin=", cin
    end if

  end subroutine moist_cape

  function moist_theta(T,q,p) result(theta)
    implicit none

    real(kind=dp), intent(in) :: T, q, p
    real(kind=dp) :: theta, r

    r = moist_q2r(q)
    theta = T*(pr/p)**kappa(r)

  end function moist_theta

  function moist_theta2T(theta,q,p) result(T)
    implicit none

    real(kind=dp), intent(in) :: theta, q, p
    real(kind=dp) :: T

    real(kind=dp) :: r

    r = moist_q2r(q)
    T = theta*(p/pr)**kappa(r)

  end function moist_theta2T

  function kappa(r) result(kp)
! Bolton 1980
    implicit none
 
    real(kind=dp), intent(in) :: r
    real(kind=dp) :: kp

    kp = 0.2854_dp*(1.0_dp-0.28_dp*r)

  end function kappa

  function fthetae(Tl,r) result(ft)
    implicit none

    real(kind=dp), intent(in) :: Tl, r
    real(kind=dp) :: ft

! Betts and Dugan 1973
!    ft = 2675.*r/Tl

! Bolton 1980
    ft = (3376.0_dp/Tl-2.54_dp)*r*(1.0_dp+0.81_dp*r)

  end function fthetae

end module moist_module
