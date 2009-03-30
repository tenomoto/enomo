module stability_module
! calculates stability parameters
  use kind_module, only : i4b, dp
  use math_module, only : undef=>math_undef
  implicit none
  private

! Author: Takeshi Enomoto <eno@jamstec.go.jp>
! History:
!   2008-02-21 added CAPE
!   2007-11-23 first version

  logical, parameter :: lTv = .true., & ! use virtual temperature
                        ldbg = .false.

  public :: stability_kindex, stability_cape, stability_ssi

contains

  subroutine stability_kindex(T850, T700, T500, q850, q700, Ki, TT)
    use water_module, only : T0=>water_t0
    use moist_module, only : moist_Td
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: T850, T700, T500, q850, q700
    real(kind=dp), dimension(:,:), intent(inout) :: Ki, TT

    integer(kind=i4b) :: i, j, nx, ny
    real(kind=dp) :: Td850, Td700, e

    nx = size(T850,1)
    ny = size(T850,2)

    Ki(:,:) = 0.
    do j=1, ny
      do i=1, nx
        if ((T850(i,j)/=undef).and.(q850(i,j)/=undef).and. &
            (T700(i,j)/=undef).and.(q700(i,j)/=undef).and. &
            (T500(i,j)/=undef)) then
          Td850 = moist_Td(q850(i,j),850.d2)
          Td700 = moist_Td(q700(i,j),700.d2)
          Ki(i,j) = T850(i,j) - T500(i,j) + Td850 - (T700(i,j)-Td700) - T0
          TT(i,j) = T850(i,j) - T500(i,j) + Td850 - T500(i,j)
        else
          Ki(i,j) = undef
          TT(i,j) = undef
        end if
       end do
    end do
        
  end subroutine stability_kindex

  subroutine stability_cape(T, q, ps, ks, k500, p, ci, ca, pfc, pnb, li)
    use earth_module, only: g => earth_gravity
    use air_module, only: Rd=>air_rd
    use moist_module, only: moist_Tv, moist_Tl, moist_pl, &
        moist_theta, moist_theta2T, moist_thetae, moist_Tp, moist_q2e, moist_cape
    implicit none

! Assumptions
! * temperature decrease with height at dry adiabatic lapse rate
! * moisture interpolated with dlnp
! * parcel at 500 m above surface is lifted

    real(kind=dp), dimension(:,:,:), intent(in) :: T, q
    real(kind=dp), dimension(:,:), intent(in) :: ps
    integer(kind=i4b), dimension(:,:), intent(in) :: ks
    integer(kind=i4b), intent(in) :: k500
    real(kind=dp), dimension(:), intent(in) :: p
    real(kind=dp), dimension(:,:), intent(inout) :: ci, ca, pfc, pnb, li

!    real(kind=dp), parameter :: gamma = 0.0065

    real(kind=dp) :: T500m, q500m, p500m, rdgg, grd, gamma500, gamma
    real(kind=dp), dimension(size(p)) :: Tp
    integer(kind=i4b) :: i, j, k, k0, k1, nx, ny, nz

    nx = size(T,1)
    ny = size(T,2)
    nz = size(p)
    grd = g/Rd
!    rdgg = Rd/g*gamma
!    gamma500 = gamma*500.

    do j=1, ny
      do i=1, nx
        k0 = ks(i,j) 
        k0 = min(nz,k0)
        k1 = min(nz,k0+1)

        gamma = grd*(log(T(i,j,k1))-log(T(i,j,k0)))/(log(p(k1))-log(p(k0)))
        rdgg = gamma/grd
        gamma500 = gamma*500.

        T500m = T(i,j,k0) * (1.+rdgg*(log(ps(i,j))-log(p(k0)))) - gamma500
        p500m = ps(i,j) * (1.-grd/T(i,j,k0)*500.)
        if (k0/=k1) then
          q500m = q(i,j,k0) &
            + (q(i,j,k1)-q(i,j,k0))/(log(p(k1))-log(p(k0)))*(log(p500m)-log(p(k0)))
        else
          q500m = q(i,j,k0)
        end if

! debug
if ((i==160).and.(j==120)) then
!print *, 0.01*p(k0), T(i,j,k0)-273.15, 0.01*ps(i,j)
        call moist_cape(T500m,q500m,p500m, T(i,j,k1:nz),q(i,j,k1:nz),p(k1:nz), &
                        Tp(k1:nz),ci(i,j),ca(i,j),pfc(i,j),pnb(i,j), &
                        use_virtual=lTv, verbose=ldbg)
!stop
else
        call moist_cape(T500m,q500m,p500m, T(i,j,k1:nz),q(i,j,k1:nz),p(k1:nz), &
                        Tp(k1:nz),ci(i,j),ca(i,j), pfc(i,j), pnb(i,j), &
                        use_virtual=lTv)
end if

! calculate lifted index

        if (T(i,j,k500)/=undef) then
          li(i,j) = T(i,j,k500) - Tp(k500)
        else
          li(i,j) = undef
        end if

      end do
    end do

  end subroutine stability_cape

  subroutine stability_ssi(T,q,k850,k500,ssi)
    use moist_module, only : moist_Tl, moist_pl, moist_theta, &
      moist_thetae, moist_Tp, moist_theta2T
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: T, q
    integer(kind=i4b), intent(in) :: k850, k500
    real(kind=dp), dimension(:,:), intent(out) :: ssi

    real(kind=dp), parameter :: psrc = 850.d2
    integer(kind=i4b) :: i, j, k, nx, ny, nz
    real(kind=dp) :: Tl, pl, theta, thetae, Tp, Tsrc, qsrc

    nx = size(T,1)
    ny = size(T,2)
    nz = size(T,3)

    do j=1, ny
      do i= 1, nx
        Tsrc = T(i,j,k850)
        qsrc = q(i,j,k850)
        if ((Tsrc/=undef).and.(qsrc/=undef)) then
          Tl = moist_Tl(Tsrc,qsrc,psrc)
          theta = moist_theta(Tsrc,qsrc,psrc)
          pl = moist_pl(theta,qsrc,psrc)
          thetae = moist_thetae(Tsrc,qsrc,psrc,Tl)
          if (pl>500.d2) then
            Tp = moist_Tp(thetae,Tl,T(i,j,k500),500.d2)
          else
            Tp = moist_theta2T(theta,qsrc,500.d2)
          end if
          ssi(i,j) = T(i,j,k500) - Tp
        else
          ssi(i,j) = undef
        end if
      end do
    end do

  end subroutine stability_ssi
      
end module stability_module

