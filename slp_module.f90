module slp_module

  use kind_module, only : dp, i4b
  implicit none
  private

  real(kind=dp), dimension(:,:), allocatable, private :: t0, ts, gamma, x
  logical, private :: linit = .false.

  public :: slp_init, slp_clean, slp_calc

contains

  subroutine slp_init(ps)
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: ps

    integer(kind=i4b) :: nx, ny

    nx = size(ps,1)
    ny = size(ps,2)
    allocate(t0(nx,ny),ts(nx,ny),gamma(nx,ny),x(nx,ny))
    linit = .true.

  end subroutine slp_init

  subroutine slp_clean()
    implicit none

    if (linit) then
      deallocate(t0,ts,gamma,x)
    end if

  end subroutine slp_clean

  subroutine slp_calc(ps, tl, zs, sl, slp)
    use earth_module, only: g=>earth_gravity, a=>earth_radius
    use air_module, only: ra=>air_rd
    implicit none

    real(kind=dp), dimension(:,:), intent(in) :: &
      ps, & ! surface pressure Pa
      tl, & ! lowest level temperature
      zs    ! surface height
    real(kind=dp), intent(in) :: sl ! lowest level sigma
    real(kind=dp), dimension(:,:), intent(inout) :: slp

    real(kind=dp), parameter :: gamma_st = 0.0065_dp
    real(kind=dp), parameter :: &
      zsmin = 0.001_dp/g, tmin = 255.0_dp, tmax = 290.5_dp

    if (.not.linit) then
      call slp_init(ps)
    end if

    gamma = gamma_st

    ts = (1.0_dp + gamma_st*ra/g*(sl-1.0_dp))*tl
    t0 = ts + gamma_st*zs

    where (abs(zs)<zsmin)
      slp = ps
    elsewhere
      where (t0>tmax)
        where (ts<=tmax)
          gamma = (tmax-ts)/zs
        elsewhere
          gamma = 0.0_dp
          ts = 0.5_dp*(tmax + ts)
        end where
      elsewhere
        where (ts<tmin)
          ts = 0.5_dp*(tmin + ts)
        end where
      end where
      x = gamma*zs/ts
      slp = ps * exp(g*zs/(ra*ts)*(1.0_dp-(0.5_dp - x/3.0_dp)*x))
    end where

  end subroutine slp_calc
  
end module slp_module
