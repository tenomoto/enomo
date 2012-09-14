module grmsm_sfc_module
  use kind_module, only: i4b, sp, dp
  implicit none
  private

  integer(kind=i4b), parameter, private :: lsoil0 = 2
  type grmsm_sfc_type
    character(len=8), dimension(4) :: &
      label = (/"ncep    ","rsm     ","mpi     ","version "/)
    real(kind=sp) :: fhour = 0.0
    integer(kind=i4b), dimension(4) :: idate
    integer(kind=i4b) :: nx, ny, lsoil = lsoil0
    real(kind=i4b), dimension(:,:), allocatable :: &
      tsea,   & ! land and sea surface temperature (k)
      sheleg, & ! snow depth (cm)
      tg3,    & ! the lowest soil temperature at 3m (K)
      zorl,   & ! surface roughness
      cv,     & ! cloud amount
      cvb,    & ! cloud base (sigma layer number)
      cvt,    & ! cloud top (sigma layer number)
      alvsf,  & ! albedo
      alvwf,  &
      alnsf,  &
      alnwf,  &
      slmsk,  & ! sea land mask
      vfrac,  & !
      canopy, & ! surface canopy
      f10m,   & ! 10 m height factor for wind
      vtype,  & !
      stype,  & !
      facsf,  & !
      facwf,  & !
      uustar, & !
      ffmm,   &
      ffhh
    real(kind=i4b), dimension(:,:,:), allocatable :: &
      smc,     & ! two layers of soil moisture contents (0.47 - 0.1)
      stc        ! two layers of soil temperature (k)
    logical :: isinit = .false.
  end type grmsm_sfc_type

  public :: grmsm_sfc_type, grmsm_sfc_init, grmsm_sfc_clean, grmsm_sfc_read, grmsm_sfc_write

contains

  function grmsm_sfc_init(nx,ny,lsoil) result(f)
    integer(kind=i4b), intent(in) :: nx, ny
    integer(kind=i4b), intent(in), optional :: lsoil

    type(grmsm_sfc_type) :: f

    if (present(lsoil)) then
      f % lsoil = lsoil
    end if
    f % nx = nx
    f % ny = ny

    allocate( &
      f % tsea(nx,ny),   f % sheleg(nx,ny), f % tg3(nx,ny),   &
      f % zorl(nx,ny),   f % cv(nx,ny),     f % cvb(nx,ny),   &
      f % cvt(nx,ny),    f % alvsf(nx,ny),  f % alvwf(nx,ny), &
      f % alnsf(nx,ny),  f % alnwf(nx,ny),  f % slmsk(nx,ny), &
      f % vfrac(nx,ny),  f % canopy(nx,ny), f % f10m(nx,ny),  &
      f % vtype(nx,ny),  f % stype(nx,ny),  f % facsf(nx,ny), &
      f % facwf(nx,ny),  f % uustar(nx,ny), f % ffmm(nx,ny) )

    allocate(f % smc(nx,ny,f % lsoil),f % stc(nx,ny,f% lsoil))

    f % isinit = .true. 

  end function grmsm_sfc_init

  subroutine grmsm_sfc_clean(f)
    type(grmsm_sfc_type), intent(inout) :: f

    deallocate ( &
      f % tsea,   f % sheleg, f % tg3, &
      f % zorl,   f % cv,     f % cvb, &
      f % cvt,    f % alvsf,  f % alvwf, &
      f % alnsf,  f % alnwf,  f % slmsk, &
      f % vfrac,  f % canopy, f % f10m, &
      f % vtype,  f % stype,  f % facsf, &
      f % facwf,  f % uustar, f % ffmm )

    deallocate(f % smc,f % stc)

    f % lsoil = lsoil0
    f % isinit = .false.

  end subroutine grmsm_sfc_clean

  subroutine grmsm_sfc_read(un,f)
    integer(kind=i4b), intent(in) :: un
    type(grmsm_sfc_type), intent(inout) :: f

    if (.not.f % isinit) then
      print *, "call grmsm_sfc_init first"
      stop
    end if

    rewind(un)
    read(un) f % label
! rfixio writes igrd1, jgrd1 and version
! in addition to fhour and idate
    read(un) f % fhour,f % idate
    read(un) f % tsea
    read(un) f % smc
    read(un) f % sheleg
    read(un) f % stc
    read(un) f % tg3
    read(un) f % zorl
    read(un) f % cv
    read(un) f % cvb
    read(un) f % cvt
    read(un) f % alvsf, f % alvwf, f % alnsf, f % alnwf
    read(un) f % slmsk
    read(un) f % vfrac
    read(un) f % canopy
    read(un) f % f10m
    read(un) f % vtype
    read(un) f % stype
    read(un) f % facsf, f % facwf
    read(un) f % uustar
    read(un) f % ffmm
    read(un) f % ffhh
    
  end subroutine grmsm_sfc_read

  subroutine grmsm_sfc_write(un,f)
    integer(kind=i4b), intent(in) :: un
    type(grmsm_sfc_type), intent(inout) :: f

    if (.not.f % isinit) then
      print *, "call grmsm_sfc_init first"
      stop
    end if

    rewind(un)
    write(un) f % label
! rfixio writes igrd1, jgrd1 and version
! in addition to fhour and idate
    write(un) f % fhour,f % idate
    write(un) f % tsea
    write(un) f % smc
    write(un) f % sheleg
    write(un) f % stc
    write(un) f % tg3
    write(un) f % zorl
    write(un) f % cv
    write(un) f % cvb
    write(un) f % cvt
    write(un) f % alvsf, f % alvwf, f % alnsf, f % alnwf
    write(un) f % slmsk
    write(un) f % vfrac
    write(un) f % canopy
    write(un) f % f10m
    write(un) f % vtype
    write(un) f % stype
    write(un) f % facsf, f % facwf
    write(un) f % uustar
    write(un) f % ffmm
    write(un) f % ffhh
    
  end subroutine grmsm_sfc_write

end module grmsm_sfc_module
