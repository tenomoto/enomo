module shregrid_module
  use kind_module, only: i4b, dp
  implicit none
  private

! a wrapper for a spherepack3.x subroutine trssph
! grid transfer between Gaussian to equally spaced

  integer(kind=i4b), parameter, public :: &
    shregrid_f2g = 0, shregrid_g2g = 1, shregrid_g2f = 3, shregrid_f2f = 4

  integer(kind=i4b), parameter, private :: &
    iveca = 0, & ! va is a latitudinal component
    ivecb = 0    ! vb is a latitudinal component
  integer(kind=i4b), dimension(2), parameter, private :: &
    igride = (/ 1, 0/), & ! SP->NP, equally spaced, IxJ 
    igridg = (/ 2, 0/)    ! SP->NP, Gaussian, IxJ
  integer(kind=i4b), private :: intl, lsave, lsavev, lwork, lworkv, ldwork, ldworkv, &
    nlona, nlata, nlonb, nlatb, lsvmin, lwkmin, ier, nlat
  real(kind=dp), dimension(:), allocatable, private :: wsave, wsavev
  real(kind=dp), dimension(:), allocatable, private :: work, workv
  real(kind=dp), dimension(:), allocatable, private :: dwork, dworkv
  integer(kind=i4b), dimension(2), private :: igrida, igridb

! dummy arrays
  real(kind=dp), dimension(1), private :: wsave_dummy, work_dummy

  public :: shregrid_init, shregrid_clean, shregrid_regrid
  private :: regrid_scalar, regrid_vector
  interface shregrid_regrid
    module procedure regrid_scalar, regrid_vector
  end interface

contains

  subroutine shregrid_init(ua, va, ub, vb, regrid_type)
    implicit none
    real(kind=dp), dimension(:,:), intent(inout) :: ua, va, ub, vb
    integer(kind=i4b), intent(in) :: regrid_type

    select case(regrid_type)
      case (shregrid_f2g)
        igrida = igride
        igridb = igridg
      case (shregrid_g2g)
        igrida = igridg
        igridb = igridg
      case (shregrid_g2f)
        igrida = igridg
        igridb = igride
      case (shregrid_f2f)
        igrida = igride
        igridb = igride
    end select

! determine the size of data arrays
    nlona = size(ua, 1)
    nlata = size(ua, 2)
    nlonb = size(ub, 1)
    nlatb = size(ub, 2)
    nlat = max(nlata, nlatb)

! calculate size of working arrays and allocate them
    intl = 0

    ua = 0.
    va = 0.
    ub = 0.
    vb = 0
    lsavev = 0
    lworkv = 0
    ldworkv = 2*nlat*(nlat+1)+1
    allocate(dworkv(ldworkv))
    call trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
      igridb, nlonb, nlatb, ivecb, ub, vb, wsave_dummy, lsavev, lsvmin, &
      work_dummy, lworkv, lwkmin, dworkv, ldworkv, ier)
    lsavev = lsvmin
    lworkv = lwkmin
    allocate(wsavev(lsavev))
    allocate(workv(lworkv))
    call trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
      igridb, nlonb, nlatb, ivecb, ub, vb, wsavev, lsavev, lsvmin, &
      workv, lworkv, lwkmin, dworkv, ldworkv, ier)
    print *, "shregrid_init vector: ier=", ier, " size(wsavev)=", size(wsavev)

    lsave = 0
    lwork = 0
    ldwork = nlat*(nlat+4)
    allocate(dwork(ldwork))
    call trssph(intl, igrida, nlona, nlata, ua, &
      igridb, nlonb, nlatb, ub, wsave_dummy, lsave, lsvmin, &
      work_dummy, lwork, lwkmin, dwork, ldwork, ier)
    lsave = lsvmin
    lwork = lwkmin
    allocate(wsave(lsave))
    allocate(work(lwork))
    ua = 0.
    ub = 0.
    call trssph(intl, igrida, nlona, nlata, ua, &
      igridb, nlonb, nlatb, ub, wsave, lsave, lsvmin, &
      work, lwork, lwkmin, dwork, ldwork, ier)
    print *, "shregrid_init scalar: ier=", ier, " size(wsave)=", size(wsave)

    intl = 1

  end subroutine shregrid_init

  subroutine shregrid_clean
    implicit none

    deallocate(wsave, wsavev, work, workv, dwork, dworkv)

  end subroutine shregrid_clean

  subroutine regrid_scalar(da,db)
    implicit none

    real(kind=dp), dimension(nlona,nlata), intent(inout) :: da
    real(kind=dp), dimension(nlonb,nlatb), intent(inout) :: db

!print *, "intl=", intl, &
!" nlona=", nlona, " nlata=", nlata, " size(da)=", size(da), &
!" nlonb=", nlonb, " nlatb=", nlatb, " size(db)=", size(db), &
!" lsave=", lsave , " lwork=", lwork, " ldwork=", ldwork
!print *, maxval(wsave), minval(wsave)
!print *, maxval(da), minval(da)
    work = 0.
    dwork = 0.
    call trssph(intl, igrida, nlona, nlata, da, &
      igridb, nlonb, nlatb, db, wsave, lsave, lsvmin, &
      work, lwork, lwkmin, dwork, ldwork, ier)
    if (ier/=0) then
      print *, "shregrid_regrid (scalar) ier=", ier
      stop
    end if
!print *, maxval(db), minval(db)
!stop

  end subroutine regrid_scalar

  subroutine regrid_vector(ua, va, ub, vb)
    implicit none

    real(kind=dp), dimension(nlona,nlata), intent(inout) :: ua, va
    real(kind=dp), dimension(nlonb,nlatb), intent(inout) :: ub, vb

    workv = 0.
    dworkv = 0.
    call trvsph(intl, igrida, nlona, nlata, iveca, ua, va, &
      igridb, nlonb, nlatb, ivecb, ub, vb, wsavev, lsavev, lsvmin, &
      workv, lworkv, lwkmin, dworkv, ldworkv, ier)
    if (ier/=0) then
      print *, "shregrid_regrid (vector) ier=", ier
      stop
    end if

  end subroutine regrid_vector

end module shregrid_module
