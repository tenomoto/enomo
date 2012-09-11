module grmsm_module
  use kind_module, only: i4b, sp, dp
  implicit none
  private

  type grmsm_file_type
    character(len=8), dimension(4) :: &
      label = (/"ncep    ","rsm     ","mpi     ","version "/)
    real(kind=sp) :: fhour = 0.0
    integer(kind=i4b), dimension(4) :: idate
    integer(kind=i4b) :: nx, ny, nz, nt=3, levmax=100
    real(kind=sp), dimension(:), allocatable :: si, sl, ext
    real(kind=sp), dimension(:,:), allocatable :: &
      gz, q,           & ! surface geopotential and ln(p)
      fm2, fm2x, fm2y, & ! (map factor)** and its derivatives
      flat, flon         ! latitudes, longitutes
    real(kind=sp), dimension(:,:,:), allocatable :: &
      te, uu, vv,      & ! hydrostatic virtual temperature and winds
      pn, tn, wn         ! temperature, pressure and vertical velocity
    real(kind=sp), dimension(:,:,:,:), allocatable :: rq
    logical :: isinit = .false., &
      lnh = .false.,   & ! non-hydrostatic
      lext = .false.
  end type

  public :: grmsm_file_type, grmsm_file_init, grmsm_file_clean, grmsm_set_sigma, &
    grmsm_read, grmsm_set_ext, grmsm_write

contains

  function grmsm_file_init(nx,ny,nz,nt,lnh,lext) result(f)
    integer(kind=i4b) :: nx, ny, nz
    integer(kind=i4b), optional :: nt
    logical, optional :: lnh, lext

    type(grmsm_file_type) :: f

    f % nx = nx
    f % ny = ny
    f % nz = nz
    if (present(nt)) then
      f % nt = nt
    end if
    if (present(lnh)) then
      f % lnh = lnh
    end if
    if (present(lext)) then
      f % lext = lext
    end if

    allocate( &
      f % si(nz+1),f % sl(nz), &
      f % gz(nx,ny),f % q(nx,ny), &
      f % fm2(nx,ny),f % fm2x(nx,ny),f % fm2y(nx,ny), &
      f % flat(nx,ny),f % flon(nx,ny), &
      f % te(nx,ny,nz),f % uu(nx,ny,nz),f % vv(nx,ny,nz), &
      f % rq(nx,ny,nz,f % nt) )
    if (f % lnh) then
      allocate(f % pn(nx,ny,nz),f % tn(nx,ny,nz),f % wn(nx,ny,nz))
    end if
    if (f % lext) then
      allocate(f % ext(512-(6+2*f % levmax)))
      f % ext(:) = 0.0
    end if

    f % isinit = .true.

  end function grmsm_file_init

  subroutine grmsm_file_clean(f)
    type(grmsm_file_type), intent(inout) :: f

    deallocate( &
      f % si,f % sl,f % gz,f % q, &
      f % fm2,f % fm2x,f % fm2y,f % flat,f % flon, &
      f % te,f % uu,f % vv,f % rq )
    if (f % lnh) then
      deallocate(f % tn,f % pn,f % wn)
    end if
    if (f % lext) then
      deallocate(f % ext)
    end if
    
    f % isinit = .false.

  end subroutine grmsm_file_clean

  subroutine grmsm_set_sigma(dsi,f)
    use sigma_module, only: sigma_half2full

    real(kind=dp), dimension(:), intent(in) :: dsi
    type(grmsm_file_type), intent(inout) :: f

    real(kind=dp), dimension(size(dsi)-1) :: dsl

    call sigma_half2full(dsi,dsl)

    if (.not.f % isinit) then
      print *, "call grmsm_file_init first"
      stop
    end if
    f % si(:) = real(dsi(:), kind=sp)
    f % sl(:) = real(dsl(:), kind=sp)

  end subroutine grmsm_set_sigma

  subroutine grmsm_read(un,f)
    integer(kind=i4b), intent(in) :: un
    type(grmsm_file_type), intent(inout) :: f

    real(kind=sp), dimension(2*f % levmax-(f % nz+1)-f % nz) :: dummy

    integer(kind=i4b) :: k

    if (.not.f % isinit) then
      print *, "call grmsm_file_init first"
      stop
    end if

    rewind(un)
    read(un) f % label
    if (f % lext) then
      read(un) f % fhour,f % idate,f % si, f % sl, &
        dummy, f % ext
    else
      read(un) f % fhour,f % idate,f % si, f % sl
    end if
    read(un) f % gz
    read(un) f % q
    do k=1, f % nz
      read(un) f % te(:,:,k)
    end do
    do k=1, f % nz
      read(un) f % uu(:,:,k)
      read(un) f % vv(:,:,k)
    end do
    if (f % lnh) then
      do k=1, f % nz
        read(un) f % pn(:,:,k)
      end do
      do k=1, f % nz
        read(un) f % tn(:,:,k)
      end do
      do k=1, f % nz
        read(un) f % wn(:,:,k)
      end do
    end if
    read(un) f % fm2
    read(un) f % fm2x
    read(un) f % fm2y
    read(un) f % flat
    read(un) f % flon

  end subroutine grmsm_read

  subroutine grmsm_set_ext( &
    rproj,rtruth,rorient,rcenlat,rcenlon,rlftgrd,rbtmgrd,rdelx,rdely,f)
    real(kind=sp) rproj,rtruth,rorient,rcenlat,rcenlon,rlftgrd,rbtmgrd,rdelx,rdely
    type(grmsm_file_type), intent(inout) :: f

    real(kind=sp) :: iwav1, jwav1, igrd1, jgrd1

    if (.not.f % isinit) then
      print *, "call grmsm_file_init first"
      stop
    end if
    if (.not.f % lext) then
      print *, "file is not initialized for ext"
      stop
    end if

    igrd1 = f % nx + 1
    jgrd1 = f % ny + 1
    iwav1 = (f % nx-12)/3*2
    jwav1 = iwav1*f % ny/(f % nx*2)*2 + 1
    iwav1 = iwav1 + 1

    f % ext(1) = iwav1
    f % ext(2) = jwav1
    f % ext(3) = igrd1
    f % ext(4) = jgrd1
    f % ext(5) = f % nz
    f % ext(6) = 2+f % nz*4+5
    f % ext(7) = rproj
    f % ext(8) = rtruth
    f % ext(9) = rorient
    f % ext(10) = rcenlat
    f % ext(11) = rcenlon
    f % ext(12) = rlftgrd
    f % ext(13) = rbtmgrd
    f % ext(14) = rdelx
    f %  ext(15) = rdely
    if (f % lnh) then
      f % ext(16) = 1.0
    end if

  end subroutine grmsm_set_ext

  subroutine grmsm_write(un,f)
    integer(kind=i4b), intent(in) :: un
    type(grmsm_file_type), intent(inout) :: f

    real(kind=sp), dimension(2*f % levmax-(f % nz+1)-f % nz) :: dummy
    logical :: le

    integer(kind=i4b) :: k
    real(kind=sp) :: iwav1, jwav1, igrd1, jgrd1, &
      rtruth, rcenlat, rcenlon, rlftgrd, rbtmgrd, rdelx, rdely

    if (.not.f % isinit) then
      print *, "call grmsm_file_init first"
      stop
    end if

    rewind(un)
    write(un) f % label
    if (f % lext) then
      dummy(:) = 0.0
      write(un) f % fhour,f % idate,f % si, f % sl, &
        dummy, f % ext
    else
      write(un) f % fhour,f % idate,f % si, f % sl
    end if
    write(un) f % gz
    write(un) f % q
    do k=1, f % nz
      write(un) f % te(:,:,k)
    end do
    do k=1, f % nz
      write(un) f % uu(:,:,k)
      write(un) f % vv(:,:,k)
    end do
    if (f % lnh) then
      do k=1, f % nz
        write(un) f % pn(:,:,k)
      end do
      do k=1, f % nz
        write(un) f % tn(:,:,k)
      end do
      do k=1, f % nz
        write(un) f % wn(:,:,k)
      end do
    end if
    write(un) f % fm2
    write(un) f % fm2x
    write(un) f % fm2y
    write(un) f % flat
    write(un) f % flon

  end subroutine grmsm_write

end module grmsm_module
