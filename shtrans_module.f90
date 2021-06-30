module shtrans_module
  use kind_module, only : i4b, dp
  implicit none
  private

  integer(kind=i4b), parameter :: isym = 0

  logical, public :: shtrans_verbose, shtrans_isinit = .false.
  real(kind=dp), dimension(:), allocatable, public :: &
    shtrans_wshags, shtrans_wshsgs

  real(kind=dp), private :: kfil
  integer(kind=i4b), private :: nlon, nlat, mdab, ndab, idg, jdg
  real(kind=dp), dimension(:), allocatable, private :: work
  real(kind=dp), dimension(:,:,:), allocatable, private :: g
  logical, private :: lyrev

  public :: shtrans_init, shtrans_clean, shtrans_analysis, shtrans_synthesis, &
    shtrans_transf, shtrans_transb, shtrans_hoskinsfilter, shtrans_truncate, shtrans_test

  private :: checkerror, checkval

contains

  subroutine shtrans_init(nx,ny,nt,np2sp)
    implicit none

    integer(kind=i4b), intent(in) :: nx,ny,nt
    logical, optional, intent(in) :: np2sp

    integer(kind=i4b) :: l1, l2, ierror, lwork, ldwork, lshags, lshsgs
    real(kind=dp), dimension(:), allocatable :: dwork

    if (shtrans_isinit) then
      print *, "shtrans is already intialized"
      print *, "nlon =", nlon, " nlat=", nlat, " nt=", size(g,3)
      return
    end if

    nlon = nx
    nlat = ny
    idg = nlat
    jdg = nlon

    lyrev = .true.
    if (present(np2sp).and.np2sp) then
      lyrev = .false.
    end if
! init analysis
    l1 = min0(nlat,(nlon+2)/2)
    l2 = nlat/2
    lshags  = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    lwork = 4*nlat*(nlat+2)+2
    ldwork = nlat*(nlat+4)
    allocate(shtrans_wshags(lshags),work(lwork),dwork(ldwork))

    call shagsi(nlat,nlon,shtrans_wshags,lshags,work,lwork,dwork,ldwork,ierror)
    call checkerror(ierror,"shagsi")
    deallocate(work,dwork)

! init syntheis
    lshsgs = nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15
    lwork = 4*nlat*(nlat+2)+2
    ldwork = nlat*(nlat+4)
    allocate(shtrans_wshsgs(lshsgs),work(lwork),dwork(ldwork))
    call shsgsi(nlat,nlon,shtrans_wshsgs,lshsgs,work,lwork,dwork,ldwork,ierror)
    call checkerror(ierror,"shsgsi")
    deallocate(work,dwork)

    allocate(g(nlat,nlon,nt))

    lwork = nlat*nlon*(nt+1)
    allocate(work(lwork))

    mdab = min0((nlon+2)/2,nlat)
    ndab = nlat

    shtrans_isinit = .true.

  end subroutine shtrans_init

  subroutine shtrans_clean()
    implicit none

    deallocate(work,g,shtrans_wshsgs,shtrans_wshags)
    shtrans_isinit = .false.

  end subroutine shtrans_clean

  subroutine shtrans_analysis(gin,a,b)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: gin
    real(kind=dp), dimension(:,:,:), intent(inout) :: a, b

    integer(kind=i4b) :: ny, nt, j, jr, k, lwork, ierror

    ny = size(gin,2)
    nt = size(gin,3)

    call shtrans_transf(gin,g,lyrev)
    call shags(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, &
               shtrans_wshags,size(shtrans_wshags), &
               work,size(work),ierror)
    call checkerror(ierror,"shags")

    if (shtrans_verbose) then
      call checkval(gin,"gin")
      call checkval(a,"a")
      call checkval(b,"b")
    end if

  end subroutine shtrans_analysis

  subroutine shtrans_synthesis(a,b,gout)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: a, b
    real(kind=dp), dimension(:,:,:), intent(inout) :: gout

    integer(kind=i4b) :: nt, ierror

    nt = size(a,3)

    call shsgs(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, &
               shtrans_wshsgs,size(shtrans_wshsgs), &
               work,size(work),ierror)
    call checkerror(ierror,"shsgs")
    call shtrans_transb(g,gout,lyrev)

    if (shtrans_verbose) then
      call checkval(a,"a")
      call checkval(b,"b")
      call checkval(gout,"gout")
    end if

  end subroutine shtrans_synthesis

  subroutine shtrans_transf(gin,gout,lyrev)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: gin
    real(kind=dp), dimension(:,:,:), intent(inout) :: gout
    logical, optional, intent(in) :: lyrev

    integer(kind=i4b) :: j, jr, k, nx, ny, nt

    nx = size(gin,1)
    ny = size(gin,2)
    nt = size(gin,3)
    do k=1, nt
      if (lyrev) then
        do j=1, ny
          jr = ny-j+1
          gout(jr,1:nx,k) = gin(1:nx,j,k)
        end do
      else
        do j=1, ny
          gout(j,1:nx,k) = gin(1:nx,j,k)
        end do
      end if
    end do
    
  end subroutine shtrans_transf

  subroutine shtrans_transb(gin,gout,lyrev)
    implicit none   

    real(kind=dp), dimension(:,:,:), intent(in) :: gin
    real(kind=dp), dimension(:,:,:), intent(inout) :: gout
    logical, optional, intent(in) :: lyrev

    integer(kind=i4b) :: j, jr, k, nx, ny, nt

    nx = size(gout,1)
    ny = size(gout,2)
    nt = size(gout,3)
    do k=1, nt
      if (lyrev) then
        do j=1, ny
          jr = ny-j+1
          gout(1:nx,jr,k) = gin(j,1:nx,k)
        end do
      else
        do j=1, ny
          gout(1:nx,j,k) = gin(j,1:nx,k)
        end do
      end if
    end do

  end subroutine shtrans_transb

  subroutine shtrans_truncate(a,ntrunc)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(inout) :: a
    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: m, n

    m = size(a,1)
    n = size(a,2)
    a(1:m,ntrunc+2:n,:) = 0.0d0

  end subroutine shtrans_truncate

  subroutine shtrans_hoskinsfilter(a,ntrunc)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(inout) :: a
    integer(kind=i4b), intent(in) :: ntrunc

    integer(kind=i4b) :: n, m
    real(kind=dp) :: fil, kfil

    kfil = -log(0.1_dp)/(ntrunc*(ntrunc+1))**2
    do n=0, ntrunc
      fil = exp(-kfil*(n*n*(n+1)*(n+1)))
      do m=0, n
        a(m+1,n+1,:) = fil*a(m+1,n+1,:)
      end do
    end do

  end subroutine shtrans_hoskinsfilter

  subroutine checkerror(ierror,name)
    implicit none

    integer(kind=i4b), intent(in) :: ierror
    character(len=*), intent(in) :: name

    if (ierror/=0) then
      print *, "*** error in ",name,": ", ierror
      stop
    end if

  end subroutine checkerror

  subroutine checkval(x,name)
    implicit none

    real(kind=dp), dimension(:,:,:) :: x
    character(len=*), intent(in) :: name

    print *, name, " min=", minval(x), " at ", minloc(x), &
             name, " max=", maxval(x), " at ", maxloc(x)

  end subroutine checkval

  subroutine shtrans_test(nx, ny, nt)
    implicit none

    integer(kind=i4b), intent(in) :: nx, ny, nt
    real(kind=dp), dimension(nx, ny, 1) :: g
    real(kind=dp), dimension(ny, ny, 1) :: a, b

    shtrans_verbose = .true.

    call shtrans_init(nx, ny, nt)

! globally uniform value
    g = 1.0_dp
    call shtrans_analysis(g, a, b)
    call shtrans_synthesis(a, b, g)

    call shtrans_clean()

  end subroutine shtrans_test

end module shtrans_module
