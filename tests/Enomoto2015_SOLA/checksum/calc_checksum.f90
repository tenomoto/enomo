program calc_checksum
  use kind_module, only: i4b, dp
  use math_module, only: rad2deg=>math_rad2deg
  use glatwgt_module, only: glatwgt_calc
  use alf_module, only: alf_init, alf_calc, alf_clean
  use alff_module, only: alff_init, alff_calc, alff_clean
  use alfx_module, only: alfx_init, alfx_calc, alfx_calc_inline, alfx_clean
  use mpi
  implicit none

  ! mpi
  integer(kind=i4b), parameter :: root = 0
  integer :: ierr, nproc, myrank

  integer(kind=i4b), parameter :: iunit = 41
  character(len=*), parameter :: &
    ifile = "ntrunc.txt"
  integer(kind=i4b) :: ntrunc, jmax, j, jg, nlat, m
  character(len=128) :: suffix
  real(kind=dp), dimension(1) :: lat
  real(kind=dp), dimension(:), allocatable :: latg, wgtg
  real(kind=dp), dimension(:,:), allocatable :: ppw
  real(kind=dp), dimension(:,:,:), allocatable :: alf

  open(unit=iunit,file=ifile,status="old",action="read")
  read(unit=iunit,fmt=*) ntrunc
  close(unit=iunit)
  call init()

  write(unit=suffix,fmt="(a2,i0.5,a1,i0.5,a10,i0.2,a4)") &
     "_T",ntrunc,"J",nlat,"_checksum_",myrank,".dat"

  ! alf
  ppw(:,:) = 0.0_dp
  do j=1, jmax
    jg = jmax*myrank + j
    lat(1) = latg(jg)
    call alf_calc(lat,alf)
    do m=0, ntrunc
      ppw(m:ntrunc,m) = ppw(m:ntrunc,m) + wgtg(jg)*alf(m:ntrunc,m,1)**2
    end do
    if (myrank==0) then
      print "(a,f6.1,a1)", "alf:", real(j)/jmax*100, "%"
    end if
  end do
  call save_data("alf"//suffix,ppw)

  ! alfx
  ppw(:,:) = 0.0_dp
  do j=1, jmax
    jg = jmax*myrank + j
    lat(1) = latg(jg)
    call alfx_calc(lat,alf)
!    call alfx_calc_inline(lat,alf)
    do m=0, ntrunc
      ppw(m:ntrunc,m) = ppw(m:ntrunc,m) + wgtg(jg)*alf(m:ntrunc,m,1)**2
    end do
    if (myrank==0) then
      print "(a,f6.1,a1)", "alfx:", real(j)/jmax*100, "%"
    end if
  end do
  call save_data("alfx"//suffix,ppw)

  ! alff
  ppw(:,:) = 0.0_dp
  do j=1, jmax
    jg = jmax*myrank + j
    lat(1) = latg(jg)
    call alff_calc(lat,alf)
    do m=0, ntrunc
      ppw(m:ntrunc,m) = ppw(m:ntrunc,m) + wgtg(jg)*alf(m:ntrunc,m,1)**2
    end do
    if (myrank==0) then
      print "(a,f6.1,a1)", "alff:", real(j)/jmax*100, "%"
    end if
  end do
  call save_data("alff"//suffix,ppw)

  call clean()

contains

  subroutine init()

    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    nlat = (ntrunc+1)*3/2
    jmax = (nlat/2)/nproc

    allocate(latg(nlat),wgtg(nlat))
    allocate(alf(0:ntrunc,0:ntrunc,1),ppw(0:ntrunc,0:ntrunc))

    call glatwgt_calc(latg, wgtg)
    wgtg(:) = 2.0_dp*wgtg(:)

    call alf_init(ntrunc)
    call alfx_init(ntrunc)
    call alff_init(ntrunc)

  end subroutine init

  subroutine clean()

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call alf_clean()
    call alfx_clean()
    call alff_clean()
    deallocate(latg, wgtg, alf)
    call MPI_Finalize(ierr)

  end subroutine clean

  subroutine save_data(ofile, data)

    character(len=*), intent(in) :: ofile
    real(kind=dp), dimension(:,:), intent(in) :: data

    integer(kind=i4b), parameter :: ounit = 42

    print *, "Saving file to ", trim(adjustl(ofile))
    open(unit=ounit,file=trim(adjustl(ofile)),access="direct", &
      status="replace",action="write",recl=size(data)*8)
    write(unit=ounit,rec=1) data
    close(unit=ounit)

  end subroutine save_data

end program calc_checksum
