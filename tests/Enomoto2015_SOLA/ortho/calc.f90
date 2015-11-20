program calc
!$ use omp_lib
  use kind_module, only: i4b, dp, qp
  use glatwgt_module, only: glatwgt_calc
  use glatwgtq_module, only: glatwgtq_calc
!  use alfx_module, only: alfx_init, alfx_clean, alfx_calc
  use alfx_module, only: alfx_init, alfx_clean, alfx_calc_inline
!  use alff_module, only: alff_init, alff_clean, alff_calc
!  use alfq_module, only: alfq_init, alfq_clean, alfq_calc
  use alfxq_module, only: alfxq_init, alfxq_clean, alfxq_calc_inline
  implicit none
  
  integer :: iargc

  integer(kind=i4b), parameter :: unx = 41, unf = 42, unq = 43
  character(len=*), parameter :: datadir = "/Volumes/Pegasus/alfx/run/"

  integer(kind=i4b) :: ntrunc, nlat, nlath, j
  character(len=5) :: ntrunc_str
  real(kind=dp) :: lat(1)
  real(kind=dp), dimension(:), allocatable :: glat, gwgt
  real(kind=dp), dimension(:, :), allocatable :: s, c
  real(kind=dp), dimension(:,:,:), allocatable :: pnm
  real(kind=qp) :: latq(1)
  real(kind=qp), dimension(:), allocatable :: glatq, gwgtq
  real(kind=qp), dimension(:,:,:), allocatable :: pnmq

  if (iargc() < 1) then
    print *, "Usage :: calc ntrunc"
    stop
  end if

  call getarg(1,ntrunc_str)
  read(unit=ntrunc_str,fmt=*) ntrunc
  nlat = (ntrunc+1)*3/2
  nlath = nlat/2
!  print *, "ntrunc=", ntrunc, " nlat=", nlat, " jmax=", jmax

  allocate(pnm(0:ntrunc,0:ntrunc,1), pnmq(0:ntrunc,0:ntrunc,1))

!  open(unit=unx, file=datadir//"alfx_T"//trim(ntrunc_str)//".dat", &
  open(unit=unx, file=datadir//"alfxi_T"//trim(ntrunc_str)//".dat", &
    access="direct", recl=size(pnm)*8, &
    status="replace",action="write")
!  open(unit=unf, file=datadir//"alff_T"//trim(ntrunc_str)//".dat", &
!    access="direct", recl=size(pnm)*8, &
!    status="replace",action="write")
!  open(unit=unq, file=datadir//"alfq_T"//trim(ntrunc_str)//".dat", &
  open(unit=unq, file=datadir//"alfxq_T"//trim(ntrunc_str)//".dat", &
    access="direct", recl=size(pnmq)*16, &
    status="replace",action="write")

  allocate(glat(nlat), gwgt(nlat))
  call glatwgt_calc(glat,gwgt)
  allocate(glatq(nlat), gwgtq(nlat))
  call glatwgtq_calc(glatq,gwgtq)

  call alfx_init(ntrunc)
!  call alff_init(ntrunc)
!  call alfq_init(ntrunc)
  call alfxq_init(ntrunc)

!$omp parallel do private(lat, latq, pnm, pnmq)
  do j = 1, nlath
    lat(1) = glat(j)
    latq(1) = glatq(j)
!    call alfx_calc(lat,pnm)
    call alfx_calc_inline(lat,pnm)
    write(unit=unx, rec=j) pnm
!    call alff_calc(lat,pnm)
!    write(unit=unf, rec=j) pnm
!    call alfq_calc(latq,pnmq)
    call alfxq_calc_inline(latq,pnmq)
    write(unit=unq, rec=j) pnmq
  end do
!$omp end parallel do

  close(unit=unx)
!  close(unit=unf)
  close(unit=unq)

  deallocate(glat, gwgt, glatq, gwgtq)
  deallocate(pnm, pnmq)
  call alfx_clean()
!  call alff_clean()
!  call alfq_clean()
  call alfxq_clean()

end program calc
