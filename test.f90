program test
  use kind_module, only: i4b
  use xreal_module, only: xreal_test
!  use glatwgt_module, only: glatwgt_test
!  use alf_module, only: alf_test, alf_test_checksum, alf_init, alf_clean
!  use alfx_module, only: alfx_test, alfx_test_checksum, alfx_init, alfx_clean
!  use alff_module, only: alff_test, alff_test_checksum, alff_init, alff_clean
  implicit none

  integer(kind=i4b), parameter :: un=91
  integer(kind=i4b) :: ntrunc, nlat
  character(len=8) :: ntruncstr, nlatstr
  character(len=256) :: fname

  call xreal_test()

!  ntrunc = 39
!  nlat = 60
!  write(unit=ntruncstr,fmt=*) ntrunc
!  write(unit=nlatstr,fmt=*) nlat

!  fname="glatwgt_J"//trim(adjustl(nlatstr))//".txt"
!  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
!  call glatwgt_test(nlat,un)
!  close(unit=un)

!  fname="alf_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".dat"
!  open(unit=un,file=trim(adjustl(fname)), &
!    access="direct", status="replace", action="write", &
!    recl=(ntrunc+1)*(ntrunc+1)*nlat/2*8)
!  call alf_test(ntrunc,nlat,un)
!!  call alf_test(ntrunc,nlat)
!  close(unit=un)

!  fname="alf_checksum_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".txt"
!  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
!  call alf_test_checksum(ntrunc,nlat,un)
!  close(unit=un)
!
!  fname="alfx_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".dat"
!  open(unit=un,file=trim(adjustl(fname)), &
!    access="direct", status="replace", action="write", &
!    recl=(ntrunc+1)*(ntrunc+1)*nlat/2*8)
!  call alfx_test(ntrunc,nlat,un)
!!  call alfx_test(ntrunc,nlat)
!  close(unit=un)

!  fname="alfx_checksum_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".txt"
!  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
!  call alfx_test_checksum(ntrunc,nlat,un)
!  close(unit=un)
!
!  fname="alff_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".dat"
!  open(unit=un,file=trim(adjustl(fname)), &
!    access="direct", status="replace", action="write", &
!    recl=(ntrunc+1)*(ntrunc+1)*nlat/2*8)
!  call alff_test(ntrunc,nlat,un)
!!  call alff_test(ntrunc,nlat)
!  close(unit=un)
!
!  fname="alff_checksum_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".txt"
!  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
!  call alff_test_checksum(ntrunc,nlat,un)
!  close(unit=un)
!
!   call alf_init(ntrunc)
!   call alfx_init(ntrunc)
!   call alff_init(ntrunc)
!   call alff_clean()
!   call alfx_clean()
!   call alf_clean()
end program test
