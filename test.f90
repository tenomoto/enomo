program test
  use kind_module, only: i4b
  use sum_module, only: sum_test
!  use xreal_module, only: xreal_test
!  use glatwgt_module, only: glatwgt_test
!  use alf_module, only: alf_test, alf_test_checksum
!  use alfx_module, only: alfx_test, alfx_test_checksum
!  use alff_module, only: alff_test, alff_test_checksum
  implicit none

!  integer(kind=i4b), parameter :: un=91
!  integer(kind=i4b) :: ntrunc, nlat
!  character(len=8) :: ntruncstr, nlatstr
!  character(len=256) :: fname
  
  call sum_test(1000000)
!  call xreal_test()

!  nlat = 10239
!  fname="glatwgt_J"//trim(adjustl(nlatstr))//".txt"
!  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
!  call glatwgt_test(nlat,un)
!  close(unit=un)

!  ntrunc = 239
!  nlat = 360
!  write(unit=ntruncstr,fmt=*) ntrunc
!  write(unit=nlatstr,fmt=*) nlat

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
end program test
