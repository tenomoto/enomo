program test_alff
  use kind_module, only: i4b
  use alff_module, only: alff_test, alff_test_checksum
  implicit none

  integer(kind=i4b), parameter :: un=91
  integer(kind=i4b) :: ntrunc, nlat
  character(len=16) :: ntruncstr, nlatstr
  character(len=256) :: fname
  
  print *, "Enter truncation wave number:"
  read *, ntrunc
  nlat = 3 * (ntrunc + 1) / 2
  write(unit=ntruncstr,fmt=*) ntrunc
  write(unit=nlatstr,fmt=*) nlat

  fname="alff_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".dat"
  open(unit=un,file=trim(adjustl(fname)), &
    access="direct", status="replace", action="write", &
    recl=(ntrunc+1)*(ntrunc+1)*nlat/2*8)
  call alff_test(ntrunc,nlat)
  close(unit=un)

  fname="alff_checksum_T"//trim(adjustl(ntruncstr))//"J"//trim(adjustl(nlatstr))//".txt"
  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
  call alff_test_checksum(ntrunc,nlat,un)
  close(unit=un)

end program test_alff
