program test_glatwgtq
  use kind_module, only: i4b
  use glatwgtq_module, only: glatwgtq_test
  implicit none

  integer(kind=i4b), parameter :: un=91
  integer(kind=i4b) :: nlat
  character(len=16) :: nlatstr
  character(len=256) :: fname
  
  print *, "Enter # of latitudes"
  read *, nlat
  print *, "nlat = ", nlat
  write(unit=nlatstr,fmt=*) nlat
  fname="glatwgtq_J"//trim(adjustl(nlatstr))//".txt"
  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
  call glatwgtq_test(nlat,un)
  close(unit=un)

end program test_glatwgtq
