program test_gaqd
  use kind_module, only: i4b, dp
  implicit none

  integer(kind=i4b), parameter :: un=91
  integer(kind=i4b) :: nlat, ierror, lwork
  character(len=16) :: nlatstr
  character(len=256) :: fname
  real(kind=dp), dimension(:), allocatable :: theta, wts
  real(kind=dp) :: w, s 
  
  print *, "Enter # of latitudes"
  read *, nlat
  print *, "nlat = ", nlat
  allocate(theta(nlat), wts(nlat))
  write(unit=nlatstr,fmt=*) nlat
  fname="gaqd_J"//trim(adjustl(nlatstr))//".txt"
  call gaqd(nlat,theta,wts,w,lwork,ierror)
  open(unit=un,file=trim(adjustl(fname)), status="replace", action="write")
  write(unit=un,fmt=*) asin(cos(theta))
  write(unit=un,fmt=*) wts
  close(unit=un)
  s = sum(wts(1:nlat/2))
  print *, "sum of weights:", s, " error=", abs(1.0_dp - s)
  deallocate(theta, wts)

end program test_gaqd
