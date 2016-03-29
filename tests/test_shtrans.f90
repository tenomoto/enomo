program test_shtrans
  use kind_module, only: i4b
  use shtrans_module, only: shtrans_test
  implicit none

  integer(kind=i4b) :: nlat, nlon, ntrunc
  character(len=16) :: nlatstr
  
  print *, "Enter # of latitudes"
  read *, nlat
  nlon = nlat * 2
  ntrunc = (nlon - 1) / 3
  print *, "nlat = ", nlon, " nlat = ", nlat, " ntrunc = ", ntrunc
  call shtrans_test(nlon, nlat, ntrunc)

end program test_shtrans
